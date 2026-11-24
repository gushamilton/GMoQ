// Native, multi-threaded decoder for BGEN layout-2 biallelic diploid records.
//
// The Python reader thread reads each selected record's raw bytes with one pread
// (cheap: ~0.15 ms for a 400k-sample UKB record); this decodes a whole block of
// them in parallel with the GIL released: zlib inflate, then the packed
// probabilities to the allele-1 (first allele) dosage
//
//     dosage = (2 * P(11) + P(12)) = (2 * v0 + v1) / (2^B - 1)
//
// computed exactly in double precision for any bit depth B in 1..32.  The `bgen`
// package returns float32 (2 - alt_dosage carries ~1e-7 rounding); this is the
// exact value of the stored integers.  Samples flagged missing in the ploidy
// byte get NaN, which the Stage 1 kernels mean-impute.
//
// Phased records, ploidy other than 2, and zstd compression are not handled
// here: the record is flagged UNSUPPORTED and the caller falls back to `bgen`.

#ifndef GMOQ_STAGE1_BGEN_DECODE_HPP
#define GMOQ_STAGE1_BGEN_DECODE_HPP

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <zlib.h>
#ifdef GMOQ_HAVE_LIBDEFLATE
#include <libdeflate.h>
#endif
#include <unistd.h>
#include <cerrno>

#include <cstdint>
#include <cstring>
#include <limits>
#include <string>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace gmoq_bgen {

namespace py = pybind11;

enum RecordStatus : int8_t {
    REC_OK = 0,
    REC_NON_BIALLELIC = 1,
    REC_UNSUPPORTED = 2,
    REC_CORRUPT = 3,
};

struct RecordMeta {
    std::string rsid, chrom, allele0, allele1;
    uint32_t pos = 0;
    uint16_t n_alleles = 0;
    int8_t status = REC_CORRUPT;
    std::string error;
};

class Cursor {
public:
    Cursor(const uint8_t* p, size_t n) : p_(p), n_(n) {}
    bool ok() const { return ok_; }
    size_t offset() const { return off_; }
    const uint8_t* here() const { return p_ + off_; }
    template <typename T> T read() {
        T v{};
        if (!need(sizeof(T))) return v;
        std::memcpy(&v, p_ + off_, sizeof(T));
        off_ += sizeof(T);
        return v;
    }
    std::string str16() { return str(read<uint16_t>()); }
    std::string str32() { return str(read<uint32_t>()); }
    bool skip(size_t k) { if (!need(k)) return false; off_ += k; return true; }
private:
    bool need(size_t k) { if (!ok_ || off_ + k > n_) { ok_ = false; return false; } return true; }
    std::string str(size_t k) {
        if (!need(k)) return std::string();
        std::string s(reinterpret_cast<const char*>(p_ + off_), k);
        off_ += k;
        return s;
    }
    const uint8_t* p_;
    size_t n_;
    size_t off_ = 0;
    bool ok_ = true;
};

// zlib-stream inflater: libdeflate when built with it (one decompressor per
// thread), zlib otherwise.  Both are exact; libdeflate is ~2x faster.
struct Inflater {
#ifdef GMOQ_HAVE_LIBDEFLATE
    libdeflate_decompressor* d = libdeflate_alloc_decompressor();
    ~Inflater() { if (d) libdeflate_free_decompressor(d); }
#endif
    // Returns true iff exactly dst_len bytes were produced.
    bool inflate(const uint8_t* src, size_t src_len, uint8_t* dst, size_t dst_len) {
#ifdef GMOQ_HAVE_LIBDEFLATE
        if (d != nullptr) {
            size_t got = 0;
            const libdeflate_result r = libdeflate_zlib_decompress(d, src, src_len, dst, dst_len, &got);
            return r == LIBDEFLATE_SUCCESS && got == dst_len;
        }
#endif
        uLongf dest_len = static_cast<uLongf>(dst_len);
        const int rc = uncompress(dst, &dest_len, src, static_cast<uLong>(src_len));
        return rc == Z_OK && dest_len == dst_len;
    }
};

inline bool built_with_libdeflate() {
#ifdef GMOQ_HAVE_LIBDEFLATE
    return true;
#else
    return false;
#endif
}

// Read B bits (B <= 32) starting at bit `bit` of a buffer padded by >= 8 bytes.
inline uint64_t read_bits(const uint8_t* buf, uint64_t bit, uint32_t B) {
    uint64_t word;
    std::memcpy(&word, buf + (bit >> 3), sizeof(word));
    return (word >> (bit & 7u)) & ((B == 64) ? ~0ull : ((1ull << B) - 1ull));
}

// Decode one record into column `out` (length n_out) for samples `idx`.
inline void decode_record(const uint8_t* rec, size_t rec_len, int compression,
                          uint32_t n_samples_file, const int64_t* idx, int64_t n_out,
                          double* out, std::vector<uint8_t>& scratch, Inflater& inf, RecordMeta& meta) {
    Cursor c(rec, rec_len);
    c.str16();                       // variant id
    meta.rsid = c.str16();
    meta.chrom = c.str16();
    meta.pos = c.read<uint32_t>();
    meta.n_alleles = c.read<uint16_t>();
    std::vector<std::string> alleles;
    for (uint16_t a = 0; a < meta.n_alleles && c.ok(); ++a) alleles.push_back(c.str32());
    if (!c.ok()) { meta.status = REC_CORRUPT; meta.error = "truncated variant header"; return; }
    if (meta.n_alleles != 2) { meta.status = REC_NON_BIALLELIC; return; }
    meta.allele0 = alleles[0];
    meta.allele1 = alleles[1];

    const uint32_t C = c.read<uint32_t>();
    if (!c.ok() || C > rec_len - c.offset()) {
        meta.status = REC_CORRUPT; meta.error = "genotype block runs past the record"; return;
    }
    const uint8_t* data;
    size_t data_len;
    if (compression == 0) {
        scratch.resize(static_cast<size_t>(C) + 8);
        std::memcpy(scratch.data(), c.here(), C);
        data_len = C;
    } else if (compression == 1) {
        if (C < 4) { meta.status = REC_CORRUPT; meta.error = "short zlib block"; return; }
        const uint32_t D = c.read<uint32_t>();
        scratch.resize(static_cast<size_t>(D) + 8);
        if (!inf.inflate(c.here(), C - 4, scratch.data(), D)) {
            meta.status = REC_CORRUPT; meta.error = "zlib inflate failed"; return;
        }
        data_len = D;
    } else {
        meta.status = REC_UNSUPPORTED; meta.error = "compression"; return;
    }
    std::memset(scratch.data() + data_len, 0, 8);   // padding for read_bits
    data = scratch.data();

    Cursor g(data, data_len);
    const uint32_t N = g.read<uint32_t>();
    const uint16_t K = g.read<uint16_t>();
    const uint8_t pmin = g.read<uint8_t>();
    const uint8_t pmax = g.read<uint8_t>();
    if (!g.ok() || N != n_samples_file || K != 2) {
        meta.status = REC_CORRUPT; meta.error = "genotype block header mismatch"; return;
    }
    if (pmin != 2 || pmax != 2) {
        // Ploidy bytes still carry the missing flag; only exactly-diploid is handled.
        meta.status = REC_UNSUPPORTED; meta.error = "ploidy"; return;
    }
    const uint8_t* ploidy = g.here();
    if (!g.skip(N)) { meta.status = REC_CORRUPT; meta.error = "truncated ploidy"; return; }
    const uint8_t phased = g.read<uint8_t>();
    const uint8_t B = g.read<uint8_t>();
    if (!g.ok()) { meta.status = REC_CORRUPT; meta.error = "truncated genotype header"; return; }
    if (phased != 0) { meta.status = REC_UNSUPPORTED; meta.error = "phased"; return; }
    if (B < 1 || B > 32) { meta.status = REC_CORRUPT; meta.error = "bit depth"; return; }
    const uint64_t prob_bits = static_cast<uint64_t>(N) * 2u * B;
    if (g.offset() + (prob_bits + 7) / 8 > data_len) {
        meta.status = REC_CORRUPT; meta.error = "truncated probabilities"; return;
    }
    const uint8_t* probs = g.here();
    const double nan = std::numeric_limits<double>::quiet_NaN();

    if (B == 8) {
        for (int64_t m = 0; m < n_out; ++m) {
            const int64_t i = idx[m];
            out[m] = (ploidy[i] & 0x80u)
                ? nan
                : static_cast<double>(2u * probs[2 * i] + probs[2 * i + 1]) / 255.0;
        }
    } else {
        const double denom = static_cast<double>((1ull << B) - 1ull);
        for (int64_t m = 0; m < n_out; ++m) {
            const int64_t i = idx[m];
            if (ploidy[i] & 0x80u) { out[m] = nan; continue; }
            const uint64_t bit = static_cast<uint64_t>(i) * 2u * B;
            const uint64_t v0 = read_bits(probs, bit, B);
            const uint64_t v1 = read_bits(probs, bit + B, B);
            out[m] = static_cast<double>(2u * v0 + v1) / denom;
        }
    }
    meta.status = REC_OK;
}

inline py::array_t<double> checked_out(py::array out_any, int64_t n_out, int64_t n_rec) {
    // Taken as a plain array and checked, never converted: a silent copy would
    // swallow every write.
    if (!py::isinstance<py::array_t<double>>(out_any)) {
        throw std::runtime_error("bgen decode: out must be a float64 array");
    }
    auto out = py::reinterpret_borrow<py::array_t<double>>(out_any);
    if (out.ndim() != 2 || out.shape(0) != n_out || out.shape(1) < n_rec) {
        throw std::runtime_error("bgen decode: out must be (len(sample_idx), >= n_records)");
    }
    if (!out.writeable()) throw std::runtime_error("bgen decode: out must be writeable");
    if (out.strides(0) != static_cast<py::ssize_t>(sizeof(double)) ||
        (out.shape(1) > 1 && out.strides(1) != static_cast<py::ssize_t>(sizeof(double) * n_out))) {
        throw std::runtime_error("bgen decode: out must be a contiguous Fortran-order float64 array");
    }
    return out;
}

inline void check_sample_idx(const int64_t* idx, int64_t n_out, uint32_t n_samples_file) {
    for (int64_t m = 0; m < n_out; ++m) {
        if (idx[m] < 0 || idx[m] >= static_cast<int64_t>(n_samples_file)) {
            throw std::runtime_error("bgen decode: sample index out of range");
        }
    }
}

inline py::list meta_list(const std::vector<RecordMeta>& meta) {
    py::list result;
    for (const RecordMeta& m : meta) {
        result.append(py::make_tuple(static_cast<int>(m.status), m.rsid, m.chrom, m.pos,
                                     m.allele0, m.allele1, m.error));
    }
    return result;
}

// Decode a block of raw records (list of bytes) into out[:, 0:len(records)] (Fortran order).
// Returns a list of (status, rsid, chrom, pos, allele0, allele1, error) tuples.
inline py::list decode_block(py::list records,
                             py::array_t<int64_t, py::array::c_style | py::array::forcecast> sample_idx,
                             py::array out_any,
                             int compression, uint32_t n_samples_file, int n_threads) {
    const int64_t n_rec = static_cast<int64_t>(records.size());
    const int64_t n_out = sample_idx.shape(0);
    auto out = checked_out(out_any, n_out, n_rec);
    const int64_t* idx = sample_idx.data();
    check_sample_idx(idx, n_out, n_samples_file);

    std::vector<const uint8_t*> ptr(n_rec);
    std::vector<size_t> len(n_rec);
    for (int64_t r = 0; r < n_rec; ++r) {
        py::bytes b = records[r].cast<py::bytes>();
        char* p = nullptr;
        Py_ssize_t n = 0;
        if (PyBytes_AsStringAndSize(b.ptr(), &p, &n) != 0) throw py::error_already_set();
        ptr[r] = reinterpret_cast<const uint8_t*>(p);
        len[r] = static_cast<size_t>(n);
    }
    std::vector<RecordMeta> meta(n_rec);
    double* out_ptr = out.mutable_data();

    {
        py::gil_scoped_release release;
#ifdef _OPENMP
        const int nt = n_threads > 0 ? n_threads : 1;
#pragma omp parallel num_threads(nt)
#endif
        {
            std::vector<uint8_t> scratch;
            Inflater inf;
#ifdef _OPENMP
#pragma omp for schedule(dynamic, 1)
#endif
            for (int64_t r = 0; r < n_rec; ++r) {
                decode_record(ptr[r], len[r], compression, n_samples_file, idx, n_out,
                              out_ptr + static_cast<size_t>(r) * n_out, scratch, inf, meta[r]);
            }
        }
    }
    return meta_list(meta);
}

// Same, but each thread preads its records from `fd` itself, so a block keeps
// n_threads reads in flight: on a network filesystem a sparse selection of
// records (e.g. HapMap3 out of the full imputed BGEN) is otherwise
// latency-bound at one outstanding read.
inline py::list decode_file_block(int fd,
                                  py::array_t<int64_t, py::array::c_style | py::array::forcecast> offsets,
                                  py::array_t<int64_t, py::array::c_style | py::array::forcecast> sizes,
                                  py::array_t<int64_t, py::array::c_style | py::array::forcecast> sample_idx,
                                  py::array out_any,
                                  int compression, uint32_t n_samples_file, int n_threads) {
    const int64_t n_rec = offsets.shape(0);
    if (sizes.shape(0) != n_rec) throw std::runtime_error("bgen decode: offsets and sizes differ in length");
    const int64_t n_out = sample_idx.shape(0);
    auto out = checked_out(out_any, n_out, n_rec);
    const int64_t* idx = sample_idx.data();
    check_sample_idx(idx, n_out, n_samples_file);
    const int64_t* off = offsets.data();
    const int64_t* sz = sizes.data();
    std::vector<RecordMeta> meta(n_rec);
    double* out_ptr = out.mutable_data();

    {
        py::gil_scoped_release release;
#ifdef _OPENMP
        const int nt = n_threads > 0 ? n_threads : 1;
#pragma omp parallel num_threads(nt)
#endif
        {
            std::vector<uint8_t> raw, scratch;
            Inflater inf;
#ifdef _OPENMP
#pragma omp for schedule(dynamic, 1)
#endif
            for (int64_t r = 0; r < n_rec; ++r) {
                if (sz[r] <= 0 || off[r] < 0) {
                    meta[r].status = REC_CORRUPT; meta[r].error = "bad .bgi offset/size"; continue;
                }
                raw.resize(static_cast<size_t>(sz[r]));
                size_t got = 0;
                bool io_ok = true;
                while (got < raw.size()) {
                    const ssize_t k = ::pread(fd, raw.data() + got, raw.size() - got,
                                              static_cast<off_t>(off[r] + static_cast<int64_t>(got)));
                    if (k < 0 && errno == EINTR) continue;
                    if (k <= 0) { io_ok = false; break; }
                    got += static_cast<size_t>(k);
                }
                if (!io_ok) {
                    meta[r].status = REC_CORRUPT; meta[r].error = "short read"; continue;
                }
                decode_record(raw.data(), raw.size(), compression, n_samples_file, idx, n_out,
                              out_ptr + static_cast<size_t>(r) * n_out, scratch, inf, meta[r]);
            }
        }
    }
    return meta_list(meta);
}

}  // namespace gmoq_bgen

#endif  // GMOQ_STAGE1_BGEN_DECODE_HPP
