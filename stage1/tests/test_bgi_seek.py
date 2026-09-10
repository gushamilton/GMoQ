import sqlite3

import pytest

from fungwas_stage1.cli import SeekingVariantGenerator, query_bgi_records_for_rsids


def _make_bgi(tmp_path, rows):
    """rows: (chromosome, position, rsid, allele1, allele2, file_start_position)"""
    bgen_path = tmp_path / "chr22.bgen"
    bgen_path.write_bytes(b"")
    bgi_path = tmp_path / "chr22.bgen.bgi"
    conn = sqlite3.connect(bgi_path)
    try:
        conn.execute(
            "CREATE TABLE Variant ("
            "chromosome TEXT NOT NULL, position INT NOT NULL, rsid TEXT NOT NULL, "
            "number_of_alleles INT NOT NULL, allele1 TEXT NOT NULL, allele2 TEXT NULL, "
            "file_start_position INT NOT NULL, size_in_bytes INT NOT NULL, "
            "PRIMARY KEY (chromosome, position, rsid, allele1, allele2, file_start_position))"
        )
        conn.executemany(
            "INSERT INTO Variant VALUES (?, ?, ?, 2, ?, ?, ?, 100)",
            [(c, p, r, a1, a2, off) for (c, p, r, a1, a2, off) in rows],
        )
        conn.commit()
    finally:
        conn.close()
    return str(bgen_path)


ROWS = [
    # chrom, pos, rsid, allele1, allele2, file_start_position
    ("22", 300, "rs_c", "A", "G", 30),
    ("22", 100, "rs_a", "A", "G", 10),
    ("22", 200, "rs_b", "C", "A", 20),
    ("22", 200, "rs_multi", "C", "T", 25),
    ("22", 200, "rs_multi", "C", "G", 26),
    ("22", 400, "rs_other", "G", "A", 40),
]


def test_query_returns_requested_records_in_file_order(tmp_path):
    bgen = _make_bgi(tmp_path, ROWS)
    got = query_bgi_records_for_rsids(bgen, {"rs_a", "rs_b", "rs_c", "rs_multi"})

    # ordered by file_start_position, and only the requested rsIDs
    assert [r[0] for r in got] == [10, 20, 25, 26, 30]
    assert [r[3] for r in got] == ["rs_a", "rs_b", "rs_multi", "rs_multi", "rs_c"]
    assert all(r[1] == "22" for r in got)
    assert [r[2] for r in got] == [100, 200, 200, 200, 300]


def test_query_ignores_rsids_absent_from_the_index(tmp_path):
    bgen = _make_bgi(tmp_path, ROWS)
    got = query_bgi_records_for_rsids(bgen, {"rs_a", "rs_nonexistent"})
    assert [r[3] for r in got] == ["rs_a"]


def test_query_empty_request_is_empty(tmp_path):
    bgen = _make_bgi(tmp_path, ROWS)
    assert query_bgi_records_for_rsids(bgen, set()) == []


def test_query_requires_the_index(tmp_path):
    missing = tmp_path / "absent.bgen"
    missing.write_bytes(b"")
    with pytest.raises(FileNotFoundError):
        query_bgi_records_for_rsids(str(missing), {"rs_a"})


class _FakeVar:
    def __init__(self, rsid, offset):
        self.rsid = rsid
        self.fileoffset = offset


class _FakeBgen:
    """Stands in for BgenReader.fetch, returning records at one position.

    Deliberately returns them in primary-key order (by rsid/allele, not by file offset)
    and includes a record that was not requested, which is what the real index does.
    """

    def __init__(self, by_pos):
        self.by_pos = by_pos
        self.calls = []

    def fetch(self, chrom, start, stop):
        self.calls.append((chrom, start, stop))
        for rsid, off in self.by_pos.get((chrom, start), []):
            yield _FakeVar(rsid, off)


def test_generator_yields_file_order_and_filters_by_offset():
    bfile = _FakeBgen({
        ("22", 100): [("rs_a", 10)],
        # same position, out of file order, plus an untargeted neighbour
        ("22", 200): [("rs_b", 22), ("rs_bystander", 21), ("rs_b2", 20)],
        ("22", 300): [("rs_c", 30)],
    })
    records = [
        (10, "22", 100, "rs_a"),
        (20, "22", 200, "rs_b2"),
        (22, "22", 200, "rs_b"),
        (30, "22", 300, "rs_c"),
    ]
    snp_set = {"rs_a", "rs_b", "rs_b2", "rs_c"}

    got = list(SeekingVariantGenerator(bfile, records, snp_set))

    # file order preserved, bystander dropped, co-located targets sorted by offset
    assert [v.rsid for v in got] == ["rs_a", "rs_b2", "rs_b", "rs_c"]
    assert [v.fileoffset for v in got] == [10, 20, 22, 30]
    # one indexed lookup per distinct position, not one per rsID
    assert bfile.calls == [("22", 100, 100), ("22", 200, 200), ("22", 300, 300)]


def test_generator_stops_at_target_count():
    bfile = _FakeBgen({
        ("22", 100): [("rs_a", 10)],
        ("22", 200): [("rs_b", 20)],
        ("22", 300): [("rs_c", 30)],
    })
    records = [(10, "22", 100, "rs_a"), (20, "22", 200, "rs_b"), (30, "22", 300, "rs_c")]

    got = list(SeekingVariantGenerator(bfile, records, {"rs_a", "rs_b"}))

    assert [v.rsid for v in got] == ["rs_a", "rs_b"]
    # stopped before seeking to the third position
    assert bfile.calls == [("22", 100, 100), ("22", 200, 200)]


def test_generator_handles_no_records():
    bfile = _FakeBgen({})
    assert list(SeekingVariantGenerator(bfile, [], set())) == []
    assert bfile.calls == []
