import sqlite3

from fungwas_stage1.cli import find_multiallelic_rsids_from_bgi


def test_find_multiallelic_rsids_from_bgi_filters_requested_rsids(tmp_path):
    bgen_path = tmp_path / "chr22.bgen"
    bgen_path.write_bytes(b"")
    bgi_path = tmp_path / "chr22.bgen.bgi"

    conn = sqlite3.connect(bgi_path)
    try:
        conn.execute(
            "CREATE TABLE Variant ("
            "chromosome TEXT, position INTEGER, rsid TEXT, "
            "number_of_alleles INTEGER, allele1 TEXT, allele2 TEXT)"
        )
        conn.executemany(
            "INSERT INTO Variant VALUES (?, ?, ?, ?, ?, ?)",
            [
                ("22", 100, "rs_unique", 2, "A", "G"),
                ("22", 200, "rs_multi", 2, "C", "A"),
                ("22", 200, "rs_multi", 2, "C", "T"),
                ("22", 300, "rs_not_requested", 2, "G", "A"),
                ("22", 300, "rs_not_requested", 2, "G", "T"),
            ],
        )
        conn.commit()
    finally:
        conn.close()

    out = find_multiallelic_rsids_from_bgi(
        str(bgen_path),
        requested_rsids={"rs_unique", "rs_multi"},
    )

    assert out == {"rs_multi"}


def test_find_multiallelic_rsids_from_bgi_can_scan_all_rsids(tmp_path):
    bgen_path = tmp_path / "chr22.bgen"
    bgen_path.write_bytes(b"")
    bgi_path = tmp_path / "chr22.bgen.bgi"

    conn = sqlite3.connect(bgi_path)
    try:
        conn.execute(
            "CREATE TABLE Variant ("
            "chromosome TEXT, position INTEGER, rsid TEXT, "
            "number_of_alleles INTEGER, allele1 TEXT, allele2 TEXT)"
        )
        conn.executemany(
            "INSERT INTO Variant VALUES (?, ?, ?, ?, ?, ?)",
            [
                ("22", 100, "rs1", 2, "A", "G"),
                ("22", 200, "rs2", 2, "C", "A"),
                ("22", 200, "rs2", 2, "C", "T"),
            ],
        )
        conn.commit()
    finally:
        conn.close()

    assert find_multiallelic_rsids_from_bgi(str(bgen_path)) == {"rs2"}
