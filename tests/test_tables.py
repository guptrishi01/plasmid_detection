import argparse

import pytest

from plasmid_detection import tables


def test_tag_contig_virus_provirus_plasmid_chromosome():
    plasmids, viruses = {"p1"}, {"v1"}
    proviruses = {"c1": [(100, 200)]}
    assert tables.tag_contig("v1", 1, 10, plasmids, viruses, proviruses) == "virus"
    assert tables.tag_contig("c1", 150, 400, plasmids, viruses, proviruses) == "virus"
    assert tables.tag_contig("c1", 300, 400, plasmids, viruses, proviruses) == "chromosome"
    assert tables.tag_contig("p1", 1, 10, plasmids, viruses, proviruses) == "plasmid"
    assert tables.tag_contig("x", 1, 10, plasmids, viruses, proviruses) == "chromosome"


def test_amr_counts_reads_the_count_block_not_the_depth_block(tmp_path):
    # Layout verified with samtools 1.24 `bedcov -c` on two BAMs: 4 BED columns,
    # then one depth sum per BAM, then one read count per BAM.
    bedcov = tmp_path / "bedcov.txt"
    bedcov.write_text("c1\t89\t120\thitA\t20\t0\t2\t0\nc1\t489\t520\thitB\t0\t10\t0\t1\n")
    out = tmp_path / "counts.tsv"
    tables.main(
        [
            "amr-counts",
            "--bedcov",
            str(bedcov),
            "--bams",
            "x/SRR1.bam",
            "y/SRR2.bam",
            "--out",
            str(out),
        ]
    )
    rows = tables.read_tsv(out)
    assert rows == [
        {"hit_id": "hitA", "reads_SRR1": "2", "reads_SRR2": "0"},
        {"hit_id": "hitB", "reads_SRR1": "0", "reads_SRR2": "1"},
    ]


def test_mimag_tiers():
    a = argparse.Namespace(high_comp=90, high_cont=5, medium_comp=50, medium_cont=10)
    assert tables.mimag_tier(95, 1, a) == "high"
    assert tables.mimag_tier(95, 5, a) == "medium"  # contamination must be strictly < 5
    assert tables.mimag_tier(50, 9.9, a) == "medium"
    assert tables.mimag_tier(49.9, 1, a) == "low"
    assert tables.mimag_tier(99, 10, a) == "low"


def test_spearman_handles_ties_and_constant_input():
    assert tables.ranks([10, 20, 20, 30]) == [1.0, 2.5, 2.5, 4.0]
    assert tables.spearman([1, 2, 3, 4], [10, 20, 30, 40]) == pytest.approx(1.0)
    assert tables.spearman([1, 2, 3, 4], [4, 3, 2, 1]) == pytest.approx(-1.0)
    rho = tables.spearman([1, 2, 3], [5, 5, 5])
    assert rho != rho  # NaN when one profile is constant


def test_encode_censored_values():
    assert tables.encode_censored("ND", "na", 0.5) == "NA"
    assert tables.encode_censored("<LOQ", "zero", 0.5) == 0
    assert tables.encode_censored("ND", "half_min", 0.5) == 0.5
    assert tables.encode_censored("2.5", "na", 0.5) == 2.5
