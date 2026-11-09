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


def test_resistance_hit_keeps_amr_type_and_all_rgi():
    amr = {"source": "amrfinder", "detail": "AMR"}
    stress = {"source": "amrfinder", "detail": "STRESS"}
    rgi = {"source": "rgi", "detail": "Strict"}
    assert tables.is_resistance_hit(amr, ["AMR"])
    assert not tables.is_resistance_hit(stress, ["AMR"])
    assert tables.is_resistance_hit(stress, ["AMR", "STRESS"])
    assert tables.is_resistance_hit(rgi, ["AMR"])


def test_ambiguous_fraction_sums_merged_ranges():
    assert tables.ambiguous_fraction("", 1000) == 0
    assert tables.ambiguous_fraction("1-100", 1000) == 0.1
    assert tables.ambiguous_fraction("350-357,1642-3893", 4000) == (8 + 2252) / 4000


def _drep_tables(d, clusters):
    """clusters: {secondary_cluster: [winner, other, ...]} -> dRep Cdb.csv + Wdb.csv"""
    t = d / "data_tables"
    t.mkdir(parents=True)
    cdb = ["genome,secondary_cluster"] + [f"{g},{c}" for c, gs in clusters.items() for g in gs]
    wdb = ["genome,cluster"] + [f"{gs[0]},{c}" for c, gs in clusters.items()]
    (t / "Cdb.csv").write_text("\n".join(cdb) + "\n")
    (t / "Wdb.csv").write_text("\n".join(wdb) + "\n")


def test_catalog_links_chain_group_then_catalog(tmp_path):
    # upstream: a1 + b1 are one species (a1 wins); influent: c1 alone, d1 alone.
    # catalog: a1 and c1 are one species (c1 wins); d1 alone.
    _drep_tables(tmp_path / "groups/upstream/drep", {"1_1": ["a1.fa", "b1.fa"]})
    _drep_tables(tmp_path / "groups/influent/drep", {"1_1": ["c1.fa"], "2_1": ["d1.fa"]})
    _drep_tables(tmp_path / "drep", {"1_1": ["c1.fa", "a1.fa"], "2_1": ["d1.fa"]})
    links = tables.catalog_links(tmp_path)
    assert links["b1.fa"] == {"group": "upstream", "group_rep": "a1.fa", "catalog_rep": "c1.fa"}
    assert links["a1.fa"]["catalog_rep"] == "c1.fa"
    assert links["c1.fa"] == {"group": "influent", "group_rep": "c1.fa", "catalog_rep": "c1.fa"}
    assert links["d1.fa"]["catalog_rep"] == "d1.fa"
    assert len(links) == 4


def test_catalog_links_single_level_layout(tmp_path):
    _drep_tables(tmp_path / "drep", {"1_1": ["a1.fa", "b1.fa"]})
    links = tables.catalog_links(tmp_path)
    assert links["b1.fa"] == {"group": "", "group_rep": "b1.fa", "catalog_rep": "a1.fa"}
