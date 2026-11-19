"""scripts/coassembly/coasm_tables.py on a hand-built example with known answers."""

import csv
import importlib.util
from pathlib import Path

SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "coassembly" / "coasm_tables.py"
spec = importlib.util.spec_from_file_location("coasm_tables", SCRIPT)
ct = importlib.util.module_from_spec(spec)
spec.loader.exec_module(ct)


def w(path, header, rows, delim="\t"):
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as fh:
        cw = csv.writer(fh, delimiter=delim, lineterminator="\n")
        cw.writerow(header)
        cw.writerows(rows)


def r(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def build(tmp):
    # two co-assembly groups: wet (sites A, B) and dry (site C); site groups AB:A,B and C:C
    w(tmp / "samples.tsv", ["run", "location", "site", "timepoint"],
      [["R1", "L", "A", "1"], ["R2", "L", "B", "1"], ["R3", "L", "C", "1"]])
    w(tmp / "mimag/wet.tsv", ["run", "bin", "completeness", "contamination", "tier"],
      [["wet", "wet.bin.1", "95", "1", "high"], ["wet", "wet.bin.2", "60", "4", "medium"]])
    w(tmp / "mimag/dry.tsv", ["run", "bin", "completeness", "contamination", "tier"],
      [["dry", "dry.bin.1", "70", "2", "medium"]])
    ra, mean = "Relative Abundance (%)", "Mean"
    w(tmp / "coverm/wet.tsv",
      ["Genome", f"R1_1.fastq.gz {ra}", f"R1_1.fastq.gz {mean}", f"R2_1.fastq.gz {ra}", f"R2_1.fastq.gz {mean}"],
      [["unmapped", "60", "0", "90", "0"], ["wet.bin.1", "30", "5", "10", "0.5"], ["wet.bin.2", "10", "2", "0", "0"]])
    w(tmp / "coverm/dry.tsv", ["Genome", f"R3_1.fastq.gz {ra}", f"R3_1.fastq.gz {mean}"],
      [["unmapped", "75", "0"], ["dry.bin.1", "25", "3"]])
    w(tmp / "group_mags.tsv", ["group", "mag"], [["wet", "wet.bin.1"], ["wet", "wet.bin.2"], ["dry", "dry.bin.1"]])
    w(tmp / "gtdbtk/gtdbtk.bac120.summary.tsv", ["user_genome", "classification"],
      [["wet.bin.1", "d__Bacteria;p__P;c__C;o__O;f__F;g__Aeromonas;s__Aeromonas hydrophila"],
       ["wet.bin.2", "d__Bacteria;p__P;c__C;o__O;f__F;g__Zoogloea;s__"],
       ["dry.bin.1", "d__Bacteria;p__P;c__C;o__O;f__F;g__Rhodoluna;s__"]])
    w(tmp / "mag_amr.tsv", ["Name", "Element symbol", "Type"],
      [["wet.bin.1", "blaOXA", "AMR"], ["wet.bin.1", "qacE", "STRESS"]])
    # species: c1 = wet.bin.1 + site S1 (AB); c2 = wet.bin.2 alone; c3 = dry.bin.1 + S3 (AB, other group);
    # c4 = S2 (AB) alone; c5 = S4 (C) alone
    w(tmp / "cmp/Cdb.csv", ["genome", "secondary_cluster"],
      [["wet.bin.1.fa", "1_1"], ["S1.fa", "1_1"], ["wet.bin.2.fa", "2_1"], ["dry.bin.1.fa", "3_1"],
       ["S3.fa", "3_1"], ["S2.fa", "4_1"], ["S4.fa", "5_1"]], ",")
    w(tmp / "cmp/Ndb.csv", ["reference", "querry", "ani"], [["S1.fa", "wet.bin.1.fa", "0.991"]], ",")
    w(tmp / "genomes.tsv", ["genome", "set", "group"],
      [["wet.bin.1.fa", "coassembly", ""], ["wet.bin.2.fa", "coassembly", ""], ["dry.bin.1.fa", "coassembly", ""],
       ["S1.fa", "site", "AB"], ["S2.fa", "site", "AB"], ["S3.fa", "site", "AB"], ["S4.fa", "site", "C"]])
    w(tmp / "site_group_mags.tsv",
      ["group", "genome", "run", "completeness", "contamination", "tier", "group_mag", "catalog_mag", "gtdb_classification"],
      [["AB", "S1", "R1", "90", "2", "high", "S1", "S1", "d__Bacteria;g__Aeromonas;s__Aeromonas hydrophila"],
       ["AB", "S2", "R2", "55", "1", "medium", "S2", "S2", ""],
       ["AB", "S3", "R1", "80", "1", "medium", "S3", "S3", ""],
       ["C", "S4", "R3", "65", "3", "medium", "S4", "S4", ""]])
    w(tmp / "site_amr.tsv", ["Name", "Element symbol", "Type"], [["S1", "blaOXA", "AMR"], ["S2", "tet", "AMR"]])
    ct.main(["--sheet", str(tmp / "samples.tsv"), "--groups", "wet:A,B", "dry:C", "--site-groups", "AB:A,B", "C:C",
             "--mimag-dir", str(tmp / "mimag"), "--coverm-dir", str(tmp / "coverm"),
             "--mag-list", str(tmp / "group_mags.tsv"), "--gtdbtk", str(tmp / "gtdbtk"),
             "--mag-amr", str(tmp / "mag_amr.tsv"), "--compare", str(tmp / "cmp"),
             "--genomes", str(tmp / "genomes.tsv"), "--site-mags", str(tmp / "site_group_mags.tsv"),
             "--site-mag-amr", str(tmp / "site_amr.tsv"), "--min-coverage", "1", "--out-dir", str(tmp)])


def test_run_hits(tmp_path):
    build(tmp_path)
    hits = {x["run"]: x for x in r(tmp_path / "coasm_run_hits.tsv")}
    assert hits["R1"]["group"] == "wet" and hits["R1"]["pct_reads_on_group_mags"] == "40.0000"
    assert hits["R1"]["mags_hit"] == "2" and hits["R1"]["mags_hit_with_resistance"] == "1"
    assert hits["R1"]["resistance_genes_in_hit_mags"] == "blaOXA"      # STRESS hit excluded
    assert hits["R2"]["mags_hit"] == "0"                               # 0.5x < 1x
    assert hits["R2"]["top_mag"] == "wet.bin.1" and hits["R2"]["top_mag_label"] == "Aeromonas hydrophila"
    assert hits["R3"]["group"] == "dry" and hits["R3"]["mags_hit"] == "1"


def test_mags_and_compare(tmp_path):
    build(tmp_path)
    mags = {x["mag"]: x for x in r(tmp_path / "coasm_mags.tsv")}
    assert mags["wet.bin.1"]["best_own_site_mag"] == "S1" and mags["wet.bin.1"]["best_own_site_ani"] == "0.9910"
    assert mags["wet.bin.1"]["runs_hit"] == "1" and mags["wet.bin.1"]["mean_relabund_pct"] == "20.000000"
    assert mags["wet.bin.2"]["new_vs_own_sites"] == "yes" and mags["wet.bin.2"]["label"] == "Zoogloea sp."
    # dry.bin.1's species was only binned in a site group of the other co-assembly group
    assert mags["dry.bin.1"]["new_vs_own_sites"] == "yes" and mags["dry.bin.1"]["other_site_groups"] == "AB"
    cmp = {x["group"]: x for x in r(tmp_path / "coasm_compare.tsv")}
    wet, dry = cmp["wet"], cmp["dry"]
    assert (wet["group_species"], wet["site_species"], wet["shared_species"]) == ("2", "3", "1")
    assert (wet["group_only_species"], wet["site_only_species"]) == ("1", "2")
    assert wet["resistance_site_mags"] == "2" and wet["resistance_site_mags_recovered"] == "1"
    assert wet["median_completeness_diff"] == "5.00" and wet["pairs_group_mag_better"] == "1"
    assert (dry["shared_species"], dry["group_only_species"], dry["group_only_species_in_other_sites"]) == ("0", "1", "1")
    site = {x["site_mag"]: x for x in r(tmp_path / "coasm_compare_mags.tsv")}
    assert site["S1"]["recovered_by_group"] == "yes" and site["S3"]["recovered_by_group"] == ""
