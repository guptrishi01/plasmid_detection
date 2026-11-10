"""Findings tables from a 92-run pipeline --outdir: where MAGs and resistance genes are,
and which of them "travel" from pre-treatment water through the plant into the creek.

    python scripts/findings/findings.py <outdir> [--config config/config.yaml]

Writes <outdir>/findings/. Site categories come from params.stats.site_groups (PI,
2026-11-10): pre_treatment, plant, creek_downstream, creek_upstream, control.

Definitions
- A MAG is present in a run if its CoverM mean coverage there is >= the threshold
  (1x; 0.5x and 2x reported for sensitivity), present in a category if present in any of
  its runs.
- A MAG travels if present in pre-treatment and creek downstream, and in no creek
  upstream run (upstream = the creek's own background); "wastewater-specific" travellers
  are also absent from the remote control streams.
- A resistance gene = AMRFinderPlus Type AMR or any RGI hit (params.rollup.amrfinder_types).
- Plasmids are compared through MOB-typer's MOB-suite reference clusters. mob_typer copies
  the cluster ids of a contig's nearest reference plasmid whatever the distance
  (mob_typer.py), so a contig counts for a cluster only when it is close to that
  reference: Mash distance <= 0.025 for the secondary cluster ("same plasmid") and <= 0.06
  for the primary cluster ("plasmid family"), MOB-suite's own cluster distances. geNomad
  plasmid contigs are fragments: a shared cluster means a shared plasmid backbone.
- A resistance plasmid cluster travels if its resistance-carrying contigs occur in
  pre-treatment and creek downstream runs, and in no upstream run.
"""

import argparse
import csv
import statistics
import sys
from collections import Counter, defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))
from plasmid_detection import config  # noqa: E402

ORDER = ["pre_treatment", "plant", "creek_downstream", "creek_upstream", "control"]
THRESHOLDS = [0.5, 1.0, 2.0]
# MOB-suite cluster levels: (name, mob_typer column, max Mash distance to the reference)
LEVELS = [("secondary", "secondary_cluster_id", 0.025), ("primary", "primary_cluster_id", 0.06)]


def read_tsv(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def write_tsv(path, rows, fields):
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


def rank(tax, r):
    for part in tax.split(";"):
        if part.startswith(f"{r}__"):
            return part[3:] or ""
    return ""


def is_resistance(hit, types):
    return hit["source"] != "amrfinder" or hit["detail"] in types


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("outdir")
    ap.add_argument("--config", default=None)
    a = ap.parse_args(argv)
    cfg = config.load(a.config) if a.config else config.load()
    out = Path(a.outdir)
    fdir = out / "findings"
    fdir.mkdir(exist_ok=True)
    site_cat = {s: spec.split(":")[0] for spec in cfg["params"]["stats"]["site_groups"]
                for s in spec.split(":")[1].split(",")}
    types = cfg["params"]["rollup"]["amrfinder_types"].split()

    sheet = {r["run"]: r for r in read_tsv(out / "metadata" / "samples.tsv")}
    cat = {r: site_cat[v["site"]] for r, v in sheet.items()}
    runs_in = {c: sorted(r for r in sheet if cat[r] == c) for c in ORDER}

    # ---- MAGs: abundance and presence per category
    mags = {r["mag"]: r for r in read_tsv(out / "tables" / "07_mags.tsv")}
    cov, rel = {}, {}
    for r in read_tsv(out / "07_mags" / "coverm_genome.tsv"):
        g = r.pop("Genome")
        if g == "unmapped":
            continue
        cov[g] = {k.split("_1.fastq")[0]: float(v) for k, v in r.items() if k.endswith(" Mean")}
        rel[g] = {k.split("_1.fastq")[0]: float(v) for k, v in r.items() if "Relative Abundance" in k}
    amr_hits = [h for h in read_tsv(out / "tables" / "08_mag_amr.tsv") if h["Type"] in types]
    mag_genes = {}
    for h in amr_hits:
        mag_genes.setdefault(h["Name"], set()).add(h["Element symbol"])

    def present(g, t):
        return {c for c in ORDER if any(cov[g][r] >= t for r in runs_in[c])}

    rows = []
    for g, m in sorted(mags.items()):
        tax = m["gtdb_classification"]
        p1 = present(g, 1.0)
        row = {"mag": g, "phylum": rank(tax, "p"), "class": rank(tax, "c"), "genus": rank(tax, "g"),
               "species": rank(tax, "s"), "tier": m["tier"], "binning_group": m["group"],
               "resistance_genes": ",".join(sorted(mag_genes.get(g, ())))}
        for c in ORDER:
            # 6 decimals: slides round again, and a 4-decimal value can flip their last digit
            row[f"mean_relabund_{c}"] = round(statistics.mean(rel[g][r] for r in runs_in[c]), 6)
            row[f"runs_present_{c}"] = sum(cov[g][r] >= 1.0 for r in runs_in[c])
        row["present_in"] = ",".join(c for c in ORDER if c in p1)
        row["travels"] = "yes" if {"pre_treatment", "creek_downstream"} <= p1 and "creek_upstream" not in p1 else ""
        row["wastewater_specific"] = "yes" if row["travels"] and "control" not in p1 else ""
        rows.append(row)
    fields = ["mag", "phylum", "class", "genus", "species", "tier", "binning_group", "resistance_genes",
              "present_in", "travels", "wastewater_specific"] + [f"runs_present_{c}" for c in ORDER] + [f"mean_relabund_{c}" for c in ORDER]
    write_tsv(fdir / "mag_presence.tsv", rows, fields)

    summary = []
    for t in THRESHOLDS:
        pres = {g: present(g, t) for g in mags}
        s = {"min_mean_coverage": t}
        for c in ORDER:
            s[f"mags_in_{c}"] = sum(c in p for p in pres.values())
        s["only_pre_treatment_or_plant"] = sum(1 for p in pres.values() if p and p <= {"pre_treatment", "plant"})
        s["travel_mags"] = sum(1 for p in pres.values() if {"pre_treatment", "creek_downstream"} <= p and "creek_upstream" not in p)
        s["travel_mags_through_plant"] = sum(1 for p in pres.values() if {"pre_treatment", "plant", "creek_downstream"} <= p and "creek_upstream" not in p)
        s["travel_mags_not_in_control"] = sum(1 for p in pres.values() if {"pre_treatment", "creek_downstream"} <= p and not p & {"creek_upstream", "control"})
        s["resistance_mags"] = len(mag_genes)
        s["resistance_mags_in_downstream"] = sum("creek_downstream" in pres[g] for g in mag_genes)
        s["resistance_mags_travel"] = sum(1 for g in mag_genes if {"pre_treatment", "creek_downstream"} <= pres[g] and "creek_upstream" not in pres[g])
        s["resistance_mags_travel_not_in_control"] = sum(1 for g in mag_genes if {"pre_treatment", "creek_downstream"} <= pres[g] and not pres[g] & {"creek_upstream", "control"})
        summary.append(s)
    write_tsv(fdir / "presence_summary.tsv", summary, list(summary[0]))

    # ---- resistance genes on contigs and plasmids, per run
    arg_plasmids = []
    per_run, clusters = [], defaultdict(lambda: {"contigs": 0, "arg_contigs": 0, "cats": set(), "arg_cats": set(),
                                                 "arg_runs": set(), "genes": Counter(), "rep": Counter(),
                                                 "mobility": Counter(), "host_range": Counter()})
    host = {(r["run"], r["plasmid_contig"]): r for r in read_tsv(out / "tables" / "08_plasmid_host_inference.tsv")}
    for run in sorted(sheet):
        d = out / "05_contigs" / run
        tag = [h for h in read_tsv(d / "amr_tagged.tsv") if is_resistance(h, types)]
        genes_on = defaultdict(set)
        for h in tag:
            genes_on[h["contig"]].add(h["gene"])
        plasmid_arg = {h["contig"] for h in tag if h["tag"] == "plasmid"}
        per_run.append({"run": run, "category": cat[run], "site": sheet[run]["site"],
                        "location": sheet[run]["location"], "timepoint": sheet[run]["timepoint"],
                        "resistance_hits": len(tag),
                        "resistance_hits_on_plasmids": sum(h["tag"] == "plasmid" for h in tag),
                        "resistance_genes_distinct": len({h["gene"] for h in tag})})
        # every plasmid contig carrying a resistance gene, with its host inferences
        for contig in sorted({h["contig"] for h in tag if h["tag"] == "plasmid"}):
            hi = host.get((run, contig), {})
            arg_plasmids.append({"run": run, "category": cat[run], "contig": contig,
                                 "genes": ",".join(sorted(genes_on[contig])),
                                 "mob_mobility": hi.get("mob_mobility", ""),
                                 "mob_host_range": hi.get("mob_host_range_name", ""),
                                 "cobin_genus": rank(hi.get("cobin_taxonomy", ""), "g"),
                                 "corr_genus": rank(hi.get("corr_best_taxonomy", ""), "g"),
                                 "corr_rho": hi.get("corr_best_rho", "")})
        mob = d / "mob_typer.tsv"
        if not mob.exists():
            continue
        for m in read_tsv(mob):
            contig = m["sample_id"].split()[0]
            try:
                dist = float(m["mash_neighbor_distance"])
            except ValueError:
                continue
            for level, col, max_dist in LEVELS:
                cid = m[col]
                if not cid or cid == "-" or dist > max_dist:
                    continue
                c = clusters[(level, cid)]
                c["contigs"] += 1
                c["cats"].add(cat[run])
                g = genes_on.get(contig, set()) if contig in plasmid_arg else set()
                if g:
                    c["arg_contigs"] += 1
                    c["arg_cats"].add(cat[run])
                    c["arg_runs"].add(run)
                    c["genes"].update(g)
                    c["rep"][m["rep_type(s)"]] += 1
                    c["mobility"][m["predicted_mobility"]] += 1
                    c["host_range"][m["predicted_host_range_overall_name"]] += 1
                    hi = host.get((run, contig), {})
                    for k in ("cobin_taxonomy", "corr_best_taxonomy"):
                        if hi.get(k):
                            c.setdefault(k, Counter())[rank(hi[k], "g")] += 1
    write_tsv(fdir / "resistance_by_run.tsv", per_run, list(per_run[0]))
    write_tsv(fdir / "resistance_plasmid_contigs.tsv", arg_plasmids, list(arg_plasmids[0]))

    crow = []
    for (level, cid), c in clusters.items():
        if not c["arg_contigs"]:
            continue
        top = lambda cnt: ";".join(f"{k}({v})" for k, v in cnt.most_common(4)) if cnt else ""  # noqa: E731
        crow.append({"level": level, "cluster": cid,
                     "plasmid_contigs": c["contigs"], "with_resistance_genes": c["arg_contigs"],
                     "categories_any": ",".join(x for x in ORDER if x in c["cats"]),
                     "categories_with_resistance": ",".join(x for x in ORDER if x in c["arg_cats"]),
                     "runs_with_resistance": len(c["arg_runs"]),
                     "travels": "yes" if {"pre_treatment", "creek_downstream"} <= c["arg_cats"] and "creek_upstream" not in c["arg_cats"] else "",
                     "genes": top(c["genes"]), "rep_types": top(c["rep"]), "mobility": top(c["mobility"]),
                     "mob_host_range": top(c["host_range"]), "cobin_host_genera": top(c.get("cobin_taxonomy", Counter())),
                     "corr_host_genera": top(c.get("corr_best_taxonomy", Counter()))})
    crow.sort(key=lambda r: (r["level"] != "secondary", -int(r["travels"] == "yes"), -r["runs_with_resistance"]))
    write_tsv(fdir / "resistance_plasmid_clusters.tsv", crow, list(crow[0]))
    lev = []
    for level, _col, max_dist in LEVELS:
        arg = [c for (lv, _), c in clusters.items() if lv == level and c["arg_contigs"]]
        lev.append({"cluster_level": level, "max_mash_distance": max_dist,
                    "clusters": sum(lv == level for lv, _ in clusters),
                    "resistance_contigs_assigned": sum(c["arg_contigs"] for c in arg),
                    "resistance_contigs_total": len(arg_plasmids),
                    "clusters_with_resistance": len(arg),
                    **{f"with_resistance_in_{c}": sum(c in x["arg_cats"] for x in arg) for c in ORDER},
                    "travel": sum(1 for x in arg if {"pre_treatment", "creek_downstream"} <= x["arg_cats"] and "creek_upstream" not in x["arg_cats"]),
                    "travel_through_plant": sum(1 for x in arg if {"pre_treatment", "plant", "creek_downstream"} <= x["arg_cats"] and "creek_upstream" not in x["arg_cats"])})
    write_tsv(fdir / "plasmid_cluster_summary.tsv", lev, list(lev[0]))
    print(f"wrote {fdir}: {len(rows)} MAGs; " + "; ".join(
        f"{x['cluster_level']}: {x['clusters_with_resistance']} resistance clusters, {x['travel']} travel" for x in lev))
    return 0


if __name__ == "__main__":
    sys.exit(main())
