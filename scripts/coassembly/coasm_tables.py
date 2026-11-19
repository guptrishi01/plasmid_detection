"""Summary tables of the co-assembly workflow (nextflow/coassembly.nf, COASM_TABLES).

    coasm_mags.tsv          one row per group MAG: quality, GTDB, resistance genes, how many
                            of the group's runs hit it, its species match among the per-site MAGs
    coasm_run_hits.tsv      one row per run: what it hits among its own group's MAGs
    coasm_compare.tsv       one row per group: group MAGs vs the per-site MAGs of its sites
    coasm_compare_mags.tsv  one row per per-site MAG: recovered by its co-assembly group or not

Definitions
- Species = dRep compare secondary cluster (params.coassembly.compare_ani, 95% ANI) over
  the group MAGs and the per-site MAGs (level-1 dRep of the 11 site groups) together.
- A run hits a MAG when its CoverM mean coverage on that MAG >= --min-coverage; each run is
  mapped only to the MAGs of its own group.
- Per-site MAGs "of a group" = those of the site groups (params.binning.groups) whose sites
  are all in that co-assembly group.
"""

import argparse
import csv
import re
import statistics
import sys
from collections import defaultdict
from pathlib import Path


def read_tsv(path, delim="\t"):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter=delim))


def write_tsv(path, rows, cols):
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", lineterminator="\n")
        w.writeheader()
        for r in rows:
            w.writerow({c: r.get(c, "") for c in cols})


def specs(items):
    """["group:SITE,SITE", ...] -> {site: group}, [group, ...] in order"""
    site_of, order = {}, []
    for s in items:
        g, sites = s.split(":", 1)
        order.append(g)
        for x in sites.split(","):
            if x in site_of:
                sys.exit(f"site {x} is in more than one group")
            site_of[x] = g
    return site_of, order


def rank(tax, r):
    for x in (tax or "").split(";"):
        if x.startswith(f"{r}__"):
            return x[3:]
    return ""


def label(tax):
    """species, else genus sp., else the deepest named rank"""
    s = rank(tax, "s")
    if s:
        return s
    for r in "gfocp":
        if rank(tax, r):
            return f"{rank(tax, r)} sp." if r == "g" else f"{rank(tax, r)} ({r})"
    return "unclassified"


def fmt(x, nd=4):
    return f"{x:.{nd}f}"


COVERM_COL = re.compile(r"^(.+?)_1\.fastq(?:\.gz)? (Relative Abundance \(%\)|Mean)$")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--sheet", required=True)
    ap.add_argument("--groups", nargs="+", required=True, help="params.coassembly.groups")
    ap.add_argument("--site-groups", nargs="+", required=True, help="params.binning.groups (11-group run)")
    ap.add_argument("--mimag-dir", required=True, help="<group>.tsv: CheckM2 MIMAG table per group")
    ap.add_argument("--coverm-dir", required=True, help="<group>.tsv: CoverM genome table per group")
    ap.add_argument("--mag-list", required=True, help="group_mags.tsv (group, mag)")
    ap.add_argument("--gtdbtk", required=True, help="GTDB-Tk output dir")
    ap.add_argument("--mag-amr", required=True, help="AMRFinderPlus on the group MAGs (08_mag_amr.tsv)")
    ap.add_argument("--compare", required=True, help="dRep compare data_tables/")
    ap.add_argument("--genomes", required=True, help="genomes.tsv from COASM_COMPARE (genome, set, group)")
    ap.add_argument("--site-mags", required=True, help="11-group tables/07_group_mags.tsv")
    ap.add_argument("--site-mag-amr", required=True, help="11-group tables/08_mag_amr.tsv (catalog MAGs)")
    ap.add_argument("--min-coverage", type=float, default=1.0)
    ap.add_argument("--out-dir", default=".")
    a = ap.parse_args(argv)
    out = Path(a.out_dir)

    group_of_site, groups = specs(a.groups)
    sitegroup_of_site, _ = specs(a.site_groups)
    # site group -> co-assembly group (all of its sites must fall in one co-assembly group)
    cg_of_sg = {}
    for site, sg in sitegroup_of_site.items():
        g = group_of_site.get(site)
        if g is None:
            sys.exit(f"site {site} (site group {sg}) is in no co-assembly group")
        if cg_of_sg.setdefault(sg, g) != g:
            sys.exit(f"site group {sg} is split across co-assembly groups")

    sheet = read_tsv(a.sheet)
    runs_of = defaultdict(list)
    meta = {}
    for r in sheet:
        meta[r["run"]] = r
    # ---- group MAGs: quality, taxonomy, AMR ----------------------------------------------------
    mag_group = {r["mag"]: r["group"] for r in read_tsv(a.mag_list)}
    qual = {}
    for g in groups:
        p = Path(a.mimag_dir, f"{g}.tsv")
        if p.exists():
            for r in read_tsv(p):
                qual[r["bin"]] = r
    tax = {}
    for p in sorted(Path(a.gtdbtk).glob("gtdbtk.*.summary.tsv")):
        for r in read_tsv(p):
            tax[r["user_genome"]] = r["classification"]
    amr = defaultdict(set)
    for r in read_tsv(a.mag_amr):
        if r.get("Type") == "AMR":
            amr[r["Name"]].add(r["Element symbol"])
    missing = [m for m in mag_group if m not in qual]
    if missing:
        sys.exit(f"{len(missing)} group MAGs without a MIMAG row, e.g. {missing[0]}")

    # ---- abundance: each run against its own group's MAGs --------------------------------------
    ra, cov = {}, {}           # (run, mag) -> value
    unmapped = {}
    for g in groups:
        p = Path(a.coverm_dir, f"{g}.tsv")
        if not p.exists():
            continue
        rows = read_tsv(p)
        cols = {}
        for c in rows[0]:
            m = COVERM_COL.match(c)
            if m:
                cols[c] = (m.group(1), "ra" if m.group(2).startswith("Relative") else "mean")
        for run, _k in cols.values():
            if group_of_site.get(meta[run]["site"]) != g:
                sys.exit(f"run {run} in the {g} CoverM table is not a {g} run")
            if run not in runs_of[g]:
                runs_of[g].append(run)
        for r in rows:
            for c, (run, k) in cols.items():
                v = float(r[c])
                if r["Genome"] == "unmapped":
                    if k == "ra":
                        unmapped[run] = v
                elif k == "ra":
                    ra[(run, r["Genome"])] = v
                else:
                    cov[(run, r["Genome"])] = v
    for g in runs_of:
        runs_of[g].sort()

    # ---- species clusters over group MAGs + per-site MAGs -------------------------------------
    cdb = {r["genome"]: r["secondary_cluster"] for r in read_tsv(Path(a.compare, "Cdb.csv"), ",")}
    gen = read_tsv(a.genomes)
    site_sg = {r["genome"][:-3]: r["group"] for r in gen if r["set"] == "site"}
    coasm_set = {r["genome"][:-3] for r in gen if r["set"] == "coassembly"}
    if coasm_set != set(mag_group):
        sys.exit("COASM_COMPARE genomes differ from the group MAG list")
    clus = {g[:-3]: c for g, c in cdb.items()}
    nocl = [m for m in list(coasm_set) + list(site_sg) if m not in clus]
    if nocl:
        sys.exit(f"{len(nocl)} genomes missing from dRep compare Cdb.csv, e.g. {nocl[0]}")
    members = defaultdict(lambda: {"coasm": set(), "site": set()})
    for m in coasm_set:
        members[clus[m]]["coasm"].add(m)
    for m in site_sg:
        members[clus[m]]["site"].add(m)
    # best ANI between a group MAG and a per-site MAG (Ndb, when dRep wrote it)
    ani = {}
    ndb = Path(a.compare, "Ndb.csv")
    if ndb.exists():
        for r in read_tsv(ndb, ","):
            q, ref = r["querry"][:-3], r["reference"][:-3]
            if q != ref:
                ani[(q, ref)] = max(ani.get((q, ref), 0.0), float(r["ani"]))
    sq = {}
    for r in read_tsv(a.site_mags):
        if r["genome"] == r["group_mag"]:
            sq[r["genome"]] = r
    site_arg = defaultdict(set)
    for r in read_tsv(a.site_mag_amr):
        if r.get("Type") == "AMR":
            site_arg[r["Name"]].add(r["Element symbol"])

    def score(q):
        return float(q["completeness"]) - 5 * float(q["contamination"])

    # ---- coasm_mags.tsv ------------------------------------------------------------------------
    mag_rows = []
    for m in sorted(mag_group, key=lambda x: (groups.index(mag_group[x]), x)):
        g = mag_group[m]
        q = qual[m]
        rr = runs_of.get(g, [])
        hits = [r for r in rr if cov.get((r, m), 0.0) >= a.min_coverage]
        vals = [ra.get((r, m), 0.0) for r in rr]
        c = clus[m]
        own = sorted(s for s in members[c]["site"] if cg_of_sg[site_sg[s]] == g)
        other = sorted(s for s in members[c]["site"] if cg_of_sg[site_sg[s]] != g)
        best = max(own, key=lambda s: (score(sq[s]), s)) if own else ""
        mag_rows.append({
            "group": g, "mag": m, "completeness": q["completeness"], "contamination": q["contamination"],
            "tier": q["tier"], "gtdb_classification": tax.get(m, ""), "label": label(tax.get(m, "")),
            "resistance_genes": ",".join(sorted(amr.get(m, ()))),
            "runs_in_group": len(rr), "runs_hit": len(hits),
            "mean_relabund_pct": fmt(statistics.mean(vals), 6) if vals else "",
            "max_relabund_pct": fmt(max(vals), 6) if vals else "",
            "species_cluster": c,
            "new_vs_own_sites": "yes" if not own else "",
            "own_site_mags": len(own), "own_site_groups": ",".join(sorted({site_sg[s] for s in own})),
            "other_site_groups": ",".join(sorted({site_sg[s] for s in other})),
            "best_own_site_mag": best,
            "best_own_site_ani": fmt(ani.get((m, best), ani.get((best, m), 0.0)), 4) if best and ani else "",
            "best_own_site_completeness": sq[best]["completeness"] if best else "",
            "best_own_site_contamination": sq[best]["contamination"] if best else "",
        })
    write_tsv(out / "coasm_mags.tsv", mag_rows, list(mag_rows[0]) if mag_rows else ["group", "mag"])

    # ---- coasm_run_hits.tsv --------------------------------------------------------------------
    run_rows = []
    for g in groups:
        gm = [m for m in mag_group if mag_group[m] == g]
        for run in runs_of.get(g, []):
            hit = [m for m in gm if cov.get((run, m), 0.0) >= a.min_coverage]
            top = max(gm, key=lambda m: (ra.get((run, m), 0.0), m)) if gm else ""
            mr = meta[run]
            run_rows.append({
                "run": run, "group": g, "site": mr["site"], "location": mr["location"], "timepoint": mr["timepoint"],
                "pct_reads_on_group_mags": fmt(100.0 - unmapped.get(run, 100.0), 4),
                "group_mags": len(gm), "mags_hit": len(hit),
                "mags_hit_high": sum(1 for m in hit if qual[m]["tier"] == "high"),
                "mags_hit_with_resistance": sum(1 for m in hit if amr.get(m)),
                "resistance_genes_in_hit_mags": ",".join(sorted(set().union(*(amr.get(m, set()) for m in hit)))) if hit else "",
                "top_mag": top, "top_mag_label": label(tax.get(top, "")) if top else "",
                "top_mag_relabund_pct": fmt(ra.get((run, top), 0.0), 4) if top else "",
            })
    write_tsv(out / "coasm_run_hits.tsv", run_rows, list(run_rows[0]) if run_rows else ["run", "group"])

    # ---- coasm_compare_mags.tsv: every per-site MAG, recovered by its own group or not ---------
    site_rows = []
    for s in sorted(site_sg, key=lambda x: (site_sg[x], x)):
        g = cg_of_sg[site_sg[s]]
        c = clus[s]
        own_c = sorted(m for m in members[c]["coasm"] if mag_group[m] == g)
        q = sq.get(s, {})
        cat = q.get("catalog_mag", "")
        site_rows.append({
            "site_group": site_sg[s], "group": g, "site_mag": s,
            "completeness": q.get("completeness", ""), "contamination": q.get("contamination", ""),
            "tier": q.get("tier", ""), "label": label(q.get("gtdb_classification", "")),
            "catalog_mag": cat, "catalog_resistance_genes": ",".join(sorted(site_arg.get(cat, ()))) if cat == s else "",
            "species_cluster": c, "recovered_by_group": "yes" if own_c else "",
            "group_mags_same_species": ",".join(own_c),
        })
    write_tsv(out / "coasm_compare_mags.tsv", site_rows, list(site_rows[0]) if site_rows else ["site_mag"])

    # ---- coasm_compare.tsv: one row per group --------------------------------------------------
    cmp_rows = []
    for g in groups:
        gm = [m for m in mag_group if mag_group[m] == g]
        sm = [s for s in site_sg if cg_of_sg[site_sg[s]] == g]
        cs = {clus[m] for m in gm}
        ss = {clus[s] for s in sm}
        pairs = [r for r in mag_rows if r["group"] == g and r["best_own_site_mag"]]
        dc = [float(r["completeness"]) - float(r["best_own_site_completeness"]) for r in pairs]
        dx = [float(r["contamination"]) - float(r["best_own_site_contamination"]) for r in pairs]
        argsite = [r for r in site_rows if r["group"] == g and r["catalog_resistance_genes"]]
        cmp_rows.append({
            "group": g, "runs": len(runs_of.get(g, [])),
            "site_groups": ",".join(sorted({site_sg[s] for s in sm})),
            "group_mags": len(gm),
            "group_mags_high": sum(1 for m in gm if qual[m]["tier"] == "high"),
            "group_species": len(cs),
            "site_mags": len(sm),
            "site_mags_high": sum(1 for s in sm if sq.get(s, {}).get("tier") == "high"),
            "site_species": len(ss),
            "shared_species": len(cs & ss),
            "group_only_species": len(cs - ss),
            "group_only_species_in_other_sites": sum(1 for c in cs - ss if members[c]["site"]),
            "site_only_species": len(ss - cs),
            "matched_pairs": len(pairs),
            "median_completeness_diff": fmt(statistics.median(dc), 2) if dc else "",
            "median_contamination_diff": fmt(statistics.median(dx), 2) if dx else "",
            "pairs_group_mag_better": sum(1 for r in pairs if score(qual[r["mag"]]) > score(sq[r["best_own_site_mag"]])),
            "resistance_site_mags": len(argsite),
            "resistance_site_mags_recovered": sum(1 for r in argsite if r["recovered_by_group"]),
            "median_pct_reads_on_group_mags": fmt(statistics.median(
                [float(r["pct_reads_on_group_mags"]) for r in run_rows if r["group"] == g]), 2)
            if any(r["group"] == g for r in run_rows) else "",
        })
    write_tsv(out / "coasm_compare.tsv", cmp_rows, list(cmp_rows[0]))
    for r in cmp_rows:
        print(f"{r['group']}: {r['group_mags']} group MAGs ({r['group_mags_high']} high) vs {r['site_mags']} per-site MAGs; "
              f"species shared {r['shared_species']}, group only {r['group_only_species']}, site only {r['site_only_species']}")


if __name__ == "__main__":
    main()
