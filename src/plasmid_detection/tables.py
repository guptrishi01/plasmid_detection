"""Per-stage summary tables built from tool outputs (stdlib only).

    python -m plasmid_detection.tables <table> [options]

Every table is a TSV with one row per sample (or per MAG / hit), written to --out.
"""

import argparse
import csv
import json
import sys
from pathlib import Path


def read_samples(path):
    return [line.strip() for line in open(path) if line.strip()]


def read_tsv(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def write_tsv(path, rows, fields):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def fastp_stats(fastp_dir, srr):
    report = json.load(open(Path(fastp_dir) / f"{srr}.fastp.json"))
    before = report["summary"]["before_filtering"]
    after = report["summary"]["after_filtering"]
    return {
        "raw_pairs": before["total_reads"] // 2,
        "trimmed_pairs": after["total_reads"] // 2,
        "q30_raw": round(before["q30_rate"], 4),
        "q30_trimmed": round(after["q30_rate"], 4),
        "duplication_rate": round(report["duplication"]["rate"], 4),
        "insert_size_peak": report["insert_size"]["peak"],
    }


def table_qc(args):
    rows = [{"run": s, **fastp_stats(args.fastp_dir, s)} for s in read_samples(args.samples)]
    fields = [
        "run",
        "raw_pairs",
        "trimmed_pairs",
        "q30_raw",
        "q30_trimmed",
        "duplication_rate",
        "insert_size_peak",
    ]
    write_tsv(args.out, rows, fields)


def table_host(args):
    rows = []
    for s in read_samples(args.samples):
        fp = fastp_stats(args.fastp_dir, s)
        counts = read_tsv(Path(args.counts_dir) / f"{s}.counts.tsv")[0]
        trimmed = int(counts["trimmed_pairs"])
        if trimmed != fp["trimmed_pairs"]:
            raise SystemExit(
                f"{s}: Bowtie2 saw {trimmed} pairs but fastp kept "
                f"{fp['trimmed_pairs']}; stage 3 did not read stage 2's output"
            )
        retained = int(counts["retained_pairs"])
        rows.append(
            {
                "run": s,
                "raw_pairs": fp["raw_pairs"],
                "trimmed_pairs": trimmed,
                "human_pairs": int(counts["human_pairs"]),
                "retained_pairs": retained,
                "retained_pct_of_raw": round(100 * retained / fp["raw_pairs"], 2),
            }
        )
    fields = [
        "run",
        "raw_pairs",
        "trimmed_pairs",
        "human_pairs",
        "retained_pairs",
        "retained_pct_of_raw",
    ]
    write_tsv(args.out, rows, fields)


# ---- stage 4 -------------------------------------------------------------------------


def quast_report(path):
    rows = {}
    for line in open(path):
        key, _, value = line.rstrip("\n").partition("\t")
        rows[key] = value
    return rows


def bowtie2_overall_rate(log):
    for line in open(log):
        if "overall alignment rate" in line:
            return float(line.split("%")[0])
    raise SystemExit(f"no overall alignment rate in {log}")


def table_assembly(args):
    rows = []
    for s in read_samples(args.samples):
        d = Path(args.assembly_dir) / s
        q = quast_report(d / "quast" / "report.tsv")
        rows.append(
            {
                "run": s,
                "contigs": int(q["# contigs"]),
                "total_bp": int(q["Total length"]),
                "n50": int(q["N50"]),
                "assembled_fraction_pct": bowtie2_overall_rate(d / "map" / f"{s}.bowtie2.log"),
                "graph_k": (d / "graph_k.txt").read_text().strip(),
            }
        )
    write_tsv(
        args.out, rows, ["run", "contigs", "total_bp", "n50", "assembled_fraction_pct", "graph_k"]
    )


# ---- stage 5 -------------------------------------------------------------------------


def genomad_file(genomad_dir, suffix):
    hits = sorted(Path(genomad_dir).glob(f"*_summary/*_{suffix}"))
    if not hits:
        raise SystemExit(f"no *_{suffix} under {genomad_dir}")
    return hits[0]


def genomad_calls(genomad_dir):
    """plasmid contig names, virus contig names, provirus regions {contig: [(start, end)]}"""
    plasmids = {r["seq_name"] for r in read_tsv(genomad_file(genomad_dir, "plasmid_summary.tsv"))}
    viruses, proviruses = set(), {}
    for r in read_tsv(genomad_file(genomad_dir, "virus_summary.tsv")):
        name = r["seq_name"]
        if "|provirus_" in name:
            contig, _, region = name.partition("|provirus_")
            start, end = (int(x) for x in region.split("_")[:2])
            proviruses.setdefault(contig, []).append((start, end))
        else:
            viruses.add(name)
    return plasmids, viruses, proviruses


def tag_contig(contig, start, stop, plasmids, viruses, proviruses):
    if contig in viruses:
        return "virus"
    if any(start <= e and stop >= s for s, e in proviruses.get(contig, [])):
        return "virus"
    if contig in plasmids:
        return "plasmid"
    return "chromosome"


def table_tag_amr(args):
    plasmids, viruses, proviruses = genomad_calls(args.genomad_dir)
    hits = []
    for r in read_tsv(args.amrfinder):
        hits.append(
            (
                "amrfinder",
                r["Contig id"],
                int(r["Start"]),
                int(r["Stop"]),
                r["Element symbol"],
                r.get("Class", ""),
                r.get("Type", ""),
            )
        )
    for r in read_tsv(args.rgi):
        hits.append(
            (
                "rgi",
                r["Contig"],
                int(r["Start"]),
                int(r["Stop"]),
                r["Best_Hit_ARO"],
                r.get("Drug Class", ""),
                r.get("Cut_Off", ""),
            )
        )
    rows = []
    for i, (source, contig, a, b, gene, drug, detail) in enumerate(hits, 1):
        start, stop = min(a, b), max(a, b)
        rows.append(
            {
                "hit_id": f"{args.run}_{source}_{i}",
                "source": source,
                "tag": tag_contig(contig, start, stop, plasmids, viruses, proviruses),
                "run": args.run,
                "contig": contig,
                "start": start,
                "stop": stop,
                "gene": gene,
                "drug_class": drug,
                "detail": detail,
            }
        )
    write_tsv(
        args.out,
        rows,
        [
            "hit_id",
            "source",
            "tag",
            "run",
            "contig",
            "start",
            "stop",
            "gene",
            "drug_class",
            "detail",
        ],
    )
    # BED of the hits for samtools bedcov (0-based start)
    with open(args.bed, "w") as fh:
        for r in rows:
            fh.write(f"{r['contig']}\t{r['start'] - 1}\t{r['stop']}\t{r['hit_id']}\n")


def table_amr_counts(args):
    """samtools bedcov -c output: 4 BED columns, one depth sum per BAM, one read count per BAM"""
    runs = [Path(b).name.removesuffix(".bam") for b in args.bams]
    n = len(runs)
    rows = []
    for line in open(args.bedcov):
        f = line.rstrip("\n").split("\t")
        counts = f[4 + n : 4 + 2 * n]
        if len(counts) != n:
            raise SystemExit(f"unexpected bedcov line (want 4 + 2x{n} columns): {line!r}")
        rows.append({"hit_id": f[3], **{f"reads_{r}": int(c) for r, c in zip(runs, counts)}})
    write_tsv(args.out, rows, ["hit_id"] + [f"reads_{r}" for r in runs])


# ---- stage 6 / 7 ---------------------------------------------------------------------


def mimag_tier(comp, cont, a):
    if comp >= a.high_comp and cont < a.high_cont:
        return "high"
    if comp >= a.medium_comp and cont < a.medium_cont:
        return "medium"
    return "low"


def table_mimag(args):
    rows = []
    for r in read_tsv(args.report):
        comp, cont = float(r["Completeness"]), float(r["Contamination"])
        rows.append(
            {
                "run": args.run,
                "bin": r["Name"],
                "completeness": comp,
                "contamination": cont,
                "tier": mimag_tier(comp, cont, args),
            }
        )
    write_tsv(args.out, rows, ["run", "bin", "completeness", "contamination", "tier"])


def all_mimag(binning_dir, samples):
    out = {}
    for s in read_samples(samples):
        for r in read_tsv(Path(binning_dir) / s / "mimag.tsv"):
            out[r["bin"]] = r
    return out


def table_drep_info(args):
    rows = [
        {
            "genome": f"{b}.fa",
            "completeness": r["completeness"],
            "contamination": r["contamination"],
        }
        for b, r in all_mimag(args.binning_dir, args.samples).items()
        if r["tier"] in ("high", "medium")
    ]
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    with open(args.out, "w", newline="") as fh:
        w = csv.DictWriter(
            fh, fieldnames=["genome", "completeness", "contamination"], lineterminator="\n"
        )
        w.writeheader()
        w.writerows(rows)


def gtdbtk_taxonomy(mags_dir):
    tax = {}
    for f in sorted(Path(mags_dir, "gtdbtk").rglob("gtdbtk.*.summary.tsv")):
        for r in read_tsv(f):
            tax.setdefault(r["user_genome"], r["classification"])
    return tax


def drep_representatives(mags_dir):
    """genome file -> representative (winner) genome file of its secondary cluster"""
    tables = Path(mags_dir, "drep", "data_tables")
    cluster = {
        r["genome"]: r["secondary_cluster"] for r in csv.DictReader(open(tables / "Cdb.csv"))
    }
    winner = {r["cluster"]: r["genome"] for r in csv.DictReader(open(tables / "Wdb.csv"))}
    return {g: winner.get(c, "") for g, c in cluster.items()}


def run_of(label, runs):
    for r in runs:
        if label.startswith(r):
            return r
    return None


def table_mags(args):
    runs = read_samples(args.samples)
    mimag = all_mimag(args.binning_dir, args.samples)
    tax = gtdbtk_taxonomy(args.mags_dir)
    abund = {}
    for r in read_tsv(Path(args.mags_dir) / "coverm_genome.tsv"):
        genome = r.pop("Genome")
        for col, value in r.items():
            run = run_of(col, runs)
            if run is None:
                continue
            kind = "relabund_pct" if "Relative Abundance" in col else "mean_cov"
            abund.setdefault(genome, {})[f"{kind}_{run}"] = value
    rows = []
    genomes = sorted(
        p.stem for p in Path(args.mags_dir, "drep", "dereplicated_genomes").glob("*.fa")
    )
    for g in genomes:
        q = mimag.get(g, {})
        rows.append(
            {
                "mag": g,
                "source_run": run_of(g, runs),
                "completeness": q.get("completeness"),
                "contamination": q.get("contamination"),
                "tier": q.get("tier"),
                "gtdb_classification": tax.get(g, "unclassified"),
                **abund.get(g, {}),
            }
        )
    ab_cols = [f"{k}_{r}" for r in runs for k in ("relabund_pct", "mean_cov")]
    write_tsv(
        args.out,
        rows,
        ["mag", "source_run", "completeness", "contamination", "tier", "gtdb_classification"]
        + ab_cols,
    )


# ---- stage 8 -------------------------------------------------------------------------


def fasta_names(path):
    return [line[1:].split()[0] for line in open(path) if line.startswith(">")]


def coverage_profiles(coverm_tsv):
    """contig -> {reads_run: mean coverage} from CoverM contig output"""
    prof = {}
    for r in read_tsv(coverm_tsv):
        contig = r.pop("Contig")
        prof[contig] = {
            col.rsplit(" Mean", 1)[0]: float(v) for col, v in r.items() if col.endswith(" Mean")
        }
    return prof


def ranks(values):
    order = sorted(range(len(values)), key=values.__getitem__)
    out = [0.0] * len(values)
    i = 0
    while i < len(order):
        j = i
        while j + 1 < len(order) and values[order[j + 1]] == values[order[i]]:
            j += 1
        for k in range(i, j + 1):
            out[order[k]] = (i + j) / 2 + 1
        i = j + 1
    return out


def spearman(x, y):
    rx, ry = ranks(x), ranks(y)
    n = len(x)
    mx, my = sum(rx) / n, sum(ry) / n
    sxy = sum((a - mx) * (b - my) for a, b in zip(rx, ry))
    sxx = sum((a - mx) ** 2 for a in rx)
    syy = sum((b - my) ** 2 for b in ry)
    return sxy / (sxx * syy) ** 0.5 if sxx and syy else float("nan")


def table_plasmid_hosts(args):
    runs = read_samples(args.samples)
    mimag = all_mimag(args.binning_dir, args.samples)
    reps = drep_representatives(args.mags_dir)
    tax = gtdbtk_taxonomy(args.mags_dir)
    in_mags, inference = [], []
    for run in runs:
        bins_dir = Path(args.binning_dir) / run / "bins"
        contig_bin = {c: f.stem for f in bins_dir.glob("*.fa") for c in fasta_names(f)}
        passing = {
            b
            for b in set(contig_bin.values())
            if mimag.get(b, {}).get("tier") in ("high", "medium")
        }
        genomad_dir = Path(args.contigs_dir) / run / "genomad"
        plasmids = read_tsv(genomad_file(genomad_dir, "plasmid_summary.tsv"))
        mob_file = Path(args.contigs_dir) / run / "mob_typer.tsv"
        mob = (
            {r["sample_id"].split()[0]: r for r in read_tsv(mob_file)} if mob_file.exists() else {}
        )
        prof = coverage_profiles(Path(args.contigs_dir) / run / "coverm_contig.tsv")
        samples_cov = sorted(next(iter(prof.values()), {}).keys())
        bin_prof = {}
        for b in passing:
            members = [c for c, bb in contig_bin.items() if bb == b and c in prof]
            bin_prof[b] = [sum(prof[c][s] for c in members) / len(members) for s in samples_cov]

        for p in plasmids:
            name = p["seq_name"]
            b = contig_bin.get(name, "")
            rep = reps.get(f"{b}.fa", "").removesuffix(".fa") if b in passing else ""
            if b:
                in_mags.append(
                    {
                        "run": run,
                        "plasmid_contig": name,
                        "length": p["length"],
                        "bin": b,
                        "bin_tier": mimag.get(b, {}).get("tier", ""),
                        "drep_representative": rep,
                    }
                )
            m = mob.get(name, {})
            best_bin, best_rho = "", ""
            n = len(samples_cov)
            if n >= args.min_samples and name in prof:
                x = [prof[name][s] for s in samples_cov]
                scored = [(spearman(x, y), bb) for bb, y in bin_prof.items()]
                scored = [t for t in scored if t[0] == t[0]]
                if scored:
                    best_rho, best_bin = max(scored)
                    best_rho = round(best_rho, 4)
            best_rep = reps.get(f"{best_bin}.fa", "").removesuffix(".fa") if best_bin else ""
            inference.append(
                {
                    "run": run,
                    "plasmid_contig": name,
                    "length": p["length"],
                    "cobin_bin": b if b in passing else "",
                    "cobin_mag": rep,
                    "cobin_taxonomy": tax.get(rep, ""),
                    "mob_host_range_rank": m.get("predicted_host_range_overall_rank", ""),
                    "mob_host_range_name": m.get("predicted_host_range_overall_name", ""),
                    "mob_mobility": m.get("predicted_mobility", ""),
                    "corr_n_samples": n,
                    "corr_best_bin": best_bin,
                    "corr_best_mag": best_rep,
                    "corr_best_rho": best_rho,
                    "corr_best_taxonomy": tax.get(best_rep, ""),
                }
            )
    write_tsv(
        Path(args.out_dir) / "plasmids_in_mags.tsv",
        in_mags,
        ["run", "plasmid_contig", "length", "bin", "bin_tier", "drep_representative"],
    )
    write_tsv(
        Path(args.out_dir) / "plasmid_host_inference.tsv",
        inference,
        [
            "run",
            "plasmid_contig",
            "length",
            "cobin_bin",
            "cobin_mag",
            "cobin_taxonomy",
            "mob_host_range_rank",
            "mob_host_range_name",
            "mob_mobility",
            "corr_n_samples",
            "corr_best_bin",
            "corr_best_mag",
            "corr_best_rho",
            "corr_best_taxonomy",
        ],
    )


# ---- stage 9 -------------------------------------------------------------------------


def table_plasmid_summary(args):
    rows = []
    for run in read_samples(args.samples):
        d = Path(args.contigs_dir) / run
        pl = read_tsv(genomad_file(d / "genomad", "plasmid_summary.tsv"))
        vi = read_tsv(genomad_file(d / "genomad", "virus_summary.tsv"))
        mob = read_tsv(d / "mob_typer.tsv") if (d / "mob_typer.tsv").exists() else []
        mobility = [r.get("predicted_mobility", "") for r in mob]
        tagged = read_tsv(d / "amr_tagged.tsv")
        amr_plasmid_contigs = {r["contig"] for r in tagged if r["tag"] == "plasmid"}
        scapp = sorted(Path(args.assembly_dir, run, "scapp").glob("*confident_cycs.fasta"))
        rows.append(
            {
                "run": run,
                "plasmid_contigs": len(pl),
                "plasmid_bp": sum(int(r["length"]) for r in pl),
                "plasmid_dtr_circular": sum(r.get("topology") == "DTR" for r in pl),
                "plasmid_with_conjugation_genes": sum(
                    r.get("conjugation_genes", "NA") not in ("NA", "") for r in pl
                ),
                "plasmid_with_amr_hits": len(amr_plasmid_contigs),
                "mob_conjugative": mobility.count("conjugative"),
                "mob_mobilizable": mobility.count("mobilizable"),
                "mob_non_mobilizable": mobility.count("non-mobilizable"),
                "virus_contigs": len(vi),
                "scapp_plasmids": len(fasta_names(scapp[0])) if scapp else "NA",
            }
        )
    write_tsv(args.out, rows, list(rows[0].keys()) if rows else ["run"])


def table_resistome(args):
    """gene x sample read counts over AMR hits (own reads, own assembly), split by tag"""
    runs = read_samples(args.samples)
    matrices = {"all": {}, "mobile": {}, "chromosome": {}}
    for run in runs:
        d = Path(args.contigs_dir) / run
        counts = {r["hit_id"]: r for r in read_tsv(d / "amr_counts.tsv")}
        for h in read_tsv(d / "amr_tagged.tsv"):
            if h["source"] != args.source:
                continue
            n = int(counts[h["hit_id"]][f"reads_{run}"])
            part = "mobile" if h["tag"] in ("plasmid", "virus") else "chromosome"
            for m in ("all", part):
                matrices[m].setdefault(h["gene"], {}).setdefault(run, 0)
                matrices[m][h["gene"]][run] += n
    for name, mat in matrices.items():
        rows = [{"gene": g, **{r: v.get(r, 0) for r in runs}} for g, v in sorted(mat.items())]
        write_tsv(Path(args.out_dir) / f"resistome_{args.source}_{name}.tsv", rows, ["gene"] + runs)


def table_rgi_compare(args):
    rows = []
    for run in read_samples(args.samples):
        d = Path(args.contigs_dir) / run
        counts = {r["hit_id"]: r for r in read_tsv(d / "amr_counts.tsv")}
        contig = {}
        for h in read_tsv(d / "amr_tagged.tsv"):
            if h["source"] == "rgi":
                c = contig.setdefault(h["gene"], [0, 0])
                c[0] += 1
                c[1] += int(counts[h["hit_id"]][f"reads_{run}"])
        reads = {}
        for r in read_tsv(Path(args.bwt_dir) / run / f"{run}.gene_mapping_data.txt"):
            reads[r["ARO Term"]] = reads.get(r["ARO Term"], 0) + float(r["All Mapped Reads"])
        for aro in sorted(set(contig) | set(reads)):
            rows.append(
                {
                    "run": run,
                    "aro_term": aro,
                    "contig_level_hits": contig.get(aro, [0, 0])[0],
                    "contig_level_reads": contig.get(aro, [0, 0])[1],
                    "read_level_reads": reads.get(aro, 0),
                    "found_by": "both"
                    if aro in contig and aro in reads
                    else ("contig_only" if aro in contig else "read_only"),
                }
            )
    write_tsv(
        args.out,
        rows,
        [
            "run",
            "aro_term",
            "contig_level_hits",
            "contig_level_reads",
            "read_level_reads",
            "found_by",
        ],
    )


def table_depth(args):
    rows = []
    for d in sorted(Path(args.depth_dir).glob("*/yield.tsv")):
        rows.extend(read_tsv(d))
    write_tsv(
        args.out,
        rows,
        ["run", "level", "pairs", "contigs", "assembly_bp", "plasmid_contigs", "plasmid_bp"],
    )


def table_depth_yield(args):
    pl = read_tsv(genomad_file(args.genomad_dir, "plasmid_summary.tsv"))
    contigs = fasta_names(args.contigs)
    bp = sum(len(line.strip()) for line in open(args.contigs) if not line.startswith(">"))
    write_tsv(
        args.out,
        [
            {
                "run": args.run,
                "level": args.level,
                "pairs": args.pairs,
                "contigs": len(contigs),
                "assembly_bp": bp,
                "plasmid_contigs": len(pl),
                "plasmid_bp": sum(int(r["length"]) for r in pl),
            }
        ],
        ["run", "level", "pairs", "contigs", "assembly_bp", "plasmid_contigs", "plasmid_bp"],
    )


ANTIBIOTICS = [
    "Ertapenem",
    "Amoxicillin",
    "Ciprofloxacin",
    "Doxycycline",
    "Azithromycin",
    "Clindamycin",
    "Sulfamethoxazole",
    "Cephalexin",
    "Trimethoprim",
    "Levofloxacin",
]


def fetch_biosamples(accessions):
    import urllib.request
    import xml.etree.ElementTree as ET

    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=biosample&id=" + ",".join(
        accessions
    )
    root = ET.fromstring(urllib.request.urlopen(url, timeout=120).read())
    out = {}
    for bs in root.iter("BioSample"):
        attrs = {a.get("attribute_name"): (a.text or "").strip() for a in bs.iter("Attribute")}
        out[bs.get("accession")] = attrs
    return out


def encode_censored(raw, how, half_min):
    if raw == "ND" or raw == "<LOQ":
        return {"na": "NA", "zero": 0, "half_min": half_min}[how]
    return float(raw)


def table_stats_samples(args):
    sheet = {r["run"]: r for r in read_tsv(args.sheet)}
    runs = read_samples(args.samples)
    groups = {}
    for spec in args.site_groups:
        group, _, sites = spec.partition(":")
        for site in sites.split(","):
            groups[site] = group
    attrs = fetch_biosamples([sheet[r]["biosample"] for r in runs])
    raw = {
        r: {ab: attrs[sheet[r]["biosample"]].get(f"{ab}(ng/L)", "NA") for ab in ANTIBIOTICS}
        for r in runs
    }
    half_min = {}
    for ab in ANTIBIOTICS:
        detected = [float(v[ab]) for v in raw.values() if v[ab] not in ("ND", "<LOQ", "NA")]
        half_min[ab] = min(detected) / 2 if detected else "NA"
    rows = []
    for r in runs:
        row = {k: sheet[r][k] for k in ("run", "biosample", "location", "site", "timepoint")}
        row["site_group"] = groups.get(sheet[r]["site"], "other")
        for ab in ANTIBIOTICS:
            v = raw[r][ab]
            row[f"{ab.lower()}_raw"] = v
            how = args.nd_encoding if v == "ND" else args.loq_encoding
            row[f"{ab.lower()}_ngl"] = "NA" if v == "NA" else encode_censored(v, how, half_min[ab])
        rows.append(row)
    fields = ["run", "biosample", "location", "site", "timepoint", "site_group"]
    fields += [f"{ab.lower()}_{k}" for ab in ANTIBIOTICS for k in ("raw", "ngl")]
    write_tsv(args.out, rows, fields)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    sub = parser.add_subparsers(dest="table", required=True)

    p = sub.add_parser("qc", help="stage 2: raw/trimmed pairs, Q30, duplication, insert size")
    p.add_argument("--samples", required=True)
    p.add_argument("--fastp-dir", required=True)
    p.add_argument("--out", required=True)
    p.set_defaults(func=table_qc)

    p = sub.add_parser("host", help="stage 3: raw, trimmed, human, retained pairs")
    p.add_argument("--samples", required=True)
    p.add_argument("--fastp-dir", required=True)
    p.add_argument("--counts-dir", required=True)
    p.add_argument("--out", required=True)
    p.set_defaults(func=table_host)

    p = sub.add_parser("assembly", help="stage 4: contigs, total bp, N50, assembled fraction")
    p.add_argument("--samples", required=True)
    p.add_argument("--assembly-dir", required=True)
    p.add_argument("--out", required=True)
    p.set_defaults(func=table_assembly)

    p = sub.add_parser("tag-amr", help="stage 5: tag AMR hits plasmid/virus/chromosome")
    p.add_argument("--run", required=True)
    p.add_argument("--genomad-dir", required=True)
    p.add_argument("--amrfinder", required=True)
    p.add_argument("--rgi", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--bed", required=True)
    p.set_defaults(func=table_tag_amr)

    p = sub.add_parser("amr-counts", help="stage 5: reads over each AMR hit, per BAM")
    p.add_argument("--bedcov", required=True)
    p.add_argument("--bams", nargs="+", required=True)
    p.add_argument("--out", required=True)
    p.set_defaults(func=table_amr_counts)

    p = sub.add_parser("mimag", help="stage 6: MIMAG tiers from CheckM2")
    p.add_argument("--report", required=True)
    p.add_argument("--run", required=True)
    for k in ("high-comp", "high-cont", "medium-comp", "medium-cont"):
        p.add_argument(f"--{k}", type=float, required=True)
    p.add_argument("--out", required=True)
    p.set_defaults(func=table_mimag)

    p = sub.add_parser("drep-info", help="stage 7: dRep --genomeInfo from MIMAG tables")
    p.add_argument("--binning-dir", required=True)
    p.add_argument("--samples", required=True)
    p.add_argument("--out", required=True)
    p.set_defaults(func=table_drep_info)

    p = sub.add_parser("mags", help="stage 7: quality, taxonomy, abundance per MAG")
    p.add_argument("--samples", required=True)
    p.add_argument("--binning-dir", required=True)
    p.add_argument("--mags-dir", required=True)
    p.add_argument("--out", required=True)
    p.set_defaults(func=table_mags)

    p = sub.add_parser("plasmid-hosts", help="stage 8: plasmids in MAGs + host inference")
    p.add_argument("--samples", required=True)
    p.add_argument("--binning-dir", required=True)
    p.add_argument("--mags-dir", required=True)
    p.add_argument("--contigs-dir", required=True)
    p.add_argument("--min-samples", type=int, required=True)
    p.add_argument("--out-dir", required=True)
    p.set_defaults(func=table_plasmid_hosts)

    p = sub.add_parser("plasmid-summary", help="stage 9: plasmid summary per sample")
    p.add_argument("--samples", required=True)
    p.add_argument("--contigs-dir", required=True)
    p.add_argument("--assembly-dir", required=True)
    p.add_argument("--out", required=True)
    p.set_defaults(func=table_plasmid_summary)

    p = sub.add_parser("resistome", help="stage 9: gene x sample matrix, mobile subset separate")
    p.add_argument("--samples", required=True)
    p.add_argument("--contigs-dir", required=True)
    p.add_argument("--source", choices=["amrfinder", "rgi"], required=True)
    p.add_argument("--out-dir", required=True)
    p.set_defaults(func=table_resistome)

    p = sub.add_parser("rgi-compare", help="stage 9: RGI read level (bwt) vs contig level")
    p.add_argument("--samples", required=True)
    p.add_argument("--contigs-dir", required=True)
    p.add_argument("--bwt-dir", required=True)
    p.add_argument("--out", required=True)
    p.set_defaults(func=table_rgi_compare)

    p = sub.add_parser("depth-yield", help="stage 9: one depth-experiment point")
    p.add_argument("--run", required=True)
    p.add_argument("--level", required=True)
    p.add_argument("--pairs", required=True)
    p.add_argument("--contigs", required=True)
    p.add_argument("--genomad-dir", required=True)
    p.add_argument("--out", required=True)
    p.set_defaults(func=table_depth_yield)

    p = sub.add_parser("depth", help="stage 9: collect depth-experiment points")
    p.add_argument("--depth-dir", required=True)
    p.add_argument("--out", required=True)
    p.set_defaults(func=table_depth)

    p = sub.add_parser("stats-samples", help="stage 9: sample metadata + BioSample antibiotics")
    p.add_argument("--samples", required=True)
    p.add_argument("--sheet", required=True)
    p.add_argument("--site-groups", nargs="+", required=True, help="group:SITE,SITE ...")
    p.add_argument("--nd-encoding", choices=["na", "zero", "half_min"], required=True)
    p.add_argument("--loq-encoding", choices=["na", "zero", "half_min"], required=True)
    p.add_argument("--out", required=True)
    p.set_defaults(func=table_stats_samples)

    args = parser.parse_args(argv)
    args.func(args)
    return 0


if __name__ == "__main__":
    sys.exit(main())
