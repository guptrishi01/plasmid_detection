#!/usr/bin/env python3
"""Fetch the NCBI SRA run table for an identifier and keep one library strategy.

fetch_runinfo.py <identifier> <library_strategy> <out.csv> [runs]
  identifier: BioProject (PRJNA...), SRA study (SRP...), or comma-separated run accessions
  runs:       optional comma-separated run list or file (one per line) to keep
Uses NCBI E-utilities (esearch with history -> efetch rettype=runinfo), stdlib only.
"""

import csv
import io
import re
import sys
import time
import urllib.parse
import urllib.request

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"


def get(endpoint, **params):
    url = EUTILS + endpoint + "?" + urllib.parse.urlencode(params)
    for attempt in range(5):
        try:
            with urllib.request.urlopen(url, timeout=120) as r:
                return r.read().decode()
        except OSError:
            time.sleep(2 ** attempt)
    raise SystemExit(f"NCBI request failed: {url}")


def term_for(identifier):
    ids = [x.strip() for x in identifier.split(",") if x.strip()]
    if all(re.fullmatch(r"[SED]RR\d+", x) for x in ids):
        return " OR ".join(f"{x}[accn]" for x in ids)
    if len(ids) != 1:
        raise SystemExit(f"give one BioProject/study, or run accessions only: {identifier}")
    return ids[0]


def read_runs(spec):
    if not spec:
        return None
    try:
        with open(spec) as fh:
            return {line.strip() for line in fh if line.strip()}
    except OSError:
        return {x.strip() for x in spec.split(",") if x.strip()}


def main(identifier, strategy, out, runs=None):
    search = get("esearch.fcgi", db="sra", term=term_for(identifier), usehistory="y", retmax=0)
    count = int(re.search(r"<Count>(\d+)</Count>", search).group(1))
    webenv = re.search(r"<WebEnv>([^<]+)</WebEnv>", search).group(1)
    key = re.search(r"<QueryKey>(\d+)</QueryKey>", search).group(1)
    if count == 0:
        raise SystemExit(f"no SRA records for {identifier}")
    rows, header = [], None
    for start in range(0, count, 500):
        text = get("efetch.fcgi", db="sra", rettype="runinfo", retmode="text",
                   WebEnv=webenv, query_key=key, retstart=start, retmax=500)
        reader = csv.DictReader(io.StringIO(text))
        header = header or reader.fieldnames
        rows += [r for r in reader if r.get("Run")]
        time.sleep(0.4)  # NCBI: <= 3 requests/s without an API key
    keep = read_runs(runs)
    rows = [r for r in rows if r["LibraryStrategy"] == strategy and (keep is None or r["Run"] in keep)]
    if keep is not None and len(rows) != len(keep):
        missing = keep - {r["Run"] for r in rows}
        raise SystemExit(f"runs not found under {identifier} as {strategy}: {sorted(missing)}")
    if not rows:
        raise SystemExit(f"no {strategy} runs for {identifier}")
    rows.sort(key=lambda r: r["Run"])
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=header, lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    print(f"{identifier}: {count} SRA records, {len(rows)} {strategy} runs kept", file=sys.stderr)


if __name__ == "__main__":
    main(*sys.argv[1:])
