#!/usr/bin/env python3
"""Sample sheet from an SRA run table: one row per run with location, site and timepoint
parsed from SampleName by a regex with named groups (location, site, timepoint).

sample_sheet.py <runinfo.csv> <regex> <out.tsv>
Runs whose SampleName doesn't match get location=all, site=NA, timepoint=NA (one group),
so any study still runs; only within-system cross-mapping and site statistics need it.
Columns match metadata/srp121672_wgs_samples.tsv.
"""

import csv
import re
import sys


def main(runinfo, regex, out):
    pat = re.compile(regex)
    rows = []
    for r in csv.DictReader(open(runinfo)):
        m = pat.fullmatch(r.get("SampleName", ""))
        g = m.groupdict() if m else {}
        rows.append({
            "run": r["Run"], "biosample": r["BioSample"], "sample_name": r.get("SampleName", ""),
            "library_name": r["LibraryName"], "spots": r["spots"],
            "location": g.get("location") or "all", "site": g.get("site") or "NA",
            "timepoint": g.get("timepoint") or "NA",
        })
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    unmatched = sum(r["location"] == "all" for r in rows)
    print(f"{len(rows)} runs; {unmatched} SampleName(s) did not match the regex", file=sys.stderr)


if __name__ == "__main__":
    main(*sys.argv[1:])
