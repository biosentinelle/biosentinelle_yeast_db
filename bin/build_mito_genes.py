#!/usr/bin/env python3
"""
Build the list of mitochondrial genes needed for respiration.

Growth on glycerol, or any non-fermentable carbon source, needs a working
respiratory chain, hence a mitochondrial genome (mtDNA) that still encodes it:

    rho+   complete mtDNA                          grows on glycerol
    rho-   mtDNA reduced to a fragment, often      does not grow on glycerol
           amplified many times
    rho0   no mtDNA at all                         does not grow on glycerol

The sample test (biosentinelle bin/mito_status.py) reads the depth of the reads
over the genes written here. Only the exons are kept: the group I and II
introns of COX1, COB and the 21S rRNA are optional, and strains differ in
which ones they carry (S288C has 13, many isolates fewer). A missing intron
leaves a hole in the coverage that is not a rho- deletion; a missing exon is.
The intron-encoded ORFs (AI1 to AI5, BI2 to BI4, SCEI) and the ORFs of unknown
function (Q0010, Q0255, ...) are left out for the same reason.

Kept: the eight protein genes of the respiratory complexes and the mito-
ribosome (COX1, COX2, COX3, COB, ATP6, ATP8, OLI1 = ATP9, VAR1), the two rRNAs,
the RNase P RNA (RPM1) and every tRNA. Losing any of them abolishes
respiration.

Coordinates are fetched from the Ensembl REST API, like the marker loci, and
the assembly must be the one the reference FASTA is.
"""

import argparse
import csv
import os
import sys

from build_marker_loci import ENSEMBL, chromosome_names, fetch_json

PROTEINS = ["COX1", "COX2", "COX3", "COB", "ATP6", "ATP8", "OLI1", "VAR1"]
RNAS = ["15S_RRNA", "21S_RRNA", "RPM1"]


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--reference", required=True, help="reference genome FASTA")
    parser.add_argument("--species", default="saccharomyces_cerevisiae")
    parser.add_argument("--assembly", default="R64-1-1",
                        help="assembly name the reference FASTA is; Ensembl must agree")
    parser.add_argument("--region", default="Mito",
                        help="Ensembl name of the mitochondrial chromosome (default Mito)")
    parser.add_argument("--out", required=True, help="mito_genes.tsv")
    args = parser.parse_args()

    names, sequences = chromosome_names(args.reference)
    chrom = names.get(args.region)
    if chrom is None:
        sys.exit("no mitochondrial chromosome in %s: a header must say "
                 "[location=mitochondrion]" % args.reference)
    length = len(sequences[chrom])
    release = fetch_json(ENSEMBL + "/info/data?content-type=application/json")["releases"]
    genes = fetch_json("%s/overlap/region/%s/%s:1-%d?feature=gene;content-type=application/json"
                       % (ENSEMBL, args.species, args.region, length))

    rows, found = [], set()
    for gene in sorted(genes, key=lambda g: g["start"]):
        name = gene.get("external_name") or gene["id"]
        if gene["biotype"] == "tRNA":
            role = "tRNA"
        elif name in PROTEINS:
            role = "respiratory protein" if name != "VAR1" else "mitoribosome protein"
        elif name in RNAS:
            role = "rRNA" if "RRNA" in name else "RNase P RNA"
        else:
            continue
        info = fetch_json("%s/lookup/id/%s?expand=1;content-type=application/json"
                          % (ENSEMBL, gene["id"]))
        if info.get("assembly_name") != args.assembly:
            sys.exit("%s: Ensembl assembly %s, expected %s"
                     % (name, info.get("assembly_name"), args.assembly))
        transcripts = info.get("Transcript", [])
        if len(transcripts) != 1:
            sys.exit("%s: %d transcripts, one expected" % (name, len(transcripts)))
        exons = sorted((e["start"] - 1, e["end"]) for e in transcripts[0].get("Exon", []))
        if not exons:
            exons = [(info["start"] - 1, info["end"])]
        found.add(name)
        rows.append({
            "gene_id": gene["id"], "gene": name, "role": role, "chrom": chrom,
            "start": info["start"] - 1, "end": info["end"],
            "strand": "+" if info["strand"] == 1 else "-",
            "exons": ";".join("%d-%d" % e for e in exons),
            "exon_bp": sum(b - a for a, b in exons),
            "source": "%s/lookup/id/%s (Ensembl release %s, %s)"
                      % (ENSEMBL, gene["id"], ",".join(map(str, release)),
                         info.get("assembly_name")),
        })
        sys.stderr.write("%-10s %-20s %s %d exon(s), %d bp\n"
                         % (name, role, chrom, len(exons), rows[-1]["exon_bp"]))

    missing = [g for g in PROTEINS + RNAS if g not in found]
    if missing:
        sys.exit("not found on %s: %s" % (args.region, ", ".join(missing)))
    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    with open(args.out, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t",
                                lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    sys.stderr.write("%d genes written to %s (%d tRNAs)\n"
                     % (len(rows), args.out, sum(r["role"] == "tRNA" for r in rows)))
    return 0


if __name__ == "__main__":
    sys.exit(main())
