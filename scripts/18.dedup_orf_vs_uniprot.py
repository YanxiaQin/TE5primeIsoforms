#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
from collections import defaultdict, OrderedDict

def read_fasta(path: str):
    header = None
    seq_chunks = []
    with open(path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq_chunks).replace(" ", "").upper()
                header = line[1:].strip()
                seq_chunks = []
            else:
                seq_chunks.append(line.strip())
        if header is not None:
            yield header, "".join(seq_chunks).replace(" ", "").upper()

def wrap_seq(seq: str, width: int = 60) -> str:
    return "\n".join(seq[i:i+width] for i in range(0, len(seq), width))

def main():
    ap = argparse.ArgumentParser(description="Dedup ORF fasta by sequence, remove any ORF sequences that appear in UniProt reviewed fasta.")
    ap.add_argument("--uniprot_fa", required=True, help="UniProt reviewed FASTA (human)")
    ap.add_argument("--orf_fa", required=True, help="DefaultPrimary.ORF.fa")
    ap.add_argument("--out_fa", default="DefaultPrimary.ORF.unique_only.fa", help="Output FASTA (unique ORF sequences not in UniProt)")
    ap.add_argument("--out_tsv", default="DefaultPrimary.ORF.unique_only.dedup_report.tsv", help="Report TSV")
    ap.add_argument("--wrap", type=int, default=60, help="FASTA sequence line width (default 60)")
    args = ap.parse_args()

    uniprot_seq_to_headers = defaultdict(list)
    uniprot_seqs = set()

    n_u = 0
    for h, s in read_fasta(args.uniprot_fa):
        n_u += 1
        uniprot_seqs.add(s)
        uniprot_seq_to_headers[s].append(h)

    orf_seq_to_headers = defaultdict(list)
    orf_records = []
    n_o = 0
    for h, s in read_fasta(args.orf_fa):
        n_o += 1
        orf_records.append((h, s))
        orf_seq_to_headers[s].append(h)

    kept = OrderedDict()
    for seq, headers in orf_seq_to_headers.items():
        if seq in uniprot_seqs:
            continue
        kept[seq] = headers[0]

    with open(args.out_fa, "w", encoding="utf-8") as outfa:
        for seq, rep_h in kept.items():
            outfa.write(f">{rep_h}\n")
            outfa.write(wrap_seq(seq, width=args.wrap) + "\n")

    with open(args.out_tsv, "w", encoding="utf-8") as outtsv:
        outtsv.write("\t".join([
            "status",
            "orf_header",
            "orf_seq_len",
            "rep_header_for_seq",
            "n_headers_with_same_seq_in_orf",
            "all_orf_headers_with_same_seq",
            "matched_uniprot",
            "uniprot_headers_matched"
        ]) + "\n")

        rep_for_seq = {seq: hdr for seq, hdr in kept.items()}

        for h, s in orf_records:
            same_headers = orf_seq_to_headers[s]
            n_same = len(same_headers)

            if s in uniprot_seqs:
                status = "DROP_MATCH_UNIPROT"
                rep_h = ""
                matched_uniprot = "YES"
                uni_h = "|".join(uniprot_seq_to_headers[s][:20])  # avoid insane length
                if len(uniprot_seq_to_headers[s]) > 20:
                    uni_h += f"|...(+{len(uniprot_seq_to_headers[s]) - 20} more)"
            else:
                matched_uniprot = "NO"
                uni_h = ""
                rep_h = rep_for_seq.get(s, "")
                if rep_h and h == rep_h:
                    status = "KEPT"
                else:
                    status = "DROP_DUP_IN_ORF"

            outtsv.write("\t".join([
                status,
                h,
                str(len(s)),
                rep_h,
                str(n_same),
                "|".join(same_headers),
                matched_uniprot,
                uni_h
            ]) + "\n")

    print(f"[OK] UniProt records: {n_u}")
    print(f"[OK] ORF records: {n_o}")
    print(f"[OK] Kept unique ORF sequences (not in UniProt): {len(kept)}")
    print(f"[OUT] {args.out_fa}")
    print(f"[OUT] {args.out_tsv}")

if __name__ == "__main__":
    main()
