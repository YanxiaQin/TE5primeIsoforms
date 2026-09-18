from pathlib import Path
import re
import pandas as pd

# =========================================================
# 0. 路径
# =========================================================
PACBIO_ROOT = Path("/home/yxqin/neoantigen/Analysis/pacbio/05_repeatmasker")
MERGE_XLSX = Path("/home/yxqin/neoantigen/Analysis/pacbio/09_common_epitopes/deduplication_isoform/get_transcript_and_orf_seq/star_validation/star_support_out/get_candidate_junctions/merge_Initial_mainRoleTE.xlsx")
OUT_TSV = Path("/home/yxqin/neoantigen/Analysis/pacbio/09_common_epitopes/deduplication_isoform/get_transcript_and_orf_seq/star_validation/star_support_out/get_candidate_junctions/candidate_junctions.tsv")

MERGE_SHEET = 0


# =========================================================
# 1. 读取 merge_Initial_mainRoleTE.xlsx
# =========================================================
merge_df = pd.read_excel(MERGE_XLSX, sheet_name=MERGE_SHEET, dtype=str)

merge_df.columns = [str(c).strip() for c in merge_df.columns]

required_cols = {"SRR_ID", "TE_Role", "transcript_id", "gene_id", "tx_exon_count", "repeat_name"}
missing = required_cols - set(merge_df.columns)
if missing:
    raise ValueError(f"merge_Initial_mainRoleTE.xlsx 缺少必要列: {missing}")

merge_df = merge_df[list(required_cols)].copy()
merge_df = merge_df.dropna(subset=["SRR_ID", "transcript_id"])
merge_df["SRR_ID"] = merge_df["SRR_ID"].astype(str).str.strip()
merge_df["transcript_id"] = merge_df["transcript_id"].astype(str).str.strip()
merge_df["gene_id"] = merge_df["gene_id"].astype(str).str.strip()
merge_df["TE_Role"] = merge_df["TE_Role"].astype(str).str.strip()
merge_df["repeat_name"] = merge_df["repeat_name"].astype(str).str.strip()
merge_df["tx_exon_count"] = pd.to_numeric(merge_df["tx_exon_count"], errors="coerce")

merge_df = merge_df.drop_duplicates(subset=["SRR_ID", "transcript_id"]).reset_index(drop=True)


# =========================================================
# 2. 找每个 SRR 对应的 final.gtf
# =========================================================
def find_gtf_for_srr(pacbio_root: Path, srr_id: str) -> Path:
    candidates = list(pacbio_root.glob(f"PRJNA*/{srr_id}/{srr_id}_final.gtf"))
    if len(candidates) == 0:
        raise FileNotFoundError(f"找不到 {srr_id}_final.gtf")
    if len(candidates) > 1:
        print(f"Warning: {srr_id} 找到多个 gtf，默认使用第一个: {candidates[0]}")
    return candidates[0]


# =========================================================
# 3. 解析 GTF attribute
# =========================================================
def parse_gtf_attr(attr_text):
    attrs = {}
    for m in re.finditer(r'(\S+)\s+"([^"]+)"', attr_text):
        attrs[m.group(1)] = m.group(2)
    return attrs


# =========================================================
# 4. 读取某个 SRR 的目标 transcript exon
# =========================================================
def load_target_exons_from_gtf(gtf_file: Path, target_tx_ids: set):
    rows = []

    with open(gtf_file, "r") as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9:
                continue

            chrom, source, feature, start, end, score, strand, frame, attrs = fields
            if feature != "exon":
                continue

            attr_dict = parse_gtf_attr(attrs)
            tx_id = attr_dict.get("transcript_id")
            gene_id = attr_dict.get("gene_id")

            if tx_id not in target_tx_ids:
                continue

            rows.append({
                "chr": chrom,
                "start": int(start),
                "end": int(end),
                "strand": strand,
                "transcript_id": tx_id,
                "gene_id_gtf": gene_id
            })

    return pd.DataFrame(rows)


# =========================================================
# 5. 按 transcript 生成 junction
# =========================================================
def build_junctions_for_transcript(exon_df, pacbio_srr, te_role, repeat_name, gene_id_input):
    """
    exon_df: 单个 transcript 的 exon 子表
    """
    if exon_df.empty:
        return []

    tx_id = exon_df["transcript_id"].iloc[0]
    chrom_set = set(exon_df["chr"])
    strand_set = set(exon_df["strand"])

    if len(chrom_set) != 1 or len(strand_set) != 1:
        raise ValueError(f"{pacbio_srr} {tx_id} 的 exon 跨多个 chr/strand，不合理")

    chrom = exon_df["chr"].iloc[0]
    strand = exon_df["strand"].iloc[0]

    exons_genomic = exon_df.sort_values(["start", "end"]).reset_index(drop=True).copy()
    exons_genomic["genomic_exon_idx"] = range(1, len(exons_genomic) + 1)

    if strand == "+":
        exons_tx = exons_genomic.copy()
        exons_tx["tx_exon_number"] = range(1, len(exons_tx) + 1)
    elif strand == "-":
        exons_tx = exons_genomic.copy()
        exons_tx["tx_exon_number"] = list(range(len(exons_tx), 0, -1))
    else:
        raise ValueError(f"{pacbio_srr} {tx_id} strand 非法: {strand}")

    exons_genomic = exons_tx.sort_values(["start", "end"]).reset_index(drop=True)

    junction_rows = []

    for i in range(len(exons_genomic) - 1):
        left = exons_genomic.iloc[i]
        right = exons_genomic.iloc[i + 1]

        intron_start = int(left["end"]) + 1
        intron_end = int(right["start"]) - 1

        if intron_start > intron_end:
            continue

        if strand == "+":
            tx_exon_from = int(left["tx_exon_number"])
            tx_exon_to = int(right["tx_exon_number"])
        else:
            tx_exon_from = int(right["tx_exon_number"])
            tx_exon_to = int(left["tx_exon_number"])

        is_first_junction = ({tx_exon_from, tx_exon_to} == {1, 2})
        junction_type = "TE_first_exon_to_exon2" if is_first_junction else "internal_junction"

        j_rank = min(tx_exon_from, tx_exon_to)

        candidate_id = f"{pacbio_srr}|{tx_id}|j{j_rank}"

        junction_rows.append({
            "candidate_id": candidate_id,
            "pacbio_srr": pacbio_srr,
            "transcript_id": tx_id,
            "gene_id": gene_id_input,
            "repeat_name": repeat_name,
            "te_role": te_role,

            "chr": chrom,
            "strand": strand,

            "tx_exon_from": tx_exon_from,
            "tx_exon_to": tx_exon_to,

            "left_exon_start": int(left["start"]),
            "left_exon_end": int(left["end"]),
            "right_exon_start": int(right["start"]),
            "right_exon_end": int(right["end"]),

            "intron_start": intron_start,
            "intron_end": intron_end,

            "junction_type": junction_type,
            "is_first_junction": is_first_junction
        })

    return junction_rows


# =========================================================
# 6. 主循环：按 SRR 分批解析 GTF
# =========================================================
all_junctions = []

for srr_id, sub_df in merge_df.groupby("SRR_ID"):
    print(f"Processing {srr_id} ...")

    gtf_file = find_gtf_for_srr(PACBIO_ROOT, srr_id)
    target_txs = set(sub_df["transcript_id"])

    exon_df = load_target_exons_from_gtf(gtf_file, target_txs)

    if exon_df.empty:
        print(f"Warning: {srr_id} 在 {gtf_file.name} 中没有找到任何目标 transcript exon")
        continue

    for _, row in sub_df.iterrows():
        tx_id = row["transcript_id"]
        tx_exons = exon_df[exon_df["transcript_id"] == tx_id].copy()

        if tx_exons.empty:
            print(f"Warning: {srr_id} {tx_id} 未在 GTF 中找到")
            continue

        junction_rows = build_junctions_for_transcript(
            exon_df=tx_exons,
            pacbio_srr=srr_id,
            te_role=row["TE_Role"],
            repeat_name=row["repeat_name"],
            gene_id_input=row["gene_id"]
        )
        all_junctions.extend(junction_rows)

cand_df = pd.DataFrame(all_junctions)

if cand_df.empty:
    raise ValueError("没有生成任何 junction，请检查 merge_Initial_mainRoleTE.xlsx 和 GTF 是否匹配")

cand_df = cand_df.sort_values(
    ["pacbio_srr", "transcript_id", "tx_exon_from", "tx_exon_to"]
).reset_index(drop=True)

cand_df.to_csv(OUT_TSV, sep="\t", index=False)

print("\nDone.")
print(f"candidate_junctions.tsv saved to: {OUT_TSV}")
print(cand_df.head())