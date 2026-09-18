from pathlib import Path
import pandas as pd

# =========================================================
# 0. 路径
# =========================================================
CANDIDATE_TSV = Path("/home/yxqin/neoantigen/Analysis/pacbio/09_common_epitopes/deduplication_isoform/get_transcript_and_orf_seq/star_validation/star_support_out/get_candidate_junctions/candidate_junctions.tsv")
STAR_RESULTS_DIR = Path("/home/yxqin/neoantigen/Analysis/pacbio/09_common_epitopes/deduplication_isoform/get_transcript_and_orf_seq/star_validation/star_twopass_mapping/results")

DETAIL_OUT = STAR_RESULTS_DIR / "candidate_junction_support.by_sample.tsv"
SUMMARY_OUT = STAR_RESULTS_DIR / "candidate_junction_support.summary.tsv"

SAMPLE_META_FILE = None


# =========================================================
# 1. 读 candidate_junctions.tsv
# =========================================================
cand = pd.read_csv(CANDIDATE_TSV, sep="\t", dtype={"chr": str})

required_cols = {
    "candidate_id", "pacbio_srr", "transcript_id", "gene_id", "repeat_name", "te_role",
    "chr", "strand", "intron_start", "intron_end", "junction_type", "is_first_junction"
}
missing = required_cols - set(cand.columns)
if missing:
    raise ValueError(f"candidate_junctions.tsv 缺少必要列: {missing}")

strand_map = {"+": 1, "-": 2}
cand["strand_code"] = cand["strand"].map(strand_map)
if cand["strand_code"].isna().any():
    bad = cand.loc[cand["strand_code"].isna(), ["candidate_id", "strand"]]
    raise ValueError(f"发现无法识别的strand:\n{bad}")

cand["intron_start"] = cand["intron_start"].astype(int)
cand["intron_end"] = cand["intron_end"].astype(int)


# =========================================================
# 2. 读取所有 SJ.out.tab
# =========================================================
sj_files = sorted(STAR_RESULTS_DIR.glob("*/*.SJ.out.tab"))
if not sj_files:
    raise ValueError(f"在 {STAR_RESULTS_DIR} 下没有找到任何 SJ.out.tab")

colnames = [
    "chr",
    "intron_start",
    "intron_end",
    "strand_code",
    "motif",
    "annotated",
    "unique_count",
    "multi_count",
    "max_overhang"
]

all_sj = []
for sj_file in sj_files:
    sample = sj_file.parent.name
    df = pd.read_csv(sj_file, sep="\t", header=None, names=colnames, dtype={"chr": str})
    df["sample"] = sample
    df["total_count"] = df["unique_count"] + df["multi_count"]
    all_sj.append(df)

sj = pd.concat(all_sj, ignore_index=True)

motif_map = {
    0: "non_canonical",
    1: "GT/AG",
    2: "CT/AC",
    3: "GC/AG",
    4: "CT/GC",
    5: "AT/AC",
    6: "GT/AT"
}
sj["motif_name"] = sj["motif"].map(motif_map).fillna("unknown")


# =========================================================
# 3. candidate junction × STAR SJ 匹配
# =========================================================
merged = cand.merge(
    sj,
    on=["chr", "intron_start", "intron_end", "strand_code"],
    how="left"
)

for col in ["annotated", "unique_count", "multi_count", "total_count", "max_overhang"]:
    merged[col] = merged[col].fillna(0)

merged["annotated"] = merged["annotated"].astype(int)
merged["unique_count"] = merged["unique_count"].astype(int)
merged["multi_count"] = merged["multi_count"].astype(int)
merged["total_count"] = merged["total_count"].astype(int)
merged["max_overhang"] = merged["max_overhang"].astype(int)
merged["motif_name"] = merged["motif_name"].fillna("not_detected")

merged["support_loose"] = merged["total_count"] >= 1
merged["support_moderate"] = (merged["unique_count"] >= 2) & (merged["max_overhang"] >= 8)
merged["support_strict"] = (merged["unique_count"] >= 3) & (merged["max_overhang"] >= 10)

if SAMPLE_META_FILE is not None:
    meta = pd.read_csv(SAMPLE_META_FILE, sep="\t", dtype=str)
    if "sample" not in meta.columns:
        raise ValueError("sample meta 必须包含 sample 列")
    merged = merged.merge(meta, on="sample", how="left")


# =========================================================
# 4. 输出样本明细
# =========================================================
detail_cols = [
    "candidate_id", "pacbio_srr", "transcript_id", "gene_id", "repeat_name", "te_role",
    "chr", "intron_start", "intron_end", "strand", "junction_type", "is_first_junction",
    "sample", "annotated", "motif_name",
    "unique_count", "multi_count", "total_count", "max_overhang",
    "support_loose", "support_moderate", "support_strict"
]
detail_df = merged[detail_cols].copy()
detail_df.to_csv(DETAIL_OUT, sep="\t", index=False)


# =========================================================
# 5. 输出每个 candidate 的汇总
# =========================================================
def first_non_na(x):
    x = x.dropna()
    return x.iloc[0] if len(x) > 0 else pd.NA

summary_df = (
    merged.groupby("candidate_id", as_index=False)
    .agg(
        pacbio_srr=("pacbio_srr", first_non_na),
        transcript_id=("transcript_id", first_non_na),
        gene_id=("gene_id", first_non_na),
        repeat_name=("repeat_name", first_non_na),
        te_role=("te_role", first_non_na),
        chr=("chr", first_non_na),
        intron_start=("intron_start", "first"),
        intron_end=("intron_end", "first"),
        strand=("strand", first_non_na),
        junction_type=("junction_type", first_non_na),
        is_first_junction=("is_first_junction", "max"),

        n_samples_detected=("support_loose", "sum"),
        n_samples_moderate=("support_moderate", "sum"),
        n_samples_strict=("support_strict", "sum"),

        best_unique_count=("unique_count", "max"),
        best_total_count=("total_count", "max"),
        best_overhang=("max_overhang", "max"),
        any_annotated=("annotated", "max")
    )
)

best_sample_df = (
    merged.sort_values(
        ["candidate_id", "unique_count", "total_count", "max_overhang"],
        ascending=[True, False, False, False]
    )
    .drop_duplicates("candidate_id")
    [["candidate_id", "sample"]]
    .rename(columns={"sample": "best_sample"})
)

summary_df = summary_df.merge(best_sample_df, on="candidate_id", how="left")
summary_df.to_csv(SUMMARY_OUT, sep="\t", index=False)

print("Done.")
print("Detail :", DETAIL_OUT)
print("Summary:", SUMMARY_OUT)