#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from pathlib import Path
import pandas as pd
import re

# =========================================================
# 0. 路径设置
# =========================================================
BASE_SQANTI = Path("/home/yxqin/neoantigen/Analysis/pacbio/04_sqanti3")

EVENT_XLSX = Path("/home/yxqin/neoantigen/Analysis/pacbio/09_common_epitopes/deduplication_isoform/get_transcript_and_orf_seq/get_ORF_origin/TE_event_organized_workbook_v2.xlsx")
EVENT_SHEET = "Raw_with_EventIDs"

MERGE_TE_XLSX = Path("/home/yxqin/neoantigen/Analysis/pacbio/09_common_epitopes/deduplication_isoform/get_transcript_and_orf_seq/get_ORF_origin/merge_Initial_mainRoleTE.xlsx")
MERGE_TE_SHEET = 0

OUT_DIR = Path("/home/yxqin/neoantigen/Analysis/pacbio/09_common_epitopes/deduplication_isoform/get_transcript_and_orf_seq/get_ORF_origin")
OUT_XLSX = OUT_DIR / "ORF_origin_exon_junction_TE.xlsx"
OUT_TSV = OUT_DIR / "ORF_origin_exon_junction_TE.tsv"


# =========================================================
# 1. 工具函数
# =========================================================
def collapse_unique(x, sep=","):
    vals = []
    for v in x:
        if pd.isna(v):
            continue
        s = str(v).strip()
        if s != "" and s.lower() != "nan":
            vals.append(s)
    vals = sorted(set(vals))
    return sep.join(vals) if vals else pd.NA


def parse_gtf_attr(attr_text):
    """
    GTF 第9列属性解析
    例如：
    transcript_id "PB.1.1"; gene_id "PB.1";
    """
    attrs = {}
    for m in re.finditer(r'(\S+)\s+"([^"]+)"', attr_text):
        attrs[m.group(1)] = m.group(2)
    return attrs


def overlap_len(a_start, a_end, b_start, b_end):
    s = max(int(a_start), int(b_start))
    e = min(int(a_end), int(b_end))
    return max(0, e - s + 1)


def find_srr_dir(srr_id):
    """
    自动在 PRJNA851801 / PRJNA723287 / PRJNA1176011 下找 SRR 目录
    """
    hits = list(BASE_SQANTI.glob(f"PRJNA*/{srr_id}"))
    if len(hits) == 0:
        return None
    if len(hits) > 1:
        print(f"[Warning] {srr_id} 找到多个目录，默认使用第一个: {hits[0]}")
    return hits[0]


def load_gtf_exons(gtf_file, target_tx_ids=None):
    """
    读取 GTF 里的 exon
    只保留目标 transcript
    """
    rows = []
    with open(gtf_file, "r") as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            arr = line.rstrip("\n").split("\t")
            if len(arr) != 9:
                continue
            chrom, source, feature, start, end, score, strand, frame, attrs = arr
            if feature != "exon":
                continue

            attr_dict = parse_gtf_attr(attrs)
            tx_id = attr_dict.get("transcript_id")
            gene_id = attr_dict.get("gene_id")

            if target_tx_ids is not None and tx_id not in target_tx_ids:
                continue

            rows.append({
                "chrom": chrom,
                "start": int(start),
                "end": int(end),
                "strand": strand,
                "transcript_id": tx_id,
                "gene_id_gtf": gene_id
            })
    return pd.DataFrame(rows)


def build_tx_exon_map(exon_df_one_tx):
    """
    给单条 transcript 的 exon 建立 transcript 坐标映射
    返回列表，每个元素对应一个 exon：
    {
      exon_number, chrom, strand,
      exon_genomic_start, exon_genomic_end,
      exon_length,
      tx_start, tx_end
    }
    """
    if exon_df_one_tx.empty:
        return []

    chroms = set(exon_df_one_tx["chrom"])
    strands = set(exon_df_one_tx["strand"])
    if len(chroms) != 1 or len(strands) != 1:
        raise ValueError("同一 transcript 的 exon 出现在多个 chrom 或多个 strand，不合理。")

    chrom = exon_df_one_tx["chrom"].iloc[0]
    strand = exon_df_one_tx["strand"].iloc[0]
    exons = exon_df_one_tx.sort_values(["start", "end"]).reset_index(drop=True).copy()

    if strand == "+":
        tx_order = exons.copy()
    elif strand == "-":
        tx_order = exons.iloc[::-1].reset_index(drop=True).copy()
    else:
        raise ValueError(f"非法 strand: {strand}")

    exon_map = []
    tx_cursor = 1
    for i, row in tx_order.iterrows():
        exon_len = int(row["end"]) - int(row["start"]) + 1
        tx_start = tx_cursor
        tx_end = tx_cursor + exon_len - 1

        exon_map.append({
            "exon_number": i + 1,
            "chrom": chrom,
            "strand": strand,
            "exon_genomic_start": int(row["start"]),
            "exon_genomic_end": int(row["end"]),
            "exon_length": exon_len,
            "tx_start": tx_start,
            "tx_end": tx_end
        })
        tx_cursor = tx_end + 1

    return exon_map


def map_cds_to_exons_and_genome(exon_map, cds_start, cds_end):
    """
    把 CDS 的 transcript 坐标映射到：
    1) 哪些 exon
    2) 每个 exon 上占了多少 nt
    3) genomic blocks

    返回：
    exon_blocks: 每个 overlap exon 的详细信息
    genomic_blocks: ORF 在 genome 上的 block 列表
    """
    exon_blocks = []
    genomic_blocks = []

    for ex in exon_map:
        ov_tx_start = max(int(cds_start), int(ex["tx_start"]))
        ov_tx_end = min(int(cds_end), int(ex["tx_end"]))

        if ov_tx_start > ov_tx_end:
            continue

        ov_len = ov_tx_end - ov_tx_start + 1

        if ex["strand"] == "+":
            g_start = ex["exon_genomic_start"] + (ov_tx_start - ex["tx_start"])
            g_end = ex["exon_genomic_start"] + (ov_tx_end - ex["tx_start"])
        else:
            g_end = ex["exon_genomic_end"] - (ov_tx_start - ex["tx_start"])
            g_start = ex["exon_genomic_end"] - (ov_tx_end - ex["tx_start"])

        exon_block = {
            "exon_number": ex["exon_number"],
            "chrom": ex["chrom"],
            "strand": ex["strand"],
            "exon_genomic_start": ex["exon_genomic_start"],
            "exon_genomic_end": ex["exon_genomic_end"],
            "exon_tx_start": ex["tx_start"],
            "exon_tx_end": ex["tx_end"],
            "cds_overlap_tx_start": ov_tx_start,
            "cds_overlap_tx_end": ov_tx_end,
            "cds_overlap_nt": ov_len,
            "cds_overlap_genomic_start": int(min(g_start, g_end)),
            "cds_overlap_genomic_end": int(max(g_start, g_end)),
        }
        exon_blocks.append(exon_block)

        genomic_blocks.append({
            "chrom": ex["chrom"],
            "strand": ex["strand"],
            "block_start": int(min(g_start, g_end)),
            "block_end": int(max(g_start, g_end)),
            "exon_number": ex["exon_number"],
            "block_nt": ov_len
        })

    return exon_blocks, genomic_blocks


def summarize_orf_origin(exon_blocks):
    """
    从 exon_blocks 里总结：
    - 来自哪个 exon
    - 是否跨 junction
    """
    if len(exon_blocks) == 0:
        return {
            "orf_start_exon": pd.NA,
            "orf_end_exon": pd.NA,
            "orf_exon_numbers": pd.NA,
            "orf_exon_nt_breakdown": pd.NA,
            "major_exon": pd.NA,
            "major_exon_nt": pd.NA,
            "crosses_junction": pd.NA,
            "junctions_crossed": pd.NA
        }

    exon_numbers = [x["exon_number"] for x in exon_blocks]
    breakdown = [f"exon{b['exon_number']}:{b['cds_overlap_nt']}nt" for b in exon_blocks]

    major = max(exon_blocks, key=lambda x: x["cds_overlap_nt"])

    return {
        "orf_start_exon": exon_blocks[0]["exon_number"],
        "orf_end_exon": exon_blocks[-1]["exon_number"],
        "orf_exon_numbers": ",".join([f"exon{n}" for n in exon_numbers]),
        "orf_exon_nt_breakdown": "; ".join(breakdown),
        "major_exon": f"exon{major['exon_number']}",
        "major_exon_nt": int(major["cds_overlap_nt"]),
        "crosses_junction": len(exon_blocks) > 1,
        "junctions_crossed": max(0, len(exon_blocks) - 1)
    }


def calc_te_overlap(genomic_blocks, te_df_one_tx):
    """
    计算 ORF genomic blocks 与 TE 的 overlap
    te_df_one_tx 至少包含：
    chrom, te_genomic_start, te_genomic_end, repeat_name, repeat_class
    """
    te_detail_rows = []
    total_overlap_nt = 0

    if te_df_one_tx.empty or len(genomic_blocks) == 0:
        return {
            "orf_overlap_te": False,
            "te_overlap_nt": 0,
            "te_overlap_fraction": 0.0,
            "overlap_repeat_names": pd.NA,
            "overlap_repeat_classes": pd.NA
        }, pd.DataFrame(te_detail_rows)

    for block in genomic_blocks:
        for _, te in te_df_one_tx.iterrows():
            if str(block["chrom"]) != str(te["chrom"]):
                continue

            ov = overlap_len(
                block["block_start"], block["block_end"],
                te["te_genomic_start"], te["te_genomic_end"]
            )
            if ov > 0:
                total_overlap_nt += ov
                te_detail_rows.append({
                    "chrom": block["chrom"],
                    "strand": block["strand"],
                    "exon_number": block["exon_number"],
                    "orf_block_start": block["block_start"],
                    "orf_block_end": block["block_end"],
                    "orf_block_nt": block["block_nt"],
                    "repeat_name": te["repeat_name"],
                    "repeat_class": te["repeat_class"],
                    "te_genomic_start": te["te_genomic_start"],
                    "te_genomic_end": te["te_genomic_end"],
                    "overlap_nt": ov
                })

    te_detail_df = pd.DataFrame(te_detail_rows)

    if te_detail_df.empty:
        return {
            "orf_overlap_te": False,
            "te_overlap_nt": 0,
            "te_overlap_fraction": 0.0,
            "overlap_repeat_names": pd.NA,
            "overlap_repeat_classes": pd.NA
        }, te_detail_df

    overlap_repeat_names = collapse_unique(te_detail_df["repeat_name"])
    overlap_repeat_classes = collapse_unique(te_detail_df["repeat_class"])

    return {
        "orf_overlap_te": True,
        "te_overlap_nt": int(total_overlap_nt),
        "te_overlap_fraction": float(total_overlap_nt),
        "overlap_repeat_names": overlap_repeat_names,
        "overlap_repeat_classes": overlap_repeat_classes
    }, te_detail_df


# =========================================================
# 2. 读取 Raw_with_EventIDs
#    这里至少需要 SRR_ID 和 transcript_id
# =========================================================
event_df = pd.read_excel(EVENT_XLSX, sheet_name=EVENT_SHEET, dtype=str)
event_df.columns = [str(c).strip() for c in event_df.columns]

required_event_cols = {"SRR_ID", "transcript_id"}
missing_event = required_event_cols - set(event_df.columns)
if missing_event:
    raise ValueError(f"{EVENT_SHEET} 缺少必要列: {missing_event}")

keep_event_cols = ["SRR_ID", "transcript_id"]
if "event_id" in event_df.columns:
    keep_event_cols.append("event_id")

event_df = event_df[keep_event_cols].copy()
event_df["SRR_ID"] = event_df["SRR_ID"].astype(str).str.strip()
event_df["transcript_id"] = event_df["transcript_id"].astype(str).str.strip()

# 同一个 SRR_ID + transcript_id 若对应多个 event_id，合并成逗号串
if "event_id" in event_df.columns:
    event_df = (
        event_df.groupby(["SRR_ID", "transcript_id"], as_index=False)
        .agg({"event_id": collapse_unique})
    )
else:
    event_df = event_df.drop_duplicates().reset_index(drop=True)


# =========================================================
# 3. 读取 TE 明细表
#    用于判断 ORF 是否 overlap TE
# =========================================================
te_df = pd.read_excel(MERGE_TE_XLSX, sheet_name=MERGE_TE_SHEET, dtype=str)
te_df.columns = [str(c).strip() for c in te_df.columns]

required_te_cols = {
    "SRR_ID", "transcript_id", "repeat_name", "repeat_class",
    "chrom", "te_genomic_start", "te_genomic_end"
}
missing_te = required_te_cols - set(te_df.columns)
if missing_te:
    raise ValueError(f"merge_Initial_mainRoleTE.xlsx 缺少必要列: {missing_te}")

te_df = te_df[list(required_te_cols)].copy()
te_df["SRR_ID"] = te_df["SRR_ID"].astype(str).str.strip()
te_df["transcript_id"] = te_df["transcript_id"].astype(str).str.strip()
te_df["repeat_name"] = te_df["repeat_name"].astype(str).str.strip()
te_df["repeat_class"] = te_df["repeat_class"].astype(str).str.strip()
te_df["chrom"] = te_df["chrom"].astype(str).str.strip()
te_df["te_genomic_start"] = pd.to_numeric(te_df["te_genomic_start"], errors="coerce")
te_df["te_genomic_end"] = pd.to_numeric(te_df["te_genomic_end"], errors="coerce")
te_df = te_df.dropna(subset=["te_genomic_start", "te_genomic_end"]).copy()


# =========================================================
# 4. 主循环：按 SRR 处理
# =========================================================
summary_rows = []
exon_detail_rows = []
te_overlap_detail_rows = []
missing_rows = []

for srr_id, sub_targets in event_df.groupby("SRR_ID"):
    srr_dir = find_srr_dir(srr_id)

    if srr_dir is None:
        for _, rr in sub_targets.iterrows():
            missing_rows.append({
                "SRR_ID": srr_id,
                "transcript_id": rr["transcript_id"],
                "event_id": rr["event_id"] if "event_id" in rr.index else pd.NA,
                "reason": "SRR directory not found under 04_sqanti3/PRJNA*/"
            })
        continue

    gtf_file = srr_dir / f"{srr_id}.filtered.gtf"
    cls_file = srr_dir / f"{srr_id}_RulesFilter_result_classification.txt"

    if not gtf_file.exists():
        for _, rr in sub_targets.iterrows():
            missing_rows.append({
                "SRR_ID": srr_id,
                "transcript_id": rr["transcript_id"],
                "event_id": rr["event_id"] if "event_id" in rr.index else pd.NA,
                "reason": f"GTF not found: {gtf_file}"
            })
        continue

    if not cls_file.exists():
        for _, rr in sub_targets.iterrows():
            missing_rows.append({
                "SRR_ID": srr_id,
                "transcript_id": rr["transcript_id"],
                "event_id": rr["event_id"] if "event_id" in rr.index else pd.NA,
                "reason": f"Classification file not found: {cls_file}"
            })
        continue

    # 目标 transcript 集合
    target_txs = set(sub_targets["transcript_id"])

    # ---------- 读 classification ----------
    cls_df = pd.read_csv(cls_file, sep="\t", dtype=str)
    cls_df.columns = [str(c).strip() for c in cls_df.columns]

    required_cls_cols = {
        "isoform", "coding", "ORF_length", "CDS_length",
        "CDS_start", "CDS_end", "CDS_genomic_start", "CDS_genomic_end", "ORF_seq"
    }
    missing_cls_cols = required_cls_cols - set(cls_df.columns)
    if missing_cls_cols:
        raise ValueError(f"{cls_file} 缺少必要列: {missing_cls_cols}")

    cls_df = cls_df[cls_df["isoform"].isin(target_txs)].copy()
    cls_df["coding"] = cls_df["coding"].astype(str).str.strip().str.lower()

    cls_df = cls_df[cls_df["coding"] == "coding"].copy()

    for col in ["ORF_length", "CDS_length", "CDS_start", "CDS_end", "CDS_genomic_start", "CDS_genomic_end"]:
        cls_df[col] = pd.to_numeric(cls_df[col], errors="coerce")

    dup_isoforms = cls_df["isoform"][cls_df["isoform"].duplicated()].unique().tolist()
    if dup_isoforms:
        print(f"[Warning] {srr_id} classification 中这些 isoform 出现多次，默认保留第一条: {dup_isoforms}")
        cls_df = cls_df.drop_duplicates(subset=["isoform"], keep="first").copy()

    cls_map = cls_df.set_index("isoform").to_dict(orient="index")

    # ---------- 读 GTF exon ----------
    gtf_exon_df = load_gtf_exons(gtf_file, target_tx_ids=target_txs)
    gtf_group = {
        tx: sub.copy()
        for tx, sub in gtf_exon_df.groupby("transcript_id")
    }

    # ---------- 逐 transcript 分析 ----------
    for _, row in sub_targets.iterrows():
        tx_id = row["transcript_id"]
        event_id = row["event_id"] if "event_id" in row.index else pd.NA

        if tx_id not in cls_map:
            missing_rows.append({
                "SRR_ID": srr_id,
                "transcript_id": tx_id,
                "event_id": event_id,
                "reason": "No retained coding ORF in _RulesFilter_result_classification.txt"
            })
            continue

        if tx_id not in gtf_group:
            missing_rows.append({
                "SRR_ID": srr_id,
                "transcript_id": tx_id,
                "event_id": event_id,
                "reason": "Transcript not found in filtered.gtf"
            })
            continue

        cls_row = cls_map[tx_id]
        cds_start = int(cls_row["CDS_start"])
        cds_end = int(cls_row["CDS_end"])
        cds_length = int(cls_row["CDS_length"])
        orf_length_aa = int(cls_row["ORF_length"])
        cds_g_start = int(cls_row["CDS_genomic_start"])
        cds_g_end = int(cls_row["CDS_genomic_end"])
        orf_seq_aa = cls_row["ORF_seq"]

        exon_map = build_tx_exon_map(gtf_group[tx_id])

        tx_len_from_gtf = sum([x["exon_length"] for x in exon_map])

        if cds_start < 1 or cds_end > tx_len_from_gtf or cds_start > cds_end:
            missing_rows.append({
                "SRR_ID": srr_id,
                "transcript_id": tx_id,
                "event_id": event_id,
                "reason": f"Invalid CDS coordinates vs GTF transcript length. CDS={cds_start}-{cds_end}, tx_len={tx_len_from_gtf}"
            })
            continue

        exon_blocks, genomic_blocks = map_cds_to_exons_and_genome(exon_map, cds_start, cds_end)
        orf_origin = summarize_orf_origin(exon_blocks)

        # exon blocks 明细
        for b in exon_blocks:
            exon_detail_rows.append({
                "SRR_ID": srr_id,
                "event_id": event_id,
                "transcript_id": tx_id,
                "chrom": b["chrom"],
                "strand": b["strand"],
                "exon_number": b["exon_number"],
                "exon_genomic_start": b["exon_genomic_start"],
                "exon_genomic_end": b["exon_genomic_end"],
                "exon_tx_start": b["exon_tx_start"],
                "exon_tx_end": b["exon_tx_end"],
                "cds_overlap_tx_start": b["cds_overlap_tx_start"],
                "cds_overlap_tx_end": b["cds_overlap_tx_end"],
                "cds_overlap_nt": b["cds_overlap_nt"],
                "cds_overlap_genomic_start": b["cds_overlap_genomic_start"],
                "cds_overlap_genomic_end": b["cds_overlap_genomic_end"]
            })

        te_sub = te_df[(te_df["SRR_ID"] == srr_id) & (te_df["transcript_id"] == tx_id)].copy()
        te_summary, te_detail_df = calc_te_overlap(genomic_blocks, te_sub)

        if cds_length > 0:
            te_summary["te_overlap_fraction"] = te_summary["te_overlap_nt"] / cds_length
        else:
            te_summary["te_overlap_fraction"] = 0.0

        if not te_detail_df.empty:
            te_detail_df.insert(0, "SRR_ID", srr_id)
            te_detail_df.insert(1, "event_id", event_id)
            te_detail_df.insert(2, "transcript_id", tx_id)
            te_overlap_detail_rows.extend(te_detail_df.to_dict(orient="records"))

        genomic_block_str = "; ".join([
            f"{g['chrom']}:{g['block_start']}-{g['block_end']}(exon{g['exon_number']},{g['block_nt']}nt)"
            for g in genomic_blocks
        ])

        summary_rows.append({
            "SRR_ID": srr_id,
            "event_id": event_id,
            "transcript_id": tx_id,

            "chrom": exon_map[0]["chrom"],
            "strand": exon_map[0]["strand"],
            "tx_exon_count_gtf": len(exon_map),
            "tx_length_from_gtf": tx_len_from_gtf,

            "ORF_length_aa": orf_length_aa,
            "CDS_length_nt": cds_length,
            "CDS_start_tx": cds_start,
            "CDS_end_tx": cds_end,
            "CDS_genomic_start_from_classification": cds_g_start,
            "CDS_genomic_end_from_classification": cds_g_end,
            "ORF_seq_aa": orf_seq_aa,

            "orf_start_exon": orf_origin["orf_start_exon"],
            "orf_end_exon": orf_origin["orf_end_exon"],
            "orf_exon_numbers": orf_origin["orf_exon_numbers"],
            "orf_exon_nt_breakdown": orf_origin["orf_exon_nt_breakdown"],
            "major_exon": orf_origin["major_exon"],
            "major_exon_nt": orf_origin["major_exon_nt"],
            "crosses_junction": orf_origin["crosses_junction"],
            "junctions_crossed": orf_origin["junctions_crossed"],

            "orf_genomic_blocks": genomic_block_str,

            "orf_overlap_te": te_summary["orf_overlap_te"],
            "te_overlap_nt": te_summary["te_overlap_nt"],
            "te_overlap_fraction": te_summary["te_overlap_fraction"],
            "overlap_repeat_names": te_summary["overlap_repeat_names"],
            "overlap_repeat_classes": te_summary["overlap_repeat_classes"]
        })


# =========================================================
# 5. 输出
# =========================================================
summary_df = pd.DataFrame(summary_rows)
exon_detail_df = pd.DataFrame(exon_detail_rows)
te_overlap_detail_df = pd.DataFrame(te_overlap_detail_rows)
missing_df = pd.DataFrame(missing_rows)

if summary_df.empty:
    raise ValueError("没有成功生成任何 ORF 来源结果，请检查输入路径和文件内容。")

summary_df = summary_df.sort_values(["SRR_ID", "transcript_id"]).reset_index(drop=True)

# TSV
summary_df.to_csv(OUT_TSV, sep="\t", index=False)

# Excel
with pd.ExcelWriter(OUT_XLSX, engine="openpyxl") as writer:
    summary_df.to_excel(writer, index=False, sheet_name="ORF_origin_summary")
    exon_detail_df.to_excel(writer, index=False, sheet_name="ORF_exon_detail")
    te_overlap_detail_df.to_excel(writer, index=False, sheet_name="ORF_TE_overlap_detail")
    missing_df.to_excel(writer, index=False, sheet_name="Missing_or_skipped")

print("Done.")
print(f"Summary TSV : {OUT_TSV}")
print(f"Summary XLSX: {OUT_XLSX}")