#!/bin/bash

# ===================================================================================================================
#                 SQANTI3 filter的结果(gtf,faa,fasta,classification.txt等文件进行初步筛选：常染色体+XYM,Isoform,orf lenth >=8
# ===================================================================================================================

# --- 配置---
PROJECT_ID="PRJNA851801"
BASE_INPUT_DIR="/home/yxqin/neoantigen/Analysis/pacbio/04_sqanti3"
BASE_OUTPUT_DIR="/home/yxqin/neoantigen/Analysis/pacbio/05_repeatmasker"
LOG_ROOT_DIR="/home/yxqin/neoantigen/Analysis/pacbio/00_log/sqanti3/sqanti3_filter/PRJNA851801"
MAX_JOBS=2


INPUT_PROJECT_DIR="${BASE_INPUT_DIR}/${PROJECT_ID}"
OUTPUT_PROJECT_DIR="${BASE_OUTPUT_DIR}/${PROJECT_ID}"
LOG_DIR="${LOG_ROOT_DIR}/${PROJECT_ID}"

mkdir -p "$OUTPUT_PROJECT_DIR"
mkdir -p "$LOG_DIR"

echo "================================================="
echo "处理项目: ${PROJECT_ID}"
echo "================================================="
echo ""

if [ ! -d "$INPUT_PROJECT_DIR" ]; then
    echo "错误: 输入项目目录不存在: '$INPUT_PROJECT_DIR'"
    exit 1
fi

for SAMPLE_DIR in "$INPUT_PROJECT_DIR"/SRR*/; do
    if [ ! -d "$SAMPLE_DIR" ]; then
        echo "--> 警告: 在 '$INPUT_PROJECT_DIR' 中没有找到任何名为 'SRR*' 的样本目录。"
        continue
    fi

    while [[ $(jobs -p | wc -l) -ge $MAX_JOBS ]]; do
        wait -n
    done

    SAMPLE_ID=$(basename "${SAMPLE_DIR%/}")
    OUTPUT_SAMPLE_DIR="${OUTPUT_PROJECT_DIR}/${SAMPLE_ID}"
    mkdir -p "$OUTPUT_SAMPLE_DIR"
    LOG_FILE="${LOG_DIR}/${SAMPLE_ID}_filter_run.log"

    echo "-------------------------------------------------"
    echo "准备启动样本: $SAMPLE_ID"
    echo "输出将保存至: $OUTPUT_SAMPLE_DIR"
    echo "日志文件: $LOG_FILE"
    
    (
        echo "开始处理样本: $SAMPLE_ID at $(date)"
        
        INPUT_CLASSIFICATION=$(ls -1 "${SAMPLE_DIR}/SRR"*"_RulesFilter_result_classification.txt" 2>/dev/null | head -n 1)
        INPUT_GTF=$(ls -1 "${SAMPLE_DIR}/SRR"*"filtered.gtf" 2>/dev/null | head -n 1)
        INPUT_FASTA=$(ls -1 "${SAMPLE_DIR}/SRR"*"filtered.fasta" 2>/dev/null | head -n 1)
        INPUT_FAA=$(ls -1 "${SAMPLE_DIR}/SRR"*"filtered.faa" 2>/dev/null | head -n 1)

        for f in "$INPUT_CLASSIFICATION" "$INPUT_GTF" "$INPUT_FASTA" "$INPUT_FAA"; do
            if [ -z "$f" ] || [ ! -f "$f" ]; then
                echo "错误: 在目录 $SAMPLE_DIR 中找不到必需的输入文件。脚本已为该样本终止。"
                exit 1
            fi
        done

        OUTPUT_ID_LIST="${OUTPUT_SAMPLE_DIR}/${SAMPLE_ID}_final_isoform_ids.txt"
        OUTPUT_CLASSIFICATION="${OUTPUT_SAMPLE_DIR}/${SAMPLE_ID}_final_filtered_classification.txt"
        OUTPUT_GTF="${OUTPUT_SAMPLE_DIR}/${SAMPLE_ID}_final.gtf"
        OUTPUT_NUC_FASTA="${OUTPUT_SAMPLE_DIR}/${SAMPLE_ID}_final_filtered.fasta"
        OUTPUT_FAA="${OUTPUT_SAMPLE_DIR}/${SAMPLE_ID}_final_filtered.faa"
        
        echo "--> 步骤1: 使用相对列位置进行筛选..."
        
        head -n 1 "$INPUT_CLASSIFICATION" > "$OUTPUT_CLASSIFICATION"
        > "$OUTPUT_ID_LIST"

        # ===========================================================================
        awk -F'\t' '
            {
                sub(/\r$/,"", $NF);
                sub(/\r$/,"", $(NF-2));
            }
            NR > 1 && $NF == "Isoform" &&
            $2 ~ /^chr([1-9]|1[0-9]|2[0-2]|X|Y|M)$/ &&
            $(NF-2) != "NA" && length($(NF-2)) >= 8 {
                print $0 >> "'"$OUTPUT_CLASSIFICATION"'"
                print $1 >> "'"$OUTPUT_ID_LIST"'"
            }
        ' "$INPUT_CLASSIFICATION"
        # ===========================================================================

        if [ ! -s "$OUTPUT_ID_LIST" ]; then
            echo "--> 警告: 样本 $SAMPLE_ID 没有找到任何符合条件的 isoform。后续步骤已跳过。"
            rm "$OUTPUT_CLASSIFICATION" "$OUTPUT_ID_LIST" 2>/dev/null
            exit 0
        fi
        echo "    找到了 $(wc -l < "$OUTPUT_ID_LIST") 个 isoform。"

        echo "--> 步骤2: 筛选 GTF, FASTA, 和 FAA 文件..."
        grep -Fwf "$OUTPUT_ID_LIST" "$INPUT_GTF" > "$OUTPUT_GTF"
        /home/yxqin/software/seqkit/seqkit grep -f "$OUTPUT_ID_LIST" "$INPUT_FASTA" -o "$OUTPUT_NUC_FASTA"
        /home/yxqin/software/seqkit/seqkit grep -f "$OUTPUT_ID_LIST" "$INPUT_FAA" -o "$OUTPUT_FAA"
        
        echo "样本 $SAMPLE_ID 处理完成 at $(date)"

    ) > "$LOG_FILE" 2>&1 &

    echo "样本 $SAMPLE_ID 已在后台启动。PID: $!"
    echo ""
    sleep 2
done

echo "所有样本任务均已提交。正在等待最后一批任务完成..."
wait
echo ""
echo "================================================="
echo "所有样本处理任务已全部完成！"
echo "================================================="