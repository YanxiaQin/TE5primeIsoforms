#!/bin/bash

# ==============================================================================
#                 RepeatMasker 批量处理脚本
#
# 功能:
# 1. 自动查找所有 SRR* 样本目录中经过筛选的 ..._final_filtered.fasta 文件。
# 2. 在每个样本目录内，创建 RM_out 子目录并运行 RepeatMasker。
# 3. 控制并发任务数量，防止服务器过载。
# 4. 为每个任务生成独立的日志文件。
# 5. 设计为使用 nohup 在后台安全运行。
# ==============================================================================

# ---配置---

# 1. ID (例如: PRJNA1176011)
PROJECT_ID="PRJNA1176011"

# 2. 包含所有SRR*样本文件夹的【输入】项目目录
BASE_INPUT_DIR="/home/yxqin/neoantigen/Analysis/pacbio/05_repeatmasker"

# 3. 日志文件存放目录 (脚本会自动在其中为该项目创建子目录)
LOG_ROOT_DIR="/home/yxqin/neoantigen/Analysis/pacbio/00_log/repeatmasker"

# 4. 最大并发任务数
MAX_JOBS=2


# --- 自动构建完整路径 ---
INPUT_PROJECT_DIR="${BASE_INPUT_DIR}/${PROJECT_ID}"
LOG_DIR="${LOG_ROOT_DIR}/${PROJECT_ID}"

# 检查并创建日志目录
mkdir -p "$LOG_DIR"

echo "================================================="
echo "RepeatMasker 批量处理启动"
echo "处理项目: ${PROJECT_ID}"
echo "输入/工作目录: ${INPUT_PROJECT_DIR}"
echo "日志目录: ${LOG_DIR}"
echo "最大并发数: ${MAX_JOBS}"
echo "================================================="
echo ""

# 检查输入项目目录是否存在
if [ ! -d "$INPUT_PROJECT_DIR" ]; then
    echo "错误: 输入项目目录不存在: '$INPUT_PROJECT_DIR'"
    exit 1
fi

for SAMPLE_DIR in "$INPUT_PROJECT_DIR"/SRR*; do
    if [ ! -d "$SAMPLE_DIR" ]; then
        echo "--> 警告: 在 '$INPUT_PROJECT_DIR' 中没有找到任何名为 'SRR*' 的样本目录。"
        continue
    fi

    # 检查当前后台运行的任务数量
    while [[ $(jobs -p | wc -l) -ge $MAX_JOBS ]]; do
        wait -n
    done

    # 从完整路径中提取样本ID (例如 SRR31072027)
    SAMPLE_ID=$(basename "${SAMPLE_DIR%/}")
    
    # 构建该样本的输入 FASTA 文件路径
    INPUT_FASTA="${SAMPLE_DIR}/${SAMPLE_ID}_final_filtered.fasta"

    if [ ! -f "$INPUT_FASTA" ]; then
        echo "--> 警告: 在目录 '$SAMPLE_DIR' 中找不到输入文件 '$INPUT_FASTA'，已跳过此样本。"
        continue
    fi
    
    # 构建日志文件路径
    LOG_FILE="${LOG_DIR}/${SAMPLE_ID}_repeatmasker.log"

    echo "-------------------------------------------------"
    echo "准备启动样本: $SAMPLE_ID"
    echo "输入文件: $INPUT_FASTA"
    echo "日志文件: $LOG_FILE"
    
    (
        echo "开始处理样本: $SAMPLE_ID at $(date)"
        
        # 定义输出目录，并确保它存在
        # RepeatMasker 的 -dir 选项会自动创建目录，但自己创建更保险
        OUTPUT_RM_DIR="${SAMPLE_DIR}/RM_out"
        mkdir -p "$OUTPUT_RM_DIR"

        echo "--> 运行 RepeatMasker..."
        
        # 执行 RepeatMasker 命令
        cd "$SAMPLE_DIR" || exit 1

        RepeatMasker "${SAMPLE_ID}_final_filtered.fasta" \
            -species "Homo sapiens" \
            -e rmblast \
            -s \
            -xsmall \
            -gff \
            -a \
            -pa 6
            
        mv ${SAMPLE_ID}_final_filtered.fasta.* RM_out/
            
        echo "样本 $SAMPLE_ID 处理完成 at $(date)"

    ) > "$LOG_FILE" 2>&1 &

    echo "样本 $SAMPLE_ID 已在后台启动。PID: $!"
    echo ""
    sleep 2

done

echo "所有样本任务均已提交。"
echo "正在等待最后一批任务完成..."
wait

echo ""
echo "================================================="
echo "所有 RepeatMasker 任务已全部完成！"
echo "请检查各个样本目录下的 'RM_out' 子目录获取结果。"
echo "运行日志保存在 '${LOG_DIR}' 目录中。"
echo "================================================="