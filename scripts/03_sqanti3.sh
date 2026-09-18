#!/bin/bash

# --- 配置参数 ---
# 包含 Iso-Seq collapse 结果 GTF 文件的输入目录
INPUT_DIR="/home/yxqin/neoantigen/Analysis/pacbio/03_collapse/PRJNA1176011"
# SQANTI3 分析结果的根输出目录
OUTPUT_BASE_DIR="/home/yxqin/neoantigen/Analysis/pacbio/04_sqanti3/PRJNA1176011"
# 用于存放 short_reads.fofn 文件的目录
SHORT_READS_FOFN_DIR="/home/yxqin/neoantigen/Analysis/pacbio/04_sqanti3/PRJNA1176011/01_short_reads_fofn"
# Fastp 处理后的 Illumina 二代测序数据目录
ILLUMINA_FASTP_DIR="/home/yxqin/neoantigen/Analysis/pacbio/00_data/PRJNA1176011/Illumina/fastp"
# STAR 索引的源路径 
STAR_INDEX_SOURCE="/home/yxqin/neoantigen/Analysis/pacbio/04_sqanti3/PRJNA851801/SRR19785215/STAR_index"
# 参考基因组 GTF 文件
REFERENCE_GTF="/home/yxqin/neoantigen/reference/gencode.v47.primary_assembly.annotation.gtf"
# 参考基因组 FASTA 文件
REFERENCE_GENOME="/home/yxqin/neoantigen/reference/GRCh38.primary_assembly.genome.fa"
# CAGE Peak 文件
CAGE_PEAK="/home/yxqin/neoantigen/reference/refTSS/refTSS_v4.1_human_coordinate.hg38.bed"
# PolyA Motif 列表文件
POLYA_MOTIF="/home/yxqin/neoantigen/reference/polyA_motifs/mouse_and_human.polyA_motif.txt"
# 运行日志输出目录
LOG_DIR="/home/yxqin/neoantigen/Analysis/pacbio/00_log/sqanti3"
# 同时运行的任务数
MAX_JOBS=2

# 配对样本列表（二代-三代）
# 格式为 "Illumina_Sample_ID PacBio_Sample_ID"
PAIRED_SAMPLES=(
    "SRR31072018 SRR31072153" "SRR31072017 SRR31072151" "SRR31072016 SRR31072149"
    "SRR31072015 SRR31072147" "SRR31072099 SRR31072144" "SRR31072097 SRR31072142"
    "SRR31072096 SRR31072141" "SRR31072095 SRR31072139" "SRR31072094 SRR31072137"
    "SRR31072093 SRR31072135" "SRR31072092 SRR31072133" "SRR31072091 SRR31072130"
    "SRR31072090 SRR31072128" "SRR31072089 SRR31072127" "SRR31072088 SRR31072065"
    "SRR31072087 SRR31072063" "SRR31072086 SRR31072061" "SRR31072085 SRR31072056"
    "SRR31072084 SRR31072058" "SRR31072083 SRR31072054" "SRR31072014 SRR31072053"
    "SRR31072013 SRR31072051" "SRR31072010 SRR31072042" "SRR31072009 SRR31072044"
    "SRR31072008 SRR31072040" "SRR31072007 SRR31072039" "SRR31072004 SRR31072121"
    "SRR31072003 SRR31072119" "SRR31072002 SRR31072116" "SRR31072001 SRR31072114"
    "SRR31072000 SRR31072113" "SRR31071999 SRR31072111" "SRR31072082 SRR31072109"
    "SRR31072081 SRR31072107" "SRR31072080 SRR31072105" "SRR31072079 SRR31072102"
    "SRR31072078 SRR31072100" "SRR31072077 SRR31072098" "SRR31072076 SRR31072037"
    "SRR31072075 SRR31072035" "SRR31072074 SRR31072033" "SRR31072073 SRR31072028"
    "SRR31072072 SRR31072030" "SRR31072071 SRR31072026" "SRR31072070 SRR31072025"
    "SRR31072069 SRR31072023" "SRR31072068 SRR31072021" "SRR31072067 SRR31072019"
    "SRR31072018 SRR31072154" "SRR31072017 SRR31072152" "SRR31072016 SRR31072150"
    "SRR31072015 SRR31072148" "SRR31072099 SRR31072145" "SRR31072097 SRR31072143"
    "SRR31072096 SRR31072146" "SRR31072095 SRR31072140" "SRR31072094 SRR31072138"
    "SRR31072093 SRR31072136" "SRR31072092 SRR31072134" "SRR31072091 SRR31072131"
    "SRR31072090 SRR31072129" "SRR31072089 SRR31072132" "SRR31072088 SRR31072066"
    "SRR31072087 SRR31072064" "SRR31072086 SRR31072062" "SRR31072085 SRR31072057"
    "SRR31072084 SRR31072059" "SRR31072083 SRR31072055" "SRR31072014 SRR31072060"
    "SRR31072013 SRR31072052" "SRR31072010 SRR31072043" "SRR31072009 SRR31072045"
    "SRR31072008 SRR31072041" "SRR31072007 SRR31072046" "SRR31072004 SRR31072122"
    "SRR31072003 SRR31072120" "SRR31072002 SRR31072117" "SRR31072001 SRR31072115"
    "SRR31072000 SRR31072118" "SRR31071999 SRR31072112" "SRR31072082 SRR31072110"
    "SRR31072081 SRR31072108" "SRR31072080 SRR31072106" "SRR31072079 SRR31072103"
    "SRR31072078 SRR31072101" "SRR31072077 SRR31072104" "SRR31072076 SRR31072038"
    "SRR31072075 SRR31072036" "SRR31072074 SRR31072034" "SRR31072073 SRR31072029"
    "SRR31072072 SRR31072031" "SRR31072071 SRR31072027" "SRR31072070 SRR31072032"
    "SRR31072069 SRR31072024" "SRR31072068 SRR31072022" "SRR31072067 SRR31072020"
)

# --- 脚本设置 ---
set -euo pipefail 

echo "$(date '+%Y-%m-%d %H:%M:%S') Starting SQANTI3 analysis tasks."

mkdir -p "$SHORT_READS_FOFN_DIR" || { echo "Error: Cannot create directory $SHORT_READS_FOFN_DIR."; exit 1; }
mkdir -p "$LOG_DIR" || { echo "Error: Cannot create directory $LOG_DIR."; exit 1; }

TEMP_FIFO=$(mktemp -u)
mkfifo "$TEMP_FIFO" || { echo "Error: Cannot create FIFO $TEMP_FIFO."; exit 1; }
exec 3<>"$TEMP_FIFO"
rm -f "$TEMP_FIFO"

for ((i=0; i<MAX_JOBS; i++)); do
    echo >&3
done

# --- 函数：生成 short_reads.fofn 文件 ---
# 参数：PacBio样本ID, Illumina样本ID
generate_short_reads_fofn() {
    local pacbio_sample=$1
    local illumina_sample=$2
    local fofn_file="$SHORT_READS_FOFN_DIR/${pacbio_sample}.short_reads.fofn"
    local log_file_path="$LOG_DIR/${pacbio_sample}_sqanti3.log" 

    # 检查是否已存在该fofn文件
    if [ -f "$fofn_file" ]; then
        echo "$(date '+%Y-%m-%d %H:%M:%S') INFO: Existing short_reads.fofn file: $fofn_file" >> "$log_file_path"
        return 0
    fi

    # 构建二代测序文件路径
    local r1_file="$ILLUMINA_FASTP_DIR/${illumina_sample}_fastp_R1.fastq.gz"
    local r2_file="$ILLUMINA_FASTP_DIR/${illumina_sample}_fastp_R2.fastq.gz"

    # 检查二代测序文件是否存在
    if [ ! -f "$r1_file" ] || [ ! -f "$r2_file" ]; then
        echo "$(date '+%Y-%m-%d %H:%M:%S') ERROR: Illumina short reads not found for $illumina_sample. R1: $r1_file, R2: $r2_file" >> "$log_file_path"
        return 1
    fi

    # 写入fofn文件
    echo "$r1_file $r2_file" > "$fofn_file"
    echo "$(date '+%Y-%m-%d %H:%M:%S') INFO: Generated short_reads.fofn file: $fofn_file" >> "$log_file_path"
    return 0
}

# --- 函数：运行 SQANTI3 分析 ---
# 参数：PacBio样本ID, Illumina样本ID
run_sqanti3() {
    local pacbio_sample=$1
    local illumina_sample=$2
    local log_file="$LOG_DIR/${pacbio_sample}_sqanti3.log"

    echo "$(date '+%Y-%m-%d %H:%M:%S') INFO: Starting SQANTI3 processing for PacBio sample $pacbio_sample (Illumina paired with $illumina_sample)." | tee -a "$log_file"

    # 创建样本特定的输出目录
    local output_dir="$OUTPUT_BASE_DIR/$pacbio_sample"
    mkdir -p "$output_dir" || { echo "$(date '+%Y-%m-%d %H:%M:%S') ERROR: Cannot create output directory $output_dir for $pacbio_sample." | tee -a "$log_file"; return 1; }

    # 复制STAR_index到样本目录
    echo "$(date '+%Y-%m-%d %H:%M:%S') INFO: Copying STAR_index from $STAR_INDEX_SOURCE to $output_dir..." | tee -a "$log_file"
    cp -r "$STAR_INDEX_SOURCE" "$output_dir/" 2>> "$log_file"
    if [ $? -ne 0 ]; then
        echo "$(date '+%Y-%m-%d %H:%M:%S') ERROR: Failed to copy STAR_index for $pacbio_sample. Skipping." | tee -a "$log_file"
        return 1
    fi
    echo "$(date '+%Y-%m-%d %H:%M:%S') INFO: STAR_index copied." | tee -a "$log_file"

    # 输入 GTF 和 Abundance 文件路径
    local input_gtf="$INPUT_DIR/${pacbio_sample}.collapsed.gff"
    local input_abundance="$INPUT_DIR/${pacbio_sample}.collapsed.abundance.txt" 
    # 检查输入GTF文件是否存在
    if [ ! -f "$input_gtf" ]; then
        echo "$(date '+%Y-%m-%d %H:%M:%S') ERROR: Input GTF file '$input_gtf' not found for $pacbio_sample. Skipping." | tee -a "$log_file"
        return 1
    fi

    # 生成 short_reads.fofn 文件
    if ! generate_short_reads_fofn "$pacbio_sample" "$illumina_sample"; then
        echo "$(date '+%Y-%m-%d %H:%M:%S') ERROR: Failed to generate short_reads.fofn for $pacbio_sample. Skipping." | tee -a "$log_file"
        return 1
    fi
    local short_reads_fofn="$SHORT_READS_FOFN_DIR/${pacbio_sample}.short_reads.fofn"

    echo "$(date '+%Y-%m-%d %H:%M:%S') INFO: Running SQANTI3 analysis with Illumina short reads for $pacbio_sample..." | tee -a "$log_file"
    
    python /home/yxqin/software/SQANTI3_v5.3.6/sqanti3/sqanti3_qc.py \
        "$input_gtf" \
        "$REFERENCE_GTF" \
        "$REFERENCE_GENOME" \
        --CAGE_peak "$CAGE_PEAK" \
        --polyA_motif_list "$POLYA_MOTIF" \
        -o "$pacbio_sample" \
        -d "$output_dir" \
        -fl "$input_abundance" \
        --short_reads "$short_reads_fofn" \
        --isoAnnotLite \
        --cpus 6 \
        --report both \
        >> "$log_file" 2>&1 

    if [ $? -ne 0 ]; then
        echo "$(date '+%Y-%m-%d %H:%M:%S') ERROR: SQANTI3 analysis failed for $pacbio_sample!" | tee -a "$log_file"
        return 1
    fi
    echo "$(date '+%Y-%m-%d %H:%M:%S') INFO: SQANTI3 analysis completed successfully for $pacbio_sample." | tee -a "$log_file"
    return 0
}

# --- 主循环：处理所有配对样本 ---
for pair in "${PAIRED_SAMPLES[@]}"; do
    read -r illumina_sample pacbio_sample <<< "$pair"
    
    read -u3
    (
        run_sqanti3 "$pacbio_sample" "$illumina_sample"

        echo >&3
    ) &

    echo "$(date '+%Y-%m-%d %H:%M:%S') Submitted paired sample: PacBio $pacbio_sample (Illumina: $illumina_sample)"
done

echo "$(date '+%Y-%m-%d %H:%M:%S') Waiting for all background SQANTI3 tasks to finish..."
wait

exec 3>&-

echo "$(date '+%Y-%m-%d %H:%M:%S') All SQANTI3 analysis tasks finished!"
