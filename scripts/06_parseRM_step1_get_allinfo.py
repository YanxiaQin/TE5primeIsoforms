import pandas as pd
import re
import os
import glob

def parse_pacbio_gtf(gtf_file):
    print(f"--- 正在解析PacBio GTF: {gtf_file} ---")
    try:
        gtf_df = pd.read_csv(gtf_file, sep='\t', comment='#', header=None, names=['chrom','s','feature','start','end','sc','strand','fr','attr'])
    except FileNotFoundError: print(f"错误: GTF文件 '{gtf_file}' 未找到。"); exit()
    transcript_db = {}
    exon_df = gtf_df[gtf_df['feature'] == 'exon'].copy()
    exon_df['transcript_id'] = exon_df['attr'].apply(lambda x: re.search(r'transcript_id "([^"]+)"', x).group(1))
    for tx_id, group in exon_df.groupby('transcript_id'):
        sorted_exons = group.sort_values(by='start').reset_index(drop=True)
        transcript_db[tx_id] = {'strand': group['strand'].iloc[0], 'exons': sorted_exons}
    print("GTF解析完成。")
    return transcript_db

def parse_repeatmasker_out_advanced(out_file):
    print(f"--- 正在解析RepeatMasker .out文件: {out_file} ---")
    data = []
    try:
        with open(out_file, 'r') as f:
            for i, line in enumerate(f, 1):
                line = line.strip()
                if not line or not line[0].isdigit(): continue
                parts = re.split(r'\s+', line)
                if len(parts) < 11: continue

                align_strand = '+' if parts[8] != 'C' else 'C'
                
                try:
                    if align_strand == 'C':
                        repeat_left = parts[11]; repeat_end = parts[12]; repeat_begin = parts[13]
                    else:
                        repeat_begin = parts[11]; repeat_end = parts[12]; repeat_left = parts[13]
                except IndexError:
                    repeat_begin, repeat_end, repeat_left = 'N/A', 'N/A', 'N/A'

                data.append({
                    'fragment_id': f"frag_line_{i}", 'transcript_id': parts[4],
                    'te_start_in_tx': int(parts[5]), 'te_end_in_tx': int(parts[6]),
                    'align_strand': align_strand, 'repeat_name': parts[9], 'repeat_class': parts[10],
                    'repeat_begin': repeat_begin, 'repeat_end': repeat_end, 'repeat_left': repeat_left
                })
    except FileNotFoundError: print(f"错误: .out文件 '{out_file}' 未找到。"); exit()
    print(f"解析完成 {len(data)} 个TE片段。")
    return pd.DataFrame(data)


# ==============================================================================
# --- 批量处理多个样本 ---
# ==============================================================================

def process_single_sample(pacbio_gtf_file, rm_out_file, output_intermediate_csv):
    """
    封装了原始main函数核心逻辑的函数，用于处理单个样本。
    接收输入和输出文件路径作为参数。
    """
    # 检查输入文件是否存在，如果不存在则跳过此样本
    if not os.path.exists(pacbio_gtf_file):
        print(f"警告: GTF文件未找到，跳过处理: {pacbio_gtf_file}")
        return
    if not os.path.exists(rm_out_file):
        print(f"警告: .out文件未找到，跳过处理: {rm_out_file}")
        return

    transcript_db = parse_pacbio_gtf(pacbio_gtf_file)
    te_df = parse_repeatmasker_out_advanced(rm_out_file)
    
    print("--- 步骤1: 核心坐标映射 ---")
    all_mapped_parts = []
    for _, te_row in te_df.iterrows():
        tx_id = te_row['transcript_id']
        if tx_id not in transcript_db: continue
        tx_info = transcript_db[tx_id]
        exons, tx_strand = tx_info['exons'], tx_info['strand']
        
        remaining_te_start, remaining_te_end = te_row['te_start_in_tx'], te_row['te_end_in_tx']
        tx_pos_counter = 0
        
        exon_iterator = exons.iterrows() if tx_strand == '+' else exons.iloc[::-1].iterrows()
        biological_exon_counter = 0

        for exon_idx, exon in exon_iterator:
            biological_exon_counter += 1
            exon_len = exon['end'] - exon['start'] + 1
            exon_tx_start, exon_tx_end = tx_pos_counter + 1, tx_pos_counter + exon_len
            overlap_start, overlap_end = max(exon_tx_start, remaining_te_start), min(exon_tx_end, remaining_te_end)

            if overlap_start <= overlap_end:
                part_offset_start = overlap_start - tx_pos_counter
                part_offset_end = overlap_end - tx_pos_counter
                g_start, g_end = (exon['start'] + part_offset_start - 1, exon['start'] + part_offset_end - 1) if tx_strand == '+' else (exon['end'] - part_offset_end + 1, exon['end'] - part_offset_start + 1)
                
                all_mapped_parts.append({
                    'fragment_id': te_row['fragment_id'], 'transcript_id': tx_id,
                    'te_start_in_tx': te_row['te_start_in_tx'], 'te_end_in_tx': te_row['te_end_in_tx'],
                    'repeat_name': te_row['repeat_name'], 'repeat_class': te_row['repeat_class'],
                    'align_strand': te_row['align_strand'], 'tx_strand': tx_strand,
                    'repeat_begin': te_row['repeat_begin'], 'repeat_end': te_row['repeat_end'], 'repeat_left': te_row['repeat_left'],
                    'te_genomic_start': min(g_start, g_end), 'te_genomic_end': max(g_start, g_end),
                    'chrom': exon['chrom'], 'host_exon_start': exon['start'], 'host_exon_end': exon['end'],
                    'host_exon_index': biological_exon_counter
                })
            tx_pos_counter += exon_len
            if remaining_te_end < exon_tx_end: break
                
    output_df = pd.DataFrame(all_mapped_parts)

    output_dir = os.path.dirname(output_intermediate_csv)
    os.makedirs(output_dir, exist_ok=True)

    print(f"--- 映射完成，保存中间结果到 {output_intermediate_csv} ---")
    output_df.to_csv(output_intermediate_csv, index=False)
    print(f"样本 {os.path.basename(pacbio_gtf_file).split('_final.gtf')[0]} 处理完成！")

def batch_process_main():
    """
    新的主函数，用于发现所有样本并对它们进行批量处理。
    """
    # **************************************************************************
    # *** 根目录 ***
    base_dir = '/home/yxqin/neoantigen/Analysis/pacbio/05_repeatmasker/PRJNA851801'
    # **************************************************************************

    print(f"开始在根目录中搜索样本: {base_dir}")
    
    search_pattern = os.path.join(base_dir, '**', 'SRR*')
    
    sample_dirs = [d for d in glob.glob(search_pattern, recursive=True) if os.path.isdir(d)]
    
    if not sample_dirs:
        print("错误：在指定目录下未找到任何SRR样本文件夹。请检查 'base_dir' 路径是否正确。")
        return

    print(f"发现 {len(sample_dirs)} 个潜在的样本目录。开始处理...")
    
    for sample_path in sample_dirs:
        srr_id = os.path.basename(sample_path)
        print(f"\n==================================================")
        print(f"正在处理样本: {srr_id}")
        print(f"样本路径: {sample_path}")
        print(f"==================================================")
        
        pacbio_gtf_file = os.path.join(sample_path, f"{srr_id}_final.gtf")
        rm_out_file = os.path.join(sample_path, 'RM_out', f"{srr_id}_final_filtered.fasta.out")
        
        output_dir = os.path.join(sample_path, 'RM_parse')
        output_intermediate_csv = os.path.join(output_dir, f"{srr_id}_parseRM_te_mapped_full.csv")
        process_single_sample(pacbio_gtf_file, rm_out_file, output_intermediate_csv)

    print("\n所有样本处理完毕！")


if __name__ == '__main__':
    batch_process_main()