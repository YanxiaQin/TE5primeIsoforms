# 功能: 生成包含所有讨论过的列的最终报告，包含外显子跳跃的合并等

import pandas as pd
import re
import os
import glob

def describe_boundary_relationship(row):
    if row['is_multi_exon']: return "跨越多个外显子"
    tolerance = 5
    te_start, te_end = row['te_genomic_start'], row['te_genomic_end']
    exon_start, exon_end = row['host_exon_start'], row['host_exon_end']
    is_start_boundary = abs(te_start - exon_start) <= tolerance
    is_end_boundary = abs(te_end - exon_end) <= tolerance
    if is_start_boundary and is_end_boundary: return "几乎完全等同于外显子"
    elif is_start_boundary: return "构成了外显子的5'端"
    elif is_end_boundary: return "构成了外显子的3'端"
    elif te_start > exon_start and te_end < exon_end: return "完全包含于外显子内部"
    else: return "部分重叠"

def get_genomic_context(row, gene_db_df):
    te_start, te_end, chrom = row['te_genomic_start'], row['te_genomic_end'], row['chrom']
    overlapping_genes = gene_db_df[(gene_db_df['chrom'] == chrom) & (gene_db_df['start'] <= te_end) & (gene_db_df['end'] >= te_start)]
    if overlapping_genes.empty:
        up_genes = gene_db_df[(gene_db_df['chrom'] == chrom) & (gene_db_df['end'] < te_start)].sort_values(by='end', ascending=False)
        down_genes = gene_db_df[(gene_db_df['chrom'] == chrom) & (gene_db_df['start'] > te_end)].sort_values(by='start', ascending=True)
        up_gene, down_gene = (up_genes['gene_id'].iloc[0] if not up_genes.empty else 'N/A'), (down_genes['gene_id'].iloc[0] if not down_genes.empty else 'N/A')
        return f"推断的基因间区 (位于基因 {up_gene} 和 {down_gene} 之间)"
    return f"位于基因 {row['gene_id']} 区域内"

def describe_te_truncation(row):
    """
    根据TE在参考序列上的比对坐标，准确判断其截断状态。
    该函数直接利用 RepeatMasker .out 文件中 begin 和 left 列的含义。
    """
    if pd.isna(row['reference_length']):
        return "无参考长度"
    
    begin_str = str(row['repeat_begin'])
    left_str = str(row['repeat_left'])
    
    try:
        begin = int(begin_str.strip('()'))
        left = int(left_str.strip('()'))
    except (ValueError, TypeError):
        return "坐标格式错误"

    tolerance = 5
    
    is_5prime_intact = (begin <= tolerance)
    
    is_3prime_intact = (left <= tolerance)

    if is_5prime_intact and is_3prime_intact:
        return "完整的"
    elif is_5prime_intact:
        return "3'端截断"
    elif is_3prime_intact:
        return "5'端截断"
    else:
        return "内部片段"

# ==============================================================================
# --- 批量处理多个样本 ---
# ==============================================================================

def process_single_sample_report(input_mapped_csv, ref_len_csv, pacbio_gtf_file, output_final_csv):
    """
    封装了原始main函数的核心逻辑，用于生成单个样本的最终报告。
    """
    print("--- 步骤1: 加载文件 ---")
    try:
        fragments_df = pd.read_csv(input_mapped_csv)
        ref_len_df = pd.read_csv(ref_len_csv)
        gtf_df = pd.read_csv(pacbio_gtf_file, sep='\t', comment='#', header=None, names=['chrom','s','feature','start','end','sc','strand','fr','attr'])
    except FileNotFoundError as e: 
        print(f"错误: 找不到必要的输入文件 {e.filename}。跳过此样本。")
        return

    print("--- 步骤2: 构建基因和转录本数据库 ---")
    gtf_df['gene_id'] = gtf_df['attr'].apply(lambda x: re.search(r'gene_id "([^"]+)"', x).group(1))
    gene_db = {gene_id: {'chrom': g['chrom'].iloc[0], 'start': g['start'].min(), 'end': g['end'].max()} for gene_id, g in gtf_df[gtf_df['feature'] == 'transcript'].groupby('gene_id')}
    gene_db_df = pd.DataFrame.from_dict(gene_db, orient='index').reset_index().rename(columns={'index':'gene_id'})
    
    exon_df = gtf_df[gtf_df['feature'] == 'exon'].copy()
    exon_df['transcript_id'] = exon_df['attr'].apply(lambda x: re.search(r'transcript_id "([^"]+)"', x).group(1))
    transcript_info = {tx_id: {'gene_id': g['gene_id'].iloc[0], 'exon_count': len(g), 'length': (g['end'] - g['start'] + 1).sum()} for tx_id, g in exon_df.groupby('transcript_id')}
    
    print("--- 步骤3: 聚合跨外显子片段并计算指标 ---")
    final_df = fragments_df.merge(ref_len_df, on='repeat_name', how='left')
    
    agg_funcs = { col: 'first' for col in final_df.columns if col not in ['fragment_id', 'te_genomic_start', 'te_genomic_end', 'host_exon_start', 'host_exon_end', 'host_exon_index'] }
    agg_funcs.update({'te_genomic_start': 'min', 'te_genomic_end': 'max', 'host_exon_start': 'min', 'host_exon_end': 'max', 'host_exon_index': lambda x: ','.join(map(str, sorted(x.unique())))})
    
    aggregated_df = final_df.groupby('fragment_id', as_index=False).agg(agg_funcs)

    aggregated_df['te_length'] = aggregated_df.apply(lambda r: (r['te_end_in_tx'] - r['te_start_in_tx'] + 1), axis=1)
    
    total_te_load_per_tx = aggregated_df.groupby('transcript_id')['te_length'].sum().to_dict()
    
    final_results = []
    for _, row in aggregated_df.iterrows():
        tx_id, row_dict = row['transcript_id'], row.to_dict()
        if tx_id not in transcript_info: continue # 安全检查，防止GTF和输入文件不匹配
        tx_info = transcript_info[tx_id]
        
        row_dict['gene_id'] = tx_info['gene_id']
        row_dict['tx_exon_count'] = tx_info['exon_count']
        
        row_dict['is_multi_exon'] = len(str(row['host_exon_index']).split(',')) > 1
        
        row_dict['genomic_context'] = get_genomic_context(row_dict, gene_db_df)
        
        idx, exon_count = str(row['host_exon_index']), tx_info['exon_count']
        if row_dict['is_multi_exon']:
            row_dict['location_in_transcript'] = f"跨越外显子 {idx}"
        else:
            exon_type = "内部" if 1 < int(idx) < exon_count else ("起始" if int(idx) == 1 else "末端")
            row_dict['location_in_transcript'] = f"位于第 {idx}/{exon_count} 个外显子 ({exon_type}外显子)"
        
        row_dict['boundary_relationship'] = describe_boundary_relationship(row_dict)
        
        host_exon_len = row['host_exon_end'] - row['host_exon_start'] + 1
        row_dict['te_coverage_on_host_exon'] = f"{row['te_length'] / host_exon_len:.2%}" if host_exon_len > 0 and not row_dict['is_multi_exon'] else '跨外显子N/A'
        
        tx_len = tx_info['length']
        row_dict['te_coverage_on_transcript'] = f"{row['te_length'] / tx_len:.2%}" if tx_len > 0 else 'N/A'
        
        total_te_load = total_te_load_per_tx.get(tx_id, 0)
        row_dict['total_te_coverage_on_transcript'] = f"{total_te_load / tx_len:.2%}" if tx_len > 0 else 'N/A'
        
        row_dict['te_truncation_status'] = describe_te_truncation(row)
        
        ref_len = row['reference_length']
        row_dict['te_completeness'] = f"{row['te_length'] / ref_len:.2%}" if pd.notna(ref_len) and ref_len > 0 else 'N/A'

        final_results.append(row_dict)
        
    final_report_df = pd.DataFrame(final_results)

    print(f"--- 步骤4: 保存最终报告 ---")
    final_report_df = final_report_df.sort_values(by=['transcript_id', 'te_start_in_tx']).reset_index(drop=True)
    
    final_cols = [
        'fragment_id', 'transcript_id', 'gene_id', 'tx_exon_count', 
        'repeat_name', 'repeat_class', 
        'genomic_context', 'location_in_transcript', 'boundary_relationship', 'te_truncation_status',
        'te_length', 'reference_length', 'te_completeness',
        'te_coverage_on_host_exon', 'te_coverage_on_transcript', 'total_te_coverage_on_transcript',
        'chrom', 'te_genomic_start', 'te_genomic_end'
    ]
    final_report_df = final_report_df.reindex(columns=final_cols)
    
    output_dir = os.path.dirname(output_final_csv)
    os.makedirs(output_dir, exist_ok=True)
    
    final_report_df.to_csv(output_final_csv, index=False)
    print(f"样本报告已生成: '{output_final_csv}'。")


def batch_process_reports():
    """
    新的主函数，用于发现所有样本并为它们生成最终报告。
    """
    # **************************************************************************
    base_dir = '/home/yxqin/neoantigen/Analysis/pacbio/05_repeatmasker/PRJNA851801'
    ref_te_dir = '/home/yxqin/neoantigen/Analysis/pacbio/05_repeatmasker'
    # **************************************************************************
    
    # 定义全局共享的TE参考长度文件
    ref_len_csv = os.path.join(ref_te_dir, 'te_reference_lengths.csv')
    if not os.path.exists(ref_len_csv):
        print(f"错误：全局TE参考长度文件未找到: {ref_len_csv}")
        return

    print(f"开始在根目录中搜索样本: {base_dir}")
    
    # 搜索所有SRR样本文件夹
    search_pattern = os.path.join(base_dir, '**', 'SRR*')
    sample_dirs = [d for d in glob.glob(search_pattern, recursive=True) if os.path.isdir(d)]
    
    if not sample_dirs:
        print("错误：在指定目录下未找到任何SRR样本文件夹。请检查 'base_dir' 路径是否正确。")
        return

    print(f"发现 {len(sample_dirs)} 个潜在的样本目录。开始处理...")
    
    # 遍历所有找到的样本目录
    for sample_path in sample_dirs:
        srr_id = os.path.basename(sample_path)
        print(f"\n==================================================")
        print(f"正在处理样本: {srr_id}")
        print(f"==================================================")
        
        # 构建每个样本所需文件的完整路径
        input_mapped_csv = os.path.join(sample_path, 'RM_parse', f"{srr_id}_parseRM_te_mapped_full.csv")
        pacbio_gtf_file = os.path.join(sample_path, f"{srr_id}_final.gtf")
        output_final_csv = os.path.join(sample_path, 'RM_parse', f"{srr_id}_parseRM_te_merge_info.csv")
        
        # 调用处理单个样本的函数
        process_single_sample_report(input_mapped_csv, ref_len_csv, pacbio_gtf_file, output_final_csv)

    print("\n所有样本处理完毕！")


if __name__ == '__main__':
    batch_process_reports()
