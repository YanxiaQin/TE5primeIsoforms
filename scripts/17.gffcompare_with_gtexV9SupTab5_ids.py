import pandas as pd
import sys
import os
import glob
import shlex


def get_cancer_specific_ids(tracking_file):
    """
    从gffcompare的tracking文件中解析出所有结构新颖的转录本ID。
    (class_code既不等于'='也不等于'c')
    """
    print(f"正在解析 gffcompare tracking 文件: {tracking_file} ...")
    
    specific_ids = set()
    try:
        with open(tracking_file, 'r') as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) < 5:
                    continue
                class_code = parts[3]
                query_info = parts[4]
                
                if class_code not in ['=', 'c']:
                    try:
                        transcript_id = query_info.split('|')[1]
                        specific_ids.add(transcript_id)
                    except IndexError:
                        continue
                        
    except FileNotFoundError:
        print(f"错误: tracking 文件 '{tracking_file}' 未找到。")
        return None
        
    print(f"找到 {len(specific_ids)} 个在正常库中结构新颖的转录本ID。")
    return list(specific_ids)

def run_gffcompare(gffcompare_exe, ref_gtf, input_gtf, output_prefix):
    """
    构造并执行gffcompare命令。
    """
    print("--- 步骤 1: 运行 gffcompare ---")
    
    if not os.path.exists(gffcompare_exe):
        print(f"错误: gffcompare 执行文件未找到: {gffcompare_exe}")
        return False
    if not os.path.exists(ref_gtf):
        print(f"错误: 参考 GTF 文件未找到: {ref_gtf}")
        return False
    if not os.path.exists(input_gtf):
        print(f"错误: 输入 GTF 文件未找到: {input_gtf}")
        return False

    cmd = (
        f"{shlex.quote(gffcompare_exe)} -r {shlex.quote(ref_gtf)} "
        f"-o {shlex.quote(output_prefix)} {shlex.quote(input_gtf)}"
    )
    
    print(f"执行命令: {cmd}")
    exit_code = os.system(cmd)
    
    if exit_code == 0:
        print("gffcompare 成功完成。")
        return True
    else:
        print(f"错误: gffcompare 执行失败，退出码: {exit_code}。")
        return False

def filter_final_report(input_csv, specific_ids, output_csv):
    print(f"--- 步骤 3: 正在筛选报告文件: {os.path.basename(input_csv)} ---")
    try:
        report_df = pd.read_csv(input_csv)
        filtered_df = report_df[report_df['transcript_id'].isin(specific_ids)]
        
        filtered_df.to_csv(output_csv, index=False)
        
        if len(filtered_df) > 0:
            print(f"筛选完成，{len(filtered_df)} 行记录已保存到: {os.path.basename(output_csv)}")
        else:
            print(f"警告: 在 {os.path.basename(input_csv)} 中未找到匹配的特异性ID，输出文件为空。")

    except FileNotFoundError:
        print(f"错误: 报告文件未找到: {input_csv}")
    except KeyError:
        print(f"错误: 报告文件 {input_csv} 中缺少 'transcript_id' 列。")

def merge_with_classification_data(merge_info_csv, classification_txt, final_output_csv):
    """
    将 classification.txt 的信息精确合并到 merge_info.csv 中。
    """
    print(f"--- 步骤 4: 正在将 classification 数据合并到报告中 ---")
    
    try:
        merge_df = pd.read_csv(merge_info_csv)
        if merge_df.empty:
            print(f"警告: {os.path.basename(merge_info_csv)} 为空，跳过合并步骤。")
            merge_df.to_csv(final_output_csv, index=False)
            return
            
        class_df = pd.read_csv(classification_txt, sep='\t')
        
    except FileNotFoundError as e:
        print(f"错误: 读取文件失败，跳过合并步骤: {e.filename}")
        return

    if 'isoform' not in class_df.columns:
        print(f"错误: {os.path.basename(classification_txt)} 中缺少 'isoform' 列，无法合并。")
        return
    class_df.rename(columns={'isoform': 'transcript_id'}, inplace=True)

    final_df = pd.merge(merge_df, class_df, on='transcript_id', how='left')
    
    final_df.to_csv(final_output_csv, index=False)
    print(f"合并成功！包含详细分类信息的新报告已保存到: {os.path.basename(final_output_csv)}")


def process_single_sample(sample_path, srr_id, gffcompare_exe, ref_gtf):
    """
    为单个样本执行完整的 gffcompare, 报告筛选, 以及数据合并流程。
    """
    rm_parse_dir = os.path.join(sample_path, 'RM_parse')
    
    cancer_gtf = os.path.join(rm_parse_dir, f"{srr_id}_isoform_with_TE.gtf")
    gffcompare_prefix = os.path.join(rm_parse_dir, srr_id)
    
    if not run_gffcompare(gffcompare_exe, ref_gtf, cancer_gtf, gffcompare_prefix):
        print(f"由于 gffcompare 失败，已跳过样本 {srr_id} 的后续处理。")
        return

    print("--- 步骤 2: 提取癌症特异性转录本 ID ---")
    tracking_file = f"{gffcompare_prefix}.tracking"
    specific_ids = get_cancer_specific_ids(tracking_file)

    if not specific_ids:
        print(f"在样本 {srr_id} 中未找到任何特异性转录本，后续步骤已跳过。")
        return

    merge_info_filtered = os.path.join(rm_parse_dir, f"{srr_id}_parseRM_te_merge_info.filtered.csv")
    output_merge_specific = os.path.join(rm_parse_dir, f"{srr_id}_cancer_specific_merge_info.csv")
    mapped_full_filtered = os.path.join(rm_parse_dir, f"{srr_id}_parseRM_te_mapped_full.filtered.csv")
    output_mapped_specific = os.path.join(rm_parse_dir, f"{srr_id}_cancer_specific_mapped_full.csv")
    
    filter_final_report(merge_info_filtered, specific_ids, output_merge_specific)
    filter_final_report(mapped_full_filtered, specific_ids, output_mapped_specific)

    classification_file = os.path.join(sample_path, f"{srr_id}_final_filtered_classification.txt")
    final_merged_output = os.path.join(rm_parse_dir, f"{srr_id}_cancer_specific_merge_info_with_classification.csv")
    
    merge_with_classification_data(output_merge_specific, classification_file, final_merged_output)


def batch_processor_main():
    """
    主函数，用于发现所有样本并启动处理流程。
    """
    root_base_dir = '/home/yxqin/neoantigen/Analysis/pacbio/05_repeatmasker'
    base_dir = '/home/yxqin/neoantigen/Analysis/pacbio/05_repeatmasker/PRJNA1176011'
    gffcompare_executable = '/home/yxqin/miniconda3/envs/repeatmasker_neoantigen/bin/gffcompare'
    reference_gtf = os.path.join(root_base_dir, 'filtered_by_gtexV9SupTab5_ids.gtf')

    print(f"开始在项目目录中搜索样本: {base_dir}")
    search_pattern = os.path.join(base_dir, 'SRR*')
    sample_dirs = [d for d in glob.glob(search_pattern) if os.path.isdir(d)]

    if not sample_dirs:
        print(f"错误: 在目录 {base_dir} 下未找到任何SRR样本文件夹。")
        return

    print(f"发现 {len(sample_dirs)} 个潜在的样本目录。开始处理...")
    
    for sample_path in sample_dirs:
        srr_id = os.path.basename(sample_path)
        print(f"\n==================================================")
        print(f"正在处理样本: {srr_id}")
        print(f"==================================================")
        
        process_single_sample(sample_path, srr_id, gffcompare_executable, reference_gtf)

    print("\n所有样本处理完毕！")

if __name__ == '__main__':
    batch_processor_main()
