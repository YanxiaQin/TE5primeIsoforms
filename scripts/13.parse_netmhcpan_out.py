# 功能: 针对指定的PRJ项目和SRR列表，自动为 "known" 和 "novel" 
#       两种类型的数据运行netMHCpan解析流程，生成定制化命名的6种报告文件。

import pandas as pd
import sys
import re
import os
import glob

def parse_netmhcpan_output(nohup_file):
    if not os.path.exists(nohup_file):
        print(f"  -> 警告: 输入文件未找到，跳过: {nohup_file}")
        return pd.DataFrame(), pd.DataFrame()
    all_binders, summary_results = [], []
    data_line_pattern = re.compile(
        r"^\s*(\d+)\s+(HLA-.+?)\s+([A-Z]+)\s+([A-Z]+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+"
        r"([A-Z\.]+)\s+([^\s]+)\s+([\d\.]+)\s+([\d\.]+)\s+([\d\.]+)\s+([\d\.]+)\s+([\d\.]+)\s+(<= (WB|SB))"
    )
    original_summary_pattern = re.compile(
        r"Protein\s+([^\s]+)\.\s+Allele\s+(HLA-.+?)\.\s+"
        r"Number of high binders\s+(\d+)\.\s+Number of weak binders\s+(\d+)\.\s+"
        r"Number of peptides\s+(\d+)"
    )
    with open(nohup_file, 'r') as f:
        for line in f:
            data_match = data_line_pattern.search(line)
            if data_match:
                groups = data_match.groups()
                all_binders.append({'Pos': int(groups[0]), 'MHC': groups[1].strip(), 'Peptide': groups[2], 'Core': groups[3],'Of': int(groups[4]), 'Gp': int(groups[5]), 'Gl': int(groups[6]), 'Ip': int(groups[7]),'Il': int(groups[8]), 'Icore': groups[9], 'Identity': groups[10],'Score_EL': float(groups[11]), 'Rank_EL': float(groups[12]),'Score_BA': float(groups[13]), 'Rank_BA': float(groups[14]),'Aff(nM)': float(groups[15]), 'BindLevel': groups[17]})
                continue
            summary_match = original_summary_pattern.search(line)
            if summary_match:
                groups = summary_match.groups()
                summary_results.append({'Protein_ID': groups[0], 'MHC_Allele': groups[1].strip(),'Num_Strong_Binders': int(groups[2]),'Num_Weak_Binders': int(groups[3]),'Total_Peptides_Tested': int(groups[4])})
    return pd.DataFrame(all_binders), pd.DataFrame(summary_results)

def generate_reports(detailed_df, original_summary_df, base_path, srr_id, data_type, suffix, is_filtered=False):
    if detailed_df.empty:
        return
    detailed_path = os.path.join(base_path, f"{srr_id}_{data_type}_netmhcpan_{suffix}_detailed.csv")
    summary_path = os.path.join(base_path, f"{srr_id}_{data_type}_netmhcpan_{suffix}_summary.csv")
    comprehensive_path = os.path.join(base_path, f"{srr_id}_{data_type}_netmhcpan_{suffix}_comprehensive.csv")
    detailed_df_sorted = detailed_df.sort_values(by='Rank_EL', ascending=True)
    detailed_df_sorted.to_csv(detailed_path, index=False)
    print(f"  -> 已保存: {os.path.basename(detailed_path)} ({len(detailed_df_sorted)} 条)")
    
    if is_filtered:
        agg_dict = {
            'BindLevel': [
                ('Num_Strong_Binders', lambda x: (x == 'SB').sum()),
                ('Num_Weak_Binders', lambda x: (x == 'WB').sum())
            ]
        }
        summary_df = detailed_df.groupby(['Identity', 'MHC'], as_index=False).agg(agg_dict)
        summary_df.columns = ['Identity', 'MHC', 'Num_Strong_Binders', 'Num_Weak_Binders']
        summary_df = summary_df.rename(columns={'Identity': 'Protein_ID', 'MHC': 'MHC_Allele'})
    else:
        summary_df = original_summary_df.copy()
    if not summary_df.empty:
        summary_df_sorted = summary_df.sort_values(by=['Num_Strong_Binders', 'Num_Weak_Binders'], ascending=False)
        summary_df_sorted.to_csv(summary_path, index=False)
        print(f"  -> 已保存: {os.path.basename(summary_path)} ({len(summary_df_sorted)} 条)")
    
    summary_cols_suffix = {'Num_Strong_Binders': f'Num_Strong_Binders_{suffix}','Num_Weak_Binders': f'Num_Weak_Binders_{suffix}'}
    summary_df_with_suffix = summary_df.rename(columns=summary_cols_suffix)
    comprehensive_df = pd.merge(detailed_df_sorted, summary_df_with_suffix,left_on=['Identity', 'MHC'], right_on=['Protein_ID', 'MHC_Allele'], how='left').drop(columns=['Protein_ID', 'MHC_Allele'], errors='ignore')
    if not original_summary_df.empty:
        comprehensive_df = pd.merge(comprehensive_df, original_summary_df[['Protein_ID', 'MHC_Allele', 'Total_Peptides_Tested']],left_on=['Identity', 'MHC'], right_on=['Protein_ID', 'MHC_Allele'], how='left').drop(columns=['Protein_ID', 'MHC_Allele'], errors='ignore')
    comprehensive_df.to_csv(comprehensive_path, index=False)
    print(f"  -> 已保存: {os.path.basename(comprehensive_path)} ({len(comprehensive_df)} 条)")

def process_single_srr(prj_id, srr_id, config):
    print(f"\n--- 开始处理样本: {srr_id} ---")
    for data_type in ["known", "novel"]:
        print(f"\n  处理类型: {data_type}")
        input_file = os.path.join(config['input_base_dir'], prj_id, srr_id, f"{srr_id}_netMHCpan.{data_type}.out")
        output_dir = os.path.join(config['output_base_dir'], prj_id, srr_id)
        os.makedirs(output_dir, exist_ok=True)
        all_binders_df, original_summary_df = parse_netmhcpan_output(input_file)
        if all_binders_df.empty:
            print(f"  在 {os.path.basename(input_file)} 中未找到结合肽，跳过此类型。")
            continue
        generate_reports(all_binders_df, original_summary_df, output_dir, srr_id, data_type, "all_out", is_filtered=False)
        df_filtered = all_binders_df.copy()
        df_filtered = df_filtered[df_filtered['Rank_EL'] <= config['filtering']['el_rank_threshold']]
        df_filtered = df_filtered[df_filtered['Rank_BA'] <= config['filtering']['ba_rank_threshold']]
        df_filtered = df_filtered[df_filtered['Aff(nM)'] <= config['filtering']['affinity_threshold_nm']]
        if not df_filtered.empty:
            generate_reports(df_filtered, original_summary_df, output_dir, srr_id, data_type, "filtered_out", is_filtered=True)
        else:
            print("\n  筛选后没有剩下任何肽段，不生成 'filtered_out' 系列报告。")

def main():
    PRJ_ID = "PRJNA1176011"
    SRR_IDS_TO_PROCESS = []
    INPUT_BASE_DIR = "/home/yxqin/neoantigen/Analysis/pacbio/06_netMHCpan/output"
    OUTPUT_BASE_DIR = "/home/yxqin/neoantigen/Analysis/pacbio/06_netMHCpan/parse"
    FILTERING_CONFIG = {"el_rank_threshold": 2.0,"ba_rank_threshold": 2.0,"affinity_threshold_nm": 500.0}
    config = {"input_base_dir": INPUT_BASE_DIR,"output_base_dir": OUTPUT_BASE_DIR,"filtering": FILTERING_CONFIG}
    print(f"====== 开始处理项目: {PRJ_ID} ======")
    prj_input_dir = os.path.join(config['input_base_dir'], PRJ_ID)
    if not os.path.isdir(prj_input_dir):
        print(f"错误: 找不到指定的项目输入目录: {prj_input_dir}")
        sys.exit(1)
    if not SRR_IDS_TO_PROCESS:
        print("未指定SRR列表，将自动发现在项目目录下的所有SRR样本...")
        srr_dirs = glob.glob(os.path.join(prj_input_dir, "SRR*"))
        srr_ids = sorted([os.path.basename(d) for d in srr_dirs])
    else:
        srr_ids = SRR_IDS_TO_PROCESS
        print(f"将处理用户指定的 {len(srr_ids)} 个SRR样本...")
    if not srr_ids:
        print(f"错误: 在 {prj_input_dir} 下未找到任何SRR样本目录或用户未指定。")
        sys.exit(1)
    print(f"目标样本列表: {', '.join(srr_ids)}")
    for srr_id in srr_ids:
        try:
            process_single_srr(PRJ_ID, srr_id, config)
        except Exception as e:
            print(f"\n!!!!!! 处理样本 {srr_id} 时发生严重错误: {e} !!!!!!")
            print("!!!!!! 继续处理下一个样本... !!!!!!")
            continue
    print(f"\n====== 项目 {PRJ_ID} 处理完成 ======")

if __name__ == '__main__':
    main()