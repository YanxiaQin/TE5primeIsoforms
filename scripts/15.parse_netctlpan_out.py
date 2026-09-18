# 文件名: run_netctlpan_parser.py
# 功能: 自动解析指定PRJ项目下所有SRR样本的netCTLpan输出，
#       提取CTL表位，并整合总结信息，生成综合性的CSV报告。

import pandas as pd
import sys
import re
import os
import glob

def parse_netctlpan_output(input_file):
    if not os.path.exists(input_file) or os.path.getsize(input_file) == 0:
        print(f"  -> 警告: 输入文件为空或不存在，跳过: {os.path.basename(input_file)}")
        return pd.DataFrame(), pd.DataFrame()

    epitopes, summaries = [], {}
    
    # 正则表达式
    epitope_line_pattern = re.compile(
        r"^\s*(\d+)\s+"                 # 1: Position
        r"([^\s]+)\s+"                 # 2: Protein_ID
        r"(HLA-.+?)\s+"                # 3: HLA_Allele
        r"([A-Z]+)\s+"                 # 4: Peptide
        r"([\d\.-]+)\s+"               # 5: MHC_Score
        r"([\d\.-]+)\s+"               # 6: TAP_Score
        r"([\d\.-]+)\s+"               # 7: Cleavage_Score
        r"([\d\.-]+)\s+"               # 8: Combined_Score
        r"([\d\.]+)\s+<-E"             # 9: MHC_Rank
    )
    summary_line_pattern = re.compile(
        r"Number of MHC ligands\s+(\d+)\s+identified\.\s+"
        r"Number of peptides\s+(\d+)\.\s+"
        r"Allele\s+(HLA-.+?)\.\s+"
        r"Protein name\s+([^\s]+)"
    )

    with open(input_file, 'r') as f:
        for line in f:
            epitope_match = epitope_line_pattern.search(line)
            if epitope_match:
                groups = epitope_match.groups()
                epitopes.append({
                    'Position': int(groups[0]), 'Protein_ID': groups[1], 'HLA_Allele': groups[2],
                    'Peptide': groups[3], 'MHC_Score': float(groups[4]), 'TAP_Score': float(groups[5]),
                    'Cleavage_Score': float(groups[6]), 'Combined_Score': float(groups[7]),
                    'MHC_Rank': float(groups[8])
                })
                continue

            summary_match = summary_line_pattern.search(line)
            if summary_match:
                groups = summary_match.groups()
                key = (groups[3], groups[2])
                summaries[key] = {
                    'Num_CTL_Epitopes': int(groups[0]),
                    'Total_Peptides_Tested': int(groups[1])
                }
                
    epitopes_df = pd.DataFrame(epitopes)
    
    if summaries:
        summary_list = [{'Protein_ID': k[0], 'HLA_Allele': k[1], **v} for k, v in summaries.items()]
        summary_df = pd.DataFrame(summary_list)
    else:
        summary_df = pd.DataFrame()
        
    return epitopes_df, summary_df


def process_srr_for_parsing(project_id, srr_id, config):
    print(f"\n--- 开始解析样本: {srr_id} ---")
    
    for data_type in ["known", "novel"]:
        print(f"\n  处理类型: {data_type}")
        
        input_file = os.path.join(config['input_base_dir'], project_id, srr_id, f"{srr_id}_netCTLpan.{data_type}.out")
        output_dir = os.path.join(config['output_base_dir'], project_id, srr_id)
        
        os.makedirs(output_dir, exist_ok=True)
        
        epitopes_df, summary_df = parse_netctlpan_output(input_file)
        
        if epitopes_df.empty:
            print("  -> 未找到任何CTL表位，不生成报告文件。")
            continue
            
        print(f"  -> 共找到 {len(epitopes_df)} 个CTL表位和 {len(summary_df)} 条总结记录。")
        
        if not summary_df.empty:
            comprehensive_df = pd.merge(
                epitopes_df,
                summary_df,
                on=['Protein_ID', 'HLA_Allele'],
                how='left'
            )
        else:
            comprehensive_df = epitopes_df
            comprehensive_df['Num_CTL_Epitopes'] = 'N/A'
            comprehensive_df['Total_Peptides_Tested'] = 'N/A'

        comprehensive_df = comprehensive_df.sort_values(by='Combined_Score', ascending=False)
        
        output_csv_path = os.path.join(output_dir, f"{srr_id}_netCTLpan_{data_type}_epitopes_report.csv")
        comprehensive_df.to_csv(output_csv_path, index=False)
        print(f"  -> 综合报告已保存到: {os.path.basename(output_csv_path)}")

def main():
    
    PRJ_ID = "PRJNA1176011"
    
    SRR_IDS_TO_PROCESS = []
    
    INPUT_BASE_DIR = "/home/yxqin/neoantigen/Analysis/pacbio/08_netCTLpan"
    OUTPUT_BASE_DIR = "/home/yxqin/neoantigen/Analysis/pacbio/08_netCTLpan/parse"

    config = {
        "input_base_dir": INPUT_BASE_DIR,
        "output_base_dir": OUTPUT_BASE_DIR,
    }
    
    print(f"====== 开始解析项目: {PRJ_ID} ======")
    
    prj_input_dir = os.path.join(config['input_base_dir'], PRJ_ID)
    if not os.path.isdir(prj_input_dir):
        print(f"错误: 找不到指定的项目输入目录: {prj_input_dir}")
        sys.exit(1)
        
    if not SRR_IDS_TO_PROCESS:
        print("未指定SRR列表，将自动发现在项目目录下的所有SRR样本...")
        search_pattern = os.path.join(prj_input_dir, "SRR*")
        srr_dirs = sorted([d for d in glob.glob(search_pattern) if os.path.isdir(d)])
        srr_ids = [os.path.basename(d) for d in srr_dirs]
    else:
        srr_ids = SRR_IDS_TO_PROCESS
        print(f"将处理用户指定的 {len(srr_ids)} 个SRR样本...")
    
    if not srr_ids:
        print(f"错误: 在 {prj_input_dir} 下未找到任何SRR样本目录或用户未指定。")
        sys.exit(1)
        
    print(f"目标样本列表: {', '.join(srr_ids)}")
    
    for srr_id in srr_ids:
        try:
            process_srr_for_parsing(PRJ_ID, srr_id, config)
        except Exception as e:
            print(f"\n!!!!!! 处理样本 {srr_id} 时发生严重错误: {e} !!!!!!")
            print("!!!!!! 继续处理下一个样本... !!!!!!")
            continue
            
    print(f"\n====== 项目 {PRJ_ID} 解析完成 ======")

if __name__ == '__main__':
    main()