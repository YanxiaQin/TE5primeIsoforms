import pandas as pd
import os
import glob

def final_filtering_and_splitting(project_id, srr_id):
    """
    读取合并后的报告文件，根据一系列复杂的标准进行筛选，
    并将结果分流到 'novel' 和 'known' 两个文件中。
    此函数现在接收 project_id 和 srr_id 作为参数。
    """
    
    # --- 路径配置  ---
    # 输入文件
    input_file = f"/home/yxqin/neoantigen/Analysis/pacbio/05_repeatmasker/{project_id}/{srr_id}/RM_parse/{srr_id}_cancer_specific_merge_info_with_classification.csv"
    
    # 输出文件
    output_base_dir = f"/home/yxqin/neoantigen/Analysis/pacbio/06_netMHCpan/data/{project_id}/{srr_id}"
    novel_output_file = os.path.join(output_base_dir, f"{srr_id}_use_netMHCpan.novel.info.csv")
    known_output_file = os.path.join(output_base_dir, f"{srr_id}_use_netMHCpan.known.info.csv")
    # ---

    print(f"--- 开始处理文件: {os.path.basename(input_file)} ---")

    # 1. 读取数据
    try:
        df = pd.read_csv(input_file)
        print(f"成功读取 {len(df)} 行数据。")
    except FileNotFoundError:
        print(f"错误: 输入文件未找到，跳过此样本: {input_file}")
        return
    except Exception as e:
        print(f"读取文件时发生错误: {e}")
        return

    # 2. 应用所有通用筛选条件
    print("\n--- 步骤 1: 应用通用筛选条件 ---")
    
    cond_location = (
        df['location_in_transcript'].str.startswith('跨越外显子 1,', na=False) | 
        df['location_in_transcript'].str.startswith('位于第 1/', na=False) |
        df['location_in_transcript'].str.contains('起始外显子', na=False)
    )

    cond_truncation = df['te_truncation_status'] != '无参考长度'
    cond_ref_length = df['ref_length'].notna()
    cond_fl = df['FL'] >= 3
    cond_iso_exp = df['iso_exp'] >= 1
    cond_min_cov = df['min_cov'] >= 5
    cond_cage_peak = df['dist_to_CAGE_peak'].notna()
    cond_coding = df['coding'] == 'coding'

    common_mask = (
        cond_location & cond_truncation & cond_ref_length & cond_fl & 
        cond_iso_exp & cond_min_cov & cond_cage_peak & cond_coding
    )
    
    filtered_df = df[common_mask].copy()
    print(f"通用筛选后，剩余 {len(filtered_df)} 行数据。")
    
    if filtered_df.empty:
        print("没有数据满足通用筛选条件，将创建空的输出文件。")

    # 3. 根据 'structural_category' 进行分流
    print("\n--- 步骤 2: 根据 structural_category 进行分流 ---")
    
    novel_mask = ~filtered_df['structural_category'].isin(['full-splice_match', 'incomplete-splice_match'])
    novel_df = filtered_df[novel_mask]
    
    known_mask = filtered_df['structural_category'].isin(['full-splice_match', 'incomplete-splice_match'])
    known_df = filtered_df[known_mask]

    # 4. 保存结果
    print("\n--- 步骤 3: 保存筛选结果 ---")
    
    print(f"确保输出目录存在: {output_base_dir}")
    os.makedirs(output_base_dir, exist_ok=True)
    
    novel_df.to_csv(novel_output_file, index=False)
    print(f"已将 {len(novel_df)} 行 'Novel' 数据写入到: {novel_output_file}")
    
    known_df.to_csv(known_output_file, index=False)
    print(f"已将 {len(known_df)} 行 'Known' 数据写入到: {known_output_file}")


def batch_process_by_project():
    """
    主函数，用于处理指定项目下的所有SRR样本。
    """
    # ======================================================================
    # ===> 在这里指定想处理的PRJ项目ID <===
    target_project_id = "PRJNA851801"
    # ======================================================================

    input_project_dir = f"/home/yxqin/neoantigen/Analysis/pacbio/05_repeatmasker/{target_project_id}"
    
    print(f"开始在项目目录中搜索样本: {input_project_dir}")

    if not os.path.isdir(input_project_dir):
        print(f"错误: 项目目录不存在: {input_project_dir}")
        return
        
    search_pattern = os.path.join(input_project_dir, 'SRR*')
    sample_dirs = [d for d in glob.glob(search_pattern) if os.path.isdir(d)]

    if not sample_dirs:
        print(f"警告: 在目录 {input_project_dir} 下未找到任何SRR样本文件夹。")
        return

    print(f"发现 {len(sample_dirs)} 个样本目录。开始处理...")
    
    for sample_path in sample_dirs:
        srr_id = os.path.basename(sample_path)
        print(f"\n==================================================")
        print(f"正在处理样本: {srr_id} (项目: {target_project_id})")
        print(f"==================================================")
        
        final_filtering_and_splitting(target_project_id, srr_id)

    print("\n指定项目的所有样本处理完毕！")


if __name__ == '__main__':
    batch_process_by_project()
