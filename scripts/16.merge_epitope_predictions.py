import pandas as pd
import os
import glob

def find_common_epitopes(project_id, srr_id):
    """
    为单个样本整合netCTLpan和netMHCpan的预测结果。
    """
    print(f"--- 正在处理样本: {srr_id} ---")
    
    # --- 路径定义 ---
    ctl_file = f"/home/yxqin/neoantigen/Analysis/pacbio/08_netCTLpan/parse/{project_id}/{srr_id}/{srr_id}_netCTLpan_novel_epitopes_report.csv"
    mhc_file = f"/home/yxqin/neoantigen/Analysis/pacbio/06_netMHCpan/parse/{project_id}/{srr_id}/{srr_id}_novel_netmhcpan_all_out_comprehensive.csv"
    output_dir = f"/home/yxqin/neoantigen/Analysis/pacbio/09_common_epitopes/{project_id}/{srr_id}"
    output_file = os.path.join(output_dir, f"{srr_id}_common_novel_epitopes_report.csv")
    # ---

    try:
        print(f"正在读取: {os.path.basename(ctl_file)}")
        ctl_df = pd.read_csv(ctl_file)
        
        print(f"正在读取: {os.path.basename(mhc_file)}")
        mhc_df = pd.read_csv(mhc_file)
    except FileNotFoundError as e:
        print(f"错误: 必需的输入文件未找到，跳过此样本: {e.filename}")
        return
        
    print(f"netCTLpan报告包含 {len(ctl_df)} 行，netMHCpan报告包含 {len(mhc_df)} 行。")

    print("正在标准化两个文件中的Protein_ID列以确保格式一致...")
    
    if 'Protein_ID' in ctl_df.columns:
        ctl_df['Protein_ID'] = ctl_df['Protein_ID'].str.replace(' ', '', regex=False).str.replace('.', '_', regex=False)
    else:
        print(f"警告: {os.path.basename(ctl_file)} 中缺少 'Protein_ID' 列。")

    if 'Identity' in mhc_df.columns:
        mhc_df['Identity'] = mhc_df['Identity'].str.replace(' ', '', regex=False).str.replace('.', '_', regex=False)
    else:
        print(f"警告: {os.path.basename(mhc_file)} 中缺少 'Identity' 列。")

    print("正在重命名netMHCpan的列以进行匹配: 'Identity' -> 'Protein_ID', 'MHC' -> 'HLA_Allele'")
    mhc_df_renamed = mhc_df.rename(columns={
        'Identity': 'Protein_ID',
        'MHC': 'HLA_Allele'
    })
    

    print("正在根据 Protein_ID, HLA_Allele, 和 Peptide 进行合并...")
    common_epitopes_df = pd.merge(
        ctl_df, 
        mhc_df_renamed, 
        on=['Protein_ID', 'HLA_Allele', 'Peptide'], 
        how='inner'
    )
    

    if common_epitopes_df.empty:
        print("完成，但未发现任何同时被两个工具支持的抗原表位。")
    else:
        print(f"完成！发现 {len(common_epitopes_df)} 个共同支持的抗原表位。")
        
        os.makedirs(output_dir, exist_ok=True)
        common_epitopes_df.to_csv(output_file, index=False)
        print(f"整合报告已保存到: {output_file}")


def batch_process_main():
    """
    主函数，用于批量处理指定项目下的所有SRR样本。
    """
    target_project_id = "PRJNA1176011"
    base_search_dir = f"/home/yxqin/neoantigen/Analysis/pacbio/08_netCTLpan/parse/{target_project_id}"
    
    print(f"开始在目录中搜索已处理的样本: {base_search_dir}")
    if not os.path.isdir(base_search_dir):
        print(f"错误: 项目目录不存在: {base_search_dir}")
        return
        
    search_pattern = os.path.join(base_search_dir, 'SRR*')
    sample_dirs = [d for d in glob.glob(search_pattern) if os.path.isdir(d)]

    if not sample_dirs:
        print(f"警告: 在目录 {base_search_dir} 下未找到任何SRR样本文件夹。")
        return

    print(f"发现 {len(sample_dirs)} 个样本。开始处理...")
    
    for sample_path in sample_dirs:
        srr_id = os.path.basename(sample_path)
        print(f"\n==================================================")
        print(f"开始整合样本: {srr_id} (项目: {target_project_id})")
        print(f"==================================================")
        
        find_common_epitopes(target_project_id, srr_id)

    print("\n所有样本整合完毕！")


if __name__ == '__main__':
    batch_process_main()