import pandas as pd
import os
import glob
import shlex

def extract_ids_and_run_seqkit(info_csv_path, faa_db_path, id_list_path, output_faa_path, seqkit_exe):
    """
    从info.csv文件中提取transcript_id，保存为列表，然后调用seqkit筛选序列。
    """
    print(f"--- 步骤 1: 从 {os.path.basename(info_csv_path)} 提取ID并运行 seqkit ---")
    
    if not os.path.exists(info_csv_path):
        print(f"警告: Info文件未找到，跳过: {info_csv_path}")
        return False
    if not os.path.exists(faa_db_path):
        print(f"警告: Fasta数据库未找到，跳过: {faa_db_path}")
        return False
        
    try:
        df = pd.read_csv(info_csv_path)
        if 'transcript_id' not in df.columns:
            print(f"错误: {info_csv_path} 中缺少 'transcript_id' 列。")
            return False
        
        unique_ids = df['transcript_id'].unique()
        
        if len(unique_ids) == 0:
            print(f"在 {os.path.basename(info_csv_path)} 中未找到任何ID，无需提取序列。")
            open(output_faa_path, 'w').close()
            return True 

        with open(id_list_path, 'w') as f:
            for item_id in unique_ids:
                f.write(f"{item_id}\n")
        print(f"提取了 {len(unique_ids)} 个不重复的ID到 {os.path.basename(id_list_path)}")

        cmd = (
            f"{shlex.quote(seqkit_exe)} grep --pattern-file {shlex.quote(id_list_path)} "
            f"{shlex.quote(faa_db_path)} -o {shlex.quote(output_faa_path)}"
        )
        
        print(f"执行seqkit命令: {cmd}")
        exit_code = os.system(cmd)
        
        if exit_code != 0:
            print(f"错误: seqkit 执行失败，退出码: {exit_code}")
            return False
            
        print("seqkit 成功完成。")
        return True

    except Exception as e:
        print(f"处理 {os.path.basename(info_csv_path)} 时发生错误: {e}")
        return False

def run_netmhcpan(netmhcpan_exe, hla_string, input_faa, output_netmhcpan):
    """
    构建并执行netMHCpan命令。
    """
    print(f"--- 步骤 2: 运行 netMHCpan ---")

    if not os.path.exists(input_faa) or os.path.getsize(input_faa) == 0:
        print(f"警告: 输入的fasta文件为空或不存在: {os.path.basename(input_faa)}。跳过netMHCpan分析。")
        open(output_netmhcpan, 'w').close()
        return

    cmd = (
        f"{shlex.quote(netmhcpan_exe)} -a {hla_string} -BA -l 9 "
        f"-f {shlex.quote(input_faa)} > {shlex.quote(output_netmhcpan)}"
    )

    print(f"执行netMHCpan命令: {cmd}")
    exit_code = os.system(cmd)
    
    if exit_code == 0:
        print(f"netMHCpan 成功完成。结果保存在: {os.path.basename(output_netmhcpan)}")
    else:
        print(f"错误: netMHCpan 执行失败，退出码: {exit_code}。")


def process_single_srr(project_id, srr_id, config):
    """
    为单个SRR样本执行完整的ID提取、seqkit和netMHCpan流程。
    """
    novel_info_csv = f"/home/yxqin/neoantigen/Analysis/pacbio/06_netMHCpan/data/{project_id}/{srr_id}/{srr_id}_use_netMHCpan.novel.info.csv"
    known_info_csv = f"/home/yxqin/neoantigen/Analysis/pacbio/06_netMHCpan/data/{project_id}/{srr_id}/{srr_id}_use_netMHCpan.known.info.csv"
    faa_database = f"/home/yxqin/neoantigen/Analysis/pacbio/05_repeatmasker/{project_id}/{srr_id}/{srr_id}_final_filtered.faa"
    
    temp_dir = f"/home/yxqin/neoantigen/Analysis/pacbio/06_netMHCpan/data/tmp/{project_id}/{srr_id}"
    output_dir = f"/home/yxqin/neoantigen/Analysis/pacbio/06_netMHCpan/output/{project_id}/{srr_id}"
    
    os.makedirs(temp_dir, exist_ok=True)
    os.makedirs(output_dir, exist_ok=True)
    
    # --- 处理 "Novel" 数据 ---
    print("\n--- 处理 Novel 数据 ---")
    novel_id_list = os.path.join(temp_dir, f"{srr_id}_novel_ids.txt")
    novel_faa_out = os.path.join(temp_dir, f"{srr_id}_use_netMHCpan_novel.faa")
    novel_netmhc_out = os.path.join(output_dir, f"{srr_id}_netMHCpan.novel.out")
    
    if extract_ids_and_run_seqkit(novel_info_csv, faa_database, novel_id_list, novel_faa_out, config['seqkit_exe']):
        run_netmhcpan(config['netmhcpan_exe'], config['hla_string'], novel_faa_out, novel_netmhc_out)

    # --- 处理 "Known" 数据 ---
    print("\n--- 处理 Known 数据 ---")
    known_id_list = os.path.join(temp_dir, f"{srr_id}_known_ids.txt")
    known_faa_out = os.path.join(temp_dir, f"{srr_id}_use_netMHCpan_known.faa")
    known_netmhc_out = os.path.join(output_dir, f"{srr_id}_netMHCpan.known.out")

    if extract_ids_and_run_seqkit(known_info_csv, faa_database, known_id_list, known_faa_out, config['seqkit_exe']):
        run_netmhcpan(config['netmhcpan_exe'], config['hla_string'], known_faa_out, known_netmhc_out)


def batch_processor_main():
    """
    主函数，用于批量处理指定项目下的所有SRR样本。
    """
    # ======================================================================
    # ===> 在这里指定想处理的PRJ项目ID <===
    target_project_id = "PRJNA1176011"
    # ======================================================================
    
    config = {
        'seqkit_exe': '/home/yxqin/software/seqkit/seqkit',
        'netmhcpan_exe': '/home/yxqin/software/netMHCpan/netMHCpan-4.1/netMHCpan',
        'hla_file': '/home/yxqin/neoantigen/Analysis/pacbio/06_netMHCpan/use_netMHCpan_MHCIid'
    }

    # 1. 读取HLA文件并格式化为netMHCpan所需的字符串
    try:
        with open(config['hla_file'], 'r') as f:
            hla_types = [line.strip() for line in f if line.strip()]
        config['hla_string'] = ",".join(hla_types)
        print(f"成功加载并格式化 {len(hla_types)} 个HLA类型。")
    except FileNotFoundError:
        print(f"错误: HLA文件未找到: {config['hla_file']}")
        return

    # 2. 查找并处理指定项目下的所有SRR样本
    input_data_dir = f"/home/yxqin/neoantigen/Analysis/pacbio/06_netMHCpan/data/{target_project_id}"
    
    print(f"\n开始在目录中搜索已处理的样本: {input_data_dir}")
    if not os.path.isdir(input_data_dir):
        print(f"错误: 项目目录不存在: {input_data_dir}")
        return
        
    search_pattern = os.path.join(input_data_dir, 'SRR*')
    sample_dirs = [d for d in glob.glob(search_pattern) if os.path.isdir(d)]

    if not sample_dirs:
        print(f"警告: 在目录 {input_data_dir} 下未找到任何SRR样本文件夹。")
        return

    print(f"发现 {len(sample_dirs)} 个样本。开始处理...")
    
    for sample_path in sample_dirs:
        srr_id = os.path.basename(sample_path)
        print(f"\n==================================================")
        print(f"正在处理样本: {srr_id} (项目: {target_project_id})")
        print(f"==================================================")
        
        process_single_srr(target_project_id, srr_id, config)

    print("\n所有样本处理完毕！")


if __name__ == '__main__':
    batch_processor_main()
