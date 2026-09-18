# 功能: 在运行netCTLpan前，预处理输入的fasta文件，
#       将header中的制表符(tab)转换为空格，并保存为新文件。

import os
import glob
import shlex

def create_tab_fixed_fasta(original_faa_path, fixed_faa_path):
    """
    读取原始fasta文件，修复header中的制表符，并写入一个新的文件。
    返回新文件的路径。
    """
    if not os.path.exists(original_faa_path) or os.path.getsize(original_faa_path) == 0:
        return None

    try:
        with open(original_faa_path, 'r') as f_in, open(fixed_faa_path, 'w') as f_out:
            for line in f_in:
                if line.startswith('>'):
                    f_out.write(line.replace('\t', ' '))
                else:
                    f_out.write(line)
        
        print(f"  -> 已创建tab修复后的fasta文件: {os.path.basename(fixed_faa_path)}")
        return fixed_faa_path

    except Exception as e:
        print(f"警告: 在修复fasta文件 {os.path.basename(original_faa_path)} 时发生错误: {e}")
        return None


def run_netctlpan(netctlpan_exe, hla_string, input_faa, output_netctlpan):
    """
    构建并直接执行netCTLpan命令。
    """
    print(f"--- 运行 netCTLpan ---")

    if not os.path.exists(input_faa) or os.path.getsize(input_faa) == 0:
        print(f"警告: 输入的fasta文件为空或不存在: {os.path.basename(input_faa)}。跳过netCTLpan分析。")
        open(output_netctlpan, 'w').close()
        return

    cmd = (
        f"{shlex.quote(netctlpan_exe)} -a {hla_string} "
        f"-f {shlex.quote(input_faa)} > {shlex.quote(output_netctlpan)}"
    )

    print(f"执行netCTLpan命令: {cmd}")
    exit_code = os.system(cmd)
    
    if exit_code == 0:
        print(f"netCTLpan 成功完成。结果保存在: {os.path.basename(output_netctlpan)}")
    else:
        print(f"错误: netCTLpan 执行失败，退出码: {exit_code}。")


def process_single_srr_for_netctlpan(project_id, srr_id, config):
    """
    为单个SRR样本执行完整的netCTLpan流程。
    """
    input_faa_dir = f"/home/yxqin/neoantigen/Analysis/pacbio/06_netMHCpan/data/tmp/{project_id}/{srr_id}"
    output_dir = f"/home/yxqin/neoantigen/Analysis/pacbio/08_netCTLpan/{project_id}/{srr_id}"
    
    os.makedirs(output_dir, exist_ok=True)
    
    # 循环处理 Novel 和 Known
    for data_type in ["novel", "known"]:
        print(f"\n--- 处理 {data_type.capitalize()} 数据 ---")
        
        original_input_faa_path = os.path.join(input_faa_dir, f"{srr_id}_use_netMHCpan_{data_type}.faa")
        
        fixed_input_faa_path = os.path.join(input_faa_dir, f"{srr_id}_use_netMHCpan_{data_type}.tabfixed.faa")
        
        create_tab_fixed_fasta(original_input_faa_path, fixed_input_faa_path)
        
        output_netctlpan_file = os.path.join(output_dir, f"{srr_id}_netCTLpan.{data_type}.out")
        
        run_netctlpan(
            config['netctlpan_exe'],
            config['hla_string'],
            fixed_input_faa_path,
            output_netctlpan_file
        )


def main():
    """
    主函数，用于批量处理指定项目下的所有SRR样本。
    """
    target_project_id = "PRJNA1176011"
    
    config = {
        'netctlpan_exe': '/home/yxqin/software/netCTLpan-1.1b/netCTLpan-1.1/netCTLpan',
        'hla_file': '/home/yxqin/neoantigen/Analysis/pacbio/06_netMHCpan/use_netMHCpan_MHCIid'
    }

    try:
        with open(config['hla_file'], 'r') as f:
            hla_types = [line.strip() for line in f if line.strip()]
        config['hla_string'] = ",".join(hla_types)
        print(f"成功加载并格式化 {len(hla_types)} 个HLA类型。")
    except FileNotFoundError:
        print(f"错误: HLA文件未找到: {config['hla_file']}")
        return

    base_dir_to_scan = f"/home/yxqin/neoantigen/Analysis/pacbio/06_netMHCpan/data/tmp/{target_project_id}"
    
    print(f"\n开始在目录中搜索已处理的样本: {base_dir_to_scan}")
    if not os.path.isdir(base_dir_to_scan):
        print(f"错误: 用于发现样本的项目目录不存在: {base_dir_to_scan}")
        return
        
    search_pattern = os.path.join(base_dir_to_scan, 'SRR*')
    sample_dirs = sorted([d for d in glob.glob(search_pattern) if os.path.isdir(d)])

    if not sample_dirs:
        print(f"警告: 在目录 {base_dir_to_scan} 下未找到任何SRR样本文件夹。")
        return

    print(f"发现 {len(sample_dirs)} 个样本。开始处理...")
    
    for sample_path in sample_dirs:
        srr_id = os.path.basename(sample_path)
        print(f"\n==================================================")
        print(f"正在处理样本: {srr_id} (项目: {target_project_id})")
        print(f"==================================================")
        
        process_single_srr_for_netctlpan(target_project_id, srr_id, config)

    print("\n所有样本处理完毕！")


if __name__ == '__main__':
    main()