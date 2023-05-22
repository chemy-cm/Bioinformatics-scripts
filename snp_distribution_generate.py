import sys
import pandas as pd

def calculate_snp_distribution(vcf_file, bin_size, genome_info_file, output_file):
    # 读取VCF文件
    vcf_data = pd.read_csv(vcf_file, sep='\t', comment='#', header=None)

    # 获取染色体列表
    chromosomes = vcf_data[0].unique()

    # 读取基因组信息表格
    genome_info = pd.read_csv(genome_info_file, sep='\t', header=None)
    genome_info.columns = ['Chromosome', 'Length']

    # 计算每个染色体的起始位置
    genome_info['Start'] = genome_info['Length'].cumsum() - genome_info['Length'] + 1

    # 初始化结果数据框架
    result = pd.DataFrame(columns=['CHROM', 'BIN_START', 'BIN_END', 'SNP_COUNT'])

    # 遍历每个染色体
    for chromosome in chromosomes:
        # 获取染色体的长度
        chrom_length = genome_info.loc[genome_info['Chromosome'] == chromosome, 'Length'].values[0]

        # 计算染色体的BIN数目
        bin_count = int(chrom_length / bin_size) + 1

        # 初始化BIN起始位置和结束位置
        bin_starts = [i * bin_size for i in range(bin_count)]
        bin_ends = [(i + 1) * bin_size - 1 for i in range(bin_count)]

        # 根据VCF文件计算每个BIN中的SNP数目
        snp_counts = []
        for start, end in zip(bin_starts, bin_ends):
            snp_count = ((vcf_data[0] == chromosome) & (vcf_data[1] >= start) & (vcf_data[1] <= end)).sum()
            snp_counts.append(snp_count)

        # 构建当前染色体的结果数据框
        chrom_result = pd.DataFrame({
            'CHROM': [chromosome] * bin_count,
            'BIN_START': bin_starts,
            'BIN_END': bin_ends,
            'SNP_COUNT': snp_counts
        })

        # 将当前染色体的结果添加到总结果中
        result = result.append(chrom_result)

    # 根据染色体的起始位置对BIN_START进行偏移
    for chromosome in chromosomes:
        start_offset = genome_info.loc[genome_info['Chromosome'] == chromosome, 'Start'].values[0]
        result.loc[result['CHROM'] == chromosome, 'BIN_START'] += start_offset

    # 将结果写入输出文件
    result.to_csv(output_file, sep='\t', index=False)

def print_help():
    print("Usage: python snp_distribution.py [vcf_file] [bin_size] [genome_info_file] [output_file]")
    print("Arguments:")
    print("  vcf_file            Input VCF file")
    print("  bin_size            Bin size for SNP distribution")
    print("  genome_info_file    Input genome information file")
    print("  output_file         Output file for SNP distribution results")

if __name__ == '__main__':
    # 如果提供了 -h 或 --help 参数，显示帮助信息
    if "-h" in sys.argv or "--help" in sys.argv:
        print_help()
    else:
        # 检查输入参数是否正确
        if len(sys.argv) != 5:
            print("Error: Invalid number of arguments!")
            print_help()
        else:
            # 从命令行参数获取输入
            vcf_file = sys.argv[1]
            bin_size = int(sys.argv[2])
            genome_info_file = sys.argv[3]
            output_file = sys.argv[4]

            # 调用函数生成SNP分布结果
            calculate_snp_distribution(vcf_file, bin_size, genome_info_file, output_file)
