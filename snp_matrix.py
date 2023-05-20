import pandas as pd
import sys

# 解析vcf文件，提取样本名和SNP信息
def parse_vcf(vcf_file):
    samples = []
    snps = []
    with open(vcf_file, 'r') as file:
        for line in file:
            if line.startswith("#CHROM"):
                samples = line.strip().split('\t')[9:]
                snps = [[] for _ in range(len(samples))]
            elif not line.startswith("#"):
                fields = line.strip().split('\t')
                for i, data in enumerate(fields[9:], start=0):
                    snps[i].append(data.split(':')[0])
    return samples, snps

# 计算两个样本之间的SNP数
def calculate_snp_count(sample1, sample2):
    count = 0
    for s1, s2 in zip(sample1, sample2):
        if s1 != '.' and s2 != '.' and s1 != s2:
            count += 1
    return count

# 生成SNP矩阵
def generate_snp_matrix(samples, snps):
    matrix = [[0 for _ in range(len(samples))] for _ in range(len(samples))]
    for i in range(len(samples)):
        for j in range(i+1, len(samples)):
            count = calculate_snp_count(snps[i], snps[j])
            matrix[i][j] = count
            matrix[j][i] = count
    return matrix

# 将矩阵保存为CSV文件
def save_matrix_to_csv(matrix, samples, filename):
    df = pd.DataFrame(matrix, index=samples, columns=samples)
    df.to_csv(filename, index_label='Sample')

# 主函数
def main():
    if len(sys.argv) != 2:
        print("请提供正确的输入参数。")
        print("用法: python script.py vcf_file")
        return

    vcf_file = sys.argv[1]
    csv_file = 'snp_matrix.csv'  # 保存的CSV文件名
    
    samples, snps = parse_vcf(vcf_file)
    matrix = generate_snp_matrix(samples, snps)
    save_matrix_to_csv(matrix, samples, csv_file)
    print(f"矩阵已保存到文件: {csv_file}")

if __name__ == '__main__':
    main()
