import sys
import vcf

def generate_sample_vcfs(input_file, reference_sample):
    # 读取原始VCF文件
    vcf_reader = vcf.Reader(open(input_file, 'r'))
    
    # 创建新的VCF写入器
    vcf_writers = {}
    for sample in vcf_reader.samples:
        if sample != reference_sample:
            vcf_writers[sample] = vcf.Writer(open(f'{sample}.vcf', 'w'), vcf_reader)

    # 遍历变异位点
    for record in vcf_reader:
        reference_genotype = record.genotype(reference_sample)['GT']

        # 与参考样本不同的变异位点添加到相应的样本的新VCF文件
        for sample in vcf_reader.samples:
            if sample != reference_sample:
                genotype = record.genotype(sample)['GT']

                # 如果与参考样本不同，将变异位点添加到相应的样本的新VCF文件
                if genotype != reference_genotype:
                    vcf_writers[sample].write_record(record)

    # 关闭所有VCF写入器
    for writer in vcf_writers.values():
        writer.close()

# 从命令行获取输入参数
if len(sys.argv) < 3:
    print("Error: Please provide the input VCF filename and reference sample name.")
    print("Usage: python script.py <input_file.vcf> <reference_sample>")
else:
    input_file = sys.argv[1]
    reference_sample = sys.argv[2]

    # 调用函数生成样本的新VCF文件
    generate_sample_vcfs(input_file, reference_sample)
