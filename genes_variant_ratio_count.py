#!/usr/bin/env python

import sys

def parse_vcf_file(vcf_file):
    gene_counts = {}
    with open(vcf_file, 'r') as f:
        for line in f:
            if not line.startswith('#'):
                columns = line.strip().split('\t')
                info = columns[7]
                info_fields = info.split(';')
                if len(info_fields) >= 15:
                    gene_info = info_fields[14]
                    gene_name = gene_info.split('|')[4]
                    if gene_name in gene_counts:
                        gene_counts[gene_name] += 1
                    else:
                        gene_counts[gene_name] = 1
    return gene_counts

def calculate_mutation_rate(gene_counts, gene_length_file):
    gene_lengths = {}
    with open(gene_length_file, 'r') as f:
        for line in f:
            gene, length = line.strip().split('\t')
            gene_lengths[gene] = int(length)
    
    result = []
    for gene, count in gene_counts.items():
        length = gene_lengths.get(gene, 0)
        mutation_rate = count / length if length > 0 else 0
        result.append((gene, count, length, mutation_rate))
    
    return result

if __name__ == '__main__':
    if len(sys.argv) < 3:
        print("Usage: python script.py <vcf_file> <gene_length_file>")
        sys.exit(1)
    
    vcf_file = sys.argv[1]
    gene_length_file = sys.argv[2]
    
    gene_counts = parse_vcf_file(vcf_file)
    result = calculate_mutation_rate(gene_counts, gene_length_file)
    
    print("Gene Name\tSNPs Count\tGene Length\tMutation Rate")
    for gene, count, length, mutation_rate in result:
        print(f"{gene}\t{count}\t{length}\t{mutation_rate}")
