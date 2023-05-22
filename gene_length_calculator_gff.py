#!/usr/bin/env python

import sys

def calculate_gene_lengths(gff_file, output_file):
    gene_lengths = {}
    
    with open(gff_file, 'r') as gff:
        for line in gff:
            if not line.startswith('#'):
                fields = line.strip().split('\t')
                if len(fields) >= 9 and fields[2] == 'gene':
                    start = int(fields[3])
                    end = int(fields[4])
                    gene_name = fields[8].split(';')[0].split('=')[1]
                    gene_lengths[gene_name] = end - start + 1
    
    with open(output_file, 'w') as output:
        for gene_name, length in gene_lengths.items():
            output.write(f'{gene_name}\t{length}\n')

if __name__ == '__main__':
    if len(sys.argv) != 3:
        print('Usage: python gene_length_calculator.py <gff_file> <output_file>')
    else:
        gff_file = sys.argv[1]
        output_file = sys.argv[2]
        calculate_gene_lengths(gff_file, output_file)
