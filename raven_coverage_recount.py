#!/usr/bin/env python3
import argparse
import os
import sys

def parse_args(args):
###### Command Line Argument Parser
    parser = argparse.ArgumentParser(description = "recounting coverage for raven with average read length in the dataset",
    formatter_class=argparse.RawTextHelpFormatter)

    parser._action_groups.pop()
    required_args = parser.add_argument_group('required arguments')
    required_args.add_argument('--contigs', required = True, help='raven contigs file')
    required_args.add_argument('--filtered_contigs', required = True, help='filtered contigs file')   
    required_args.add_argument('--reads', help='Path to long reads')    
    optional_args = parser.add_argument_group('optional arguments')
    optional_args.add_argument('--coverage_cutoff', default = 10, help = 'coverage cutoff, default 10 ')    
    if len(sys.argv)==1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    return parser.parse_args()


def average_readlength(f):
    from Bio import SeqIO
    format = "fasta"
    arr = os.path.splitext(f)
    period = 2
    ext = arr[-1]
    gzipped = False
    if ext == ".gz":
        ext = os.path.splitext(arr[-2])[-1]
        gzipped = True
    if ext == ".fastq" or ext == ".fq":    
        format = "fastq"
    elif ext == ".fasta" or ext == ".fa":    
        format = "fasta"
    else:
        print("unknown format " + ext)
        exit()
    read_count = 0
    read_l = 0
    if gzipped:
        import gzip
        handle = gzip.open(f, "rt")
    else:
        handle = open(f, "rt")   
    for record in SeqIO.parse(handle, format):
        read_count +=1
        read_l += len(record.seq)
    return read_l/read_count

if __name__ == "__main__":
    args = parse_args(sys.argv[1:])
    rl = average_readlength(args.reads)
#    rl = 3300.0
    print (f'average read length: {rl}')
    count = 0
    outcontigs = open(args.filtered_contigs, 'w')
    to_out = False
#>Utg61358 LN:i:45259 RC:i:17 XO:i:1
    for line in open (args.contigs, 'r'):
        if line[0] == ">":
            cov = float (line.split()[2].split(':')[2])
            length = float(line.split()[1].split(':')[2])
            recounted_cov = cov * rl / length
            print (f' old cov: {cov}, len: {length}, new cov:{recounted_cov}')

            if recounted_cov > float(args.coverage_cutoff):
                to_out = True
            else:
                to_out = False
        if to_out:
            outcontigs.write(line)

            
