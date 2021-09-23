#!/usr/bin/env python3


import sys
import os
import random
import time
import datetime
from os import listdir
from os.path import isfile, join, exists
from extract_circulars_metaflye import extract_circulars, extract_linears, extract_circulars_raven, extract_linears_raven
from extract_small_comps_gfa import  extract_paths_in_components
from pb_download import download_sra_from_list
from check_linear_viruses import prepare_index_and_depth
import argparse


#Not used
pbclip = "/home/dantipov/other_tools/pbclip/pbclip"
metaflye = "/home/dantipov/other_tools/Flye/bin/flye"




def parse_args(args):
###### Command Line Argument Parser
    parser = argparse.ArgumentParser(description=f'''Wrapper script for viralFlye pipeline \n
Usage: {sys.argv[0]} <metaflye_output_folder>
Circular complete viruses can be found in vc_circ and vc_linears subfolders in metaflye_output_folder
Viral components can be found in vv_components subfolder
All viruses (including incomplete) can be found in vv_circulars and vv_linears
CircularDisconnector results are in vv_circulars/linear_check/ subfolder
Details can be found in viralverify and viralcomplete manual''',
    formatter_class=argparse.RawTextHelpFormatter)
    parser._action_groups.pop()
    required_args = parser.add_argument_group('required arguments')
    required_args.add_argument('--dir', required = True, help='metaFlye output directory')
    required_args.add_argument('--hmm', required = True, help='Path to Pfam-A HMM database for viralVerify script')    
    required_args.add_argument('--reads', help='Path to long reads')    
    optional_args = parser.add_argument_group('optional arguments')
    optional_args.add_argument('--min_viral_length', default = 5000, help = 'minimal limit on the viral length under study, default 5k')    
    optional_args.add_argument('--ill1', default = '', help = "file with left illumina reads for polishing")
    optional_args.add_argument('--ill2', default = '', help = "file with right illumina reads for polishing")
    parser.add_argument('--raven', dest='raven', action='store_true')

    parser.set_defaults(raven=False)
    
    if len(sys.argv)==1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    return parser.parse_args()



def run_circular_vv (args):
        fullpath = args.dir
        if os.path.isdir(fullpath):
            contigs_all = args.assembly
            outdir = join(fullpath, "vv_circulars")
            if os.path.exists (contigs_all):
                os.makedirs(outdir, exist_ok=True)
                stats = join(fullpath, "assembly_info.txt")
                circulars = join(fullpath, "circulars.txt")
                circ_fasta = join(fullpath, "circulars.fasta")
                if not args.raven:
                    extract_circulars(stats, circulars, int(args.min_viral_length))
                else:
                    extract_circulars_raven(contigs_all, circulars, int(args.min_viral_length),args.rl)

                seqtk_line = (f'seqtk subseq {contigs_all} {circulars} > {circ_fasta}')
                print (seqtk_line)
                os.system(seqtk_line)
                vv_line =f"viralverify -f {circ_fasta} -o {outdir}  --hmm {args.hmm} -t 10"
                print (vv_line)
                os.system (vv_line)

def run_linear_vv (args):
        fullpath = args.dir
        hmms = args.hmm
        if os.path.isdir(fullpath):
            contigs_all = args.assembly
            outdir = join(fullpath, "vv_linears")
            if os.path.exists (contigs_all):
# and not os.path.exists(join(outdir, "Prediction_results_fasta")):
                os.makedirs(outdir, exist_ok=True)
                stats = join(fullpath, "assembly_info.txt")
                linears = join(fullpath, "linears.txt")
                linears_fasta = join(fullpath, "linears.fasta")
                if not args.raven:
                    extract_linears(stats, linears, int (args.min_viral_length))
                else:
                    extract_linears_raven(contigs_all, args.graph, linears, int (args.min_viral_length), args.rl)

                seqtk_line = (f'seqtk subseq {contigs_all} {linears} > {linears_fasta}')
                print (seqtk_line)
                os.system(seqtk_line)
                vv_line = f"viralverify -f  {linears_fasta}  -o {outdir}  --hmm {hmms}  -t 10"
                print (vv_line)
                os.system (vv_line)

def run_linear_check(args):
        fullpath = args.dir
        reads = args.reads
#     for dir in listdir(indir):
#        fullpath = join(indir, dir)    
        if os.path.isdir(fullpath):
            circular_vv = join(fullpath, "vv_circulars", "Prediction_results_fasta")
            circ_fasta = join(circular_vv, "circulars_virus.fasta")
            linear_check_res = join(circular_vv, "linear_check")
            all_contigs = args.assembly
#            reads = join(sys.argv[2], )
            if os.path.exists(circ_fasta):
                prepare_index_and_depth(circ_fasta, all_contigs, args.reads, os.path.join(os.path.dirname(circ_fasta), "linear_check"))

             #   linear_check_line = (f'./check_linear_viruses.py {circ_fasta} {args.reads} {linear_check_res}')
             #   print (linear_check_line)
             #   os.system(linear_check_line)

def run_vc(args, pref, name):
        indir = args.dir
#     for dir in listdir(indir):
        fullpath = join(indir, "vv_" + pref, "Prediction_results_fasta",name +"_virus.fasta")  
        outdir = join(indir, "vc_" + pref)    
        if os.path.exists(fullpath):
            vc_str = (f'./viralComplete/viralcomplete.py -t 10 -thr 0.5 -f {fullpath} -o {outdir}')
            print(vc_str)
            os.system(vc_str)
#/Bmo/dantipov/tools/viralComplete/viralcomplete.py -thr 0.5 -f /Iceking/dantipov/metaFlye/japanese/MO1-2_clipped/vv_linears/Prediction_results_fasta/
#linears_virus.fasta -o MO1-2_0.5_check2/

#freebayes polishing, currently not used
#bwa index ../assembly.fasta; bwa mem -t 20 ../assembly.fasta /Bmo/dantipov/data/metaFlye_viral/Aloha/SRR8811963_1.fastq.gz /Bmo/dantipov/data/metaFlye_viral/Aloha/SRR8811963_2.fastq.gz -t 20 | samtools sort -@8 -o assembly.bam
#./freebayes-parallel <(./fasta_generate_regions.py /Iceking/dantipov/metaFlye/Aloha/SRR8811961_rr/assembly.fasta.fai 500000) 16 -f /Iceking/dantipov/metaFlye/Aloha/SRR8811961_rr/assembly.fasta  /Iceking/dantipov/metaFlye/Aloha/SRR8811961_rr/pilon/assembly.bam > /Iceking/dantipov/metaFlye/Aloha/SRR8811961_rr/freebayes/assembly.vcf
#samtools faidx ../assembly.fasta; samtools index assembly.bam; freebayes-parallel <(fasta_generate_regions.py ../assembly.fasta.fai 500000) 16 -f ../assembly.fasta assembly.bam > assembly.vcf
#bgzip assembly.vcf;  ~/other_tools/bcftools/bcftools index assembly.vcf.gz;  ~/other_tools/bcftools/bcftools consensus -f ../assembly.fasta -o assembly.fasta assembly.vcf.gz

def run_freebayes(args):
    #bgzip samtools freebayes bcftools
    args.bam = join(args.dir, "assembly.bam")
    args.vcf = join(args.dir, "assembly.vcf")
    bwa_line = f'bwa index {args.assembly}; bwa mem -t 20 {args.assembly} {args.ill1} {args.ill2} -t 20 | samtools sort -@8 -o {args.bam}'
    print (bwa_line)
    os.system(bwa_line)
    samtools_line =f'samtools faidx {args.assembly}; samtools index {args.bam}'
    print(samtools_line)
    os.system(samtools_line)
    freebayes_line = f'freebayes-parallel <(./fasta_generate_regions.py {args.assembly}.fai 500000) 16 -f {args.assembly} {args.bam} > {args.vcf}'
    import subprocess
    subprocess.call(['bash', '-c', freebayes_line])
    print(freebayes_line)
 #   os.system(freebayes_line)  
    bcftools_line = f'bgzip {args.vcf} --force;  bcftools index {args.vcf}.gz;  bcftools consensus -f {args.assembly} -o {join(args.dir, "assembly.cor.fasta")} {args.vcf}.gz'
    print (bcftools_line)
    os.system(bcftools_line)
    args.assembly = join(args.dir, "assembly.cor.fasta")

    
def download_and_run(read_dirs, output_dirs):
    for dir in listdir(read_dirs):
#FAKO_1
        if dir.find("_") == -1 and dir != "illumina":
            sra_list = join(read_dirs, dir, 'sra.list')
            full_in_dir =  join(read_dirs, dir)
            if not exists(sra_list):
                print(f'No sra list in {dir}')
                continue
            out_d = dir + "_clipped"
            if not exists(join(output_dirs, out_d)):
                print(f'will metaflye {dir}')
                if not exists(join(read_dirs, dir, "clipped.fasta")):
                    print(f' will download and pbclip {dir}')
                    download_sra_from_list(sra_list, full_in_dir)
                    pb_clip_line = f'{pbclip} {full_in_dir}/*.gz > {full_in_dir}/clipped.fasta'
                    print (pb_clip_line)
                    os.system(pb_clip_line)
                metaflye_line = f'{metaflye} --meta --pacbio-raw {full_in_dir}/clipped.fasta --threads 30 --o {join(output_dirs,out_d)} --genome-size 500M'            
                print (metaflye_line)
                os.system (metaflye_line)
            
        
def run_on_components (args):
        fullpath = args.dir
        if os.path.isdir(fullpath):
            graph_all = args.graph
            outdir = join(fullpath, "vv_components")
            if os.path.exists (graph_all):
# and not os.path.exists(outdir):
#                os.mkdir(outdir)
                comp_fasta = join(fullpath, "components.fasta")
                extract_paths_in_components(graph_all, 5000, 1000000, 10, comp_fasta)
                vv_line =f"viralverify -f {comp_fasta} -o {outdir}  --hmm {args.hmm} -t 10"
                print (vv_line)
                os.system (vv_line)
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

def runall(args):
#    download_and_run(sys.argv[1], sys.argv[2])
#    print (args.raven)
    if (args.raven):
        args.rl = average_readlength(args.reads)
        print (args.rl)
    args.assembly = join(args.dir, "assembly.fasta")
    args.graph = join(args.dir, "assembly_graph.gfa")
    if args.ill1!= '':
        run_freebayes(args)
    
    run_linear_vv(args)
    run_circular_vv(args)  

    run_vc(args, "linears", "linears")
    run_vc(args, "circulars", "circulars")

    run_on_components(args)

    run_linear_check(args)


if __name__ == "__main__":
    args = parse_args(sys.argv[1:])
    if len(sys.argv) < 4:
        parser.print_help(sys.stderr)
        sys.exit(1)
    indir = sys.argv[1]
    runall(args)
#    for dir in os.listdir(indir):
#        runall(join(indir,dir))        
