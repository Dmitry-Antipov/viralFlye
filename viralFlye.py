#!/usr/bin/env python3


import sys
import os
import random
import time
import datetime
from os import listdir
from os.path import isfile, join, exists
from extract_circulars_metaflye import extract_circulars, extract_linears
from extract_small_comps_gfa import  extract_paths_in_components
from pb_download import download_sra_from_list
from check_linear_viruses import *
import argparse


#Not used
pbclip = "/home/dantipov/other_tools/pbclip/pbclip"
metaflye = "/home/dantipov/other_tools/Flye/bin/flye"



#freebayes polishing, currently not used
#bwa index ../assembly.fasta; bwa mem -t 20 ../assembly.fasta /Bmo/dantipov/data/metaFlye_viral/Aloha/SRR8811963_1.fastq.gz /Bmo/dantipov/data/metaFlye_viral/Aloha/SRR8811963_2.fastq.gz -t 20 | samtools sort -@8 -o assembly.bam
#./freebayes-parallel <(./fasta_generate_regions.py /Iceking/dantipov/metaFlye/Aloha/SRR8811961_rr/assembly.fasta.fai 500000) 16 -f /Iceking/dantipov/metaFlye/Aloha/SRR8811961_rr/assembly.fasta  /Iceking/dantipov/metaFlye/Aloha/SRR8811961_rr/pilon/assembly.bam > /Iceking/dantipov/metaFlye/Aloha/SRR8811961_rr/freebayes/assembly.vcf
#samtools faidx ../assembly.fasta; samtools index assembly.bam; freebayes-parallel <(fasta_generate_regions.py ../assembly.fasta.fai 500000) 16 -f ../assembly.fasta assembly.bam > assembly.vcf
#bgzip assembly.vcf;  ~/other_tools/bcftools/bcftools index assembly.vcf.gz;  ~/other_tools/bcftools/bcftools consensus -f ../assembly.fasta -o assembly.fasta assembly.vcf.gz

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

    
    if len(sys.argv)==1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    return parser.parse_args()



def run_circular_vv (args):
        fullpath = args.dir
        if os.path.isdir(fullpath):
            contigs_all = join(fullpath, "assembly.fasta")
            outdir = join(fullpath, "vv_circulars")
            if os.path.exists (contigs_all) and not os.path.exists(outdir):
                os.mkdir(outdir)
                stats = join(fullpath, "assembly_info.txt")
                circulars = join(fullpath, "circulars.txt")
                circ_fasta = join(fullpath, "circulars.fasta")
                extract_circulars(stats, circulars, int(args.min_viral_length))
                seqtk_line = (f'seqtk subseq {contigs_all} {circulars} > {circ_fasta}')
                print (seqtk_line)
                os.system(seqtk_line)
                vv_line =f"viralverify.py -f {circ_fasta} -o {outdir}  --hmm {args.hmm} -t 10"
                print (vv_line)
                os.system (vv_line)

def run_linear_vv (args):
        fullpath = args.dir
        hmms = args.hmm
        if os.path.isdir(fullpath):
            contigs_all = join(fullpath, "assembly.fasta")
            outdir = join(fullpath, "vv_linears")
            if os.path.exists (contigs_all):
# and not os.path.exists(join(outdir, "Prediction_results_fasta")):
                os.makedirs(outdir, exist_ok=True)
                stats = join(fullpath, "assembly_info.txt")
                linears = join(fullpath, "linears.txt")
                linears_fasta = join(fullpath, "linears.fasta")
                extract_linears(stats, linears, int (args.min_viral_length))
                seqtk_line = (f'seqtk subseq {contigs_all} {linears} > {linears_fasta}')
                print (seqtk_line)
                os.system(seqtk_line)
                vv_line = f"viralverify.py -f  {linears_fasta}  -o {outdir}  --hmm {hmms}  -t 10"
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
            all_contigs = join(fullpath, "assembly.fasta")
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
            vc_str = (f'viralcomplete.py -t 10 -thr 0.5 -f {fullpath} -o {outdir}')
            print(vc_str)
            os.system(vc_str)
#/Bmo/dantipov/tools/viralComplete/viralcomplete.py -thr 0.5 -f /Iceking/dantipov/metaFlye/japanese/MO1-2_clipped/vv_linears/Prediction_results_fasta/
#linears_virus.fasta -o MO1-2_0.5_check2/

    
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
            graph_all = join(fullpath, "assembly_graph.gfa")
            outdir = join(fullpath, "vv_components")
            if os.path.exists (graph_all):
# and not os.path.exists(outdir):
#                os.mkdir(outdir)
                comp_fasta = join(fullpath, "components.fasta")
                extract_paths_in_components(graph_all, 5000, 1000000, 10, comp_fasta)
                vv_line =f"viralverify.py -f {comp_fasta} -o {outdir}  --hmm {args.hmm} -t 10"
                print (vv_line)
                os.system (vv_line)

def runall(args):
#    download_and_run(sys.argv[1], sys.argv[2])

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
