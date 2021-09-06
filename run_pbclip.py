#!/usr/bin/python3   
import sys
import os
import statistics
import joblib
from joblib import Parallel, delayed

from os.path import join

def run_sample_pbclip(subdir, name):
    pbclip_path = "/home/dantipov/other_tools/pbclip/pbclip"
    pbclip_string = pbclip_path
    for f in (os.listdir(subdir)):
        if f.find("fastq.gz")!= -1:
            pbclip_string += " " + join(subdir, f)
    pbclip_string += " > " + join(subdir, name)
    print (pbclip_string)
    os.system(pbclip_string)
def run_all_pbclip(dir, name):
    for subdir in os.listdir(dir):
        if os.path.isdir(join(dir, subdir)):
            run_sample_pbclip(join(dir, subdir), name)

def run_all_assemblies(reads, outdir, type):
#    type = 'nano-raw'
#./build/bin/canu -pacbio-raw /Nancy/data/input/Meta/infested_sheep_gut/SRR10963010.fastq.gz -maxThreads=20 -d /Iceking/dantipov/viralFlye/Canu/SRR10963010/ -genomeSize=1610M -p SRR10963010 maxInputCoverage=10000 corOutCoverage=10000 corMhapSensitivity=high corMinCoverage=0 redMemory=32 oeaMemory=32 batMemory=200
    canu_path = "/home/dantipov/other_tools/canu/build/bin/canu"
    metaFlye_path = "/home/dantipov/other_tools/flye2.8.3/Flye-flye/bin/flye"
    metaFlye_V_path = "/home/dantipov/other_tools/Flye/bin/flye"
    canu_run_str = f'{canu_path} -{type} {reads}  -maxThreads=8 -d {join(outdir, "canu")} -genomeSize=500M -p canu_asm maxInputCoverage=10000 corOutCoverage=10000 corMhapSensitivity=high corMinCoverage=0 redMemory=32 oeaMemory=32 batMemory=200'
    print (canu_run_str)     
#--meta --pacbio-raw /Iceking/dantipov/data/japanese/ES1-2/clipped.fasta /Iceking/dantipov/data/japanese/ES9-1/clipped.fasta --threads 30 --o /Iceking/dantipov/metaFlye/japanese/ES_all_clipped --genome-size 500M
    flye_suff = f' --meta --{type} {reads} --threads 5 --o '
    os.makedirs(join(outdir, "metaflye"), exist_ok=True)
    os.makedirs(join(outdir, "metaflye_v"), exist_ok=True)
    
    metaFlye_str = metaFlye_path + flye_suff + join(outdir, "metaflye")

    metaFlye_V_str = metaFlye_V_path + flye_suff + join(outdir, "metaflye_v")
    print(metaFlye_str)
    print (metaFlye_V_str)
#    os.system(canu_run_str)
    os.system(metaFlye_str)
    os.system(metaFlye_V_str)
def run_dir_all_assemblies(dir, name, outdir):
    params = []
    for subdir in os.listdir(dir):
         if os.path.isdir(join(dir, subdir)):
            reads = join(dir, subdir, name)
            if os.path.exists(reads):
                outdir_p = join(outdir, subdir)
                params.append([reads, outdir_p])
    Parallel(n_jobs=5)(delayed(run_all_assemblies)(param[0], param[1], "pacbio-raw")
             for param in params)

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print ("Usage: " + sys.argv[0] + "<dir with samples (.fastq.gz)> <final name> <outdir>")
#        exit()

#    run_all_pbclip(sys.argv[1], sys.argv[2])
    run_dir_all_assemblies(sys.argv[1], sys.argv[2], sys.argv[3])
#    run_sample_pbclip(sys.argv[1],sys.argv[2])
    for line in open (sys.argv[1]):
        arr = line.split()
        run_all_assemblies(arr[0], arr[1], arr[2])
