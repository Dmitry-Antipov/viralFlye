import sys
import os
import random
import subprocess
import shutil
import re
from genericpath import isdir, exists
from os.path import join
from joblib import Parallel, delayed

all_refs=[]
mash_bin="/home/dantipov/other_tools/mash/mash"
plasmid_sketch = "/home/dantipov/other_tools/mash-Linux64-v2.0/plasmid_refseq.msh"
#bacterial_sketch = "/home/dantipov/other_tools/mash/plasmid_refseq.msh"
def parse_coords_filtered(threshold, inputfile, ref_len):
    with open(inputfile) as input:
        align_data = re.compile("\d+ \d+ \| \d+ \d+ \| \d+ (\d+) \| [\d.]+ \| [^ ]+ ([^\s]+length_(\d+)[^\s]+)")
        contig = None
        contig_len = 0
        align_len = 0
        total = 0
        def print_contig():
            per = 100.0 * align_len / contig_len
            if (per > threshold):
                print(contig, per)
                return align_len * 100.0 / ref_len
            else:
                return 0
        for line in input:
            res = align_data.search(line)
            if res is None:
                continue
            new_contig = res.group(2)
            if contig != new_contig:
                if contig is not None:
                    total += print_contig()
                contig = new_contig
                contig_len = int(res.group(3))
                align_len = 0
            #Assuming that all alignments of the same contig are consequent
            align_len += int(res.group(1))
        #Print the last contig separately
        total += print_contig()
        return total


             

def parse_quast(out_dir, name, check_coords=False):
    out_path = os.path.join(out_dir, "report.tsv")
    mis = 0
    fraction = 0
    ref_len = 0
    with open(out_path) as input:
        for line in input:
            arr = line.strip().split()
            if len(arr) < 2:
                continue
            if arr[0] == "Genome":
                fraction = float(arr[-1])
            if arr[0] == "Reference" and arr[1] == "length":
                ref_len = int(arr[-1])
            if arr[1] == "misassemblies":
                mis = int(arr[-1])
            if arr[0] == "Assembly":
                name_ass = arr[1].strip()
    if ref_len < 1000:
        print name + " is bad, ref_len: " + str(ref_len)
        return False
    elif fraction < 90:
        print name + " is bad, fraction " + str(fraction)
        return False
    elif mis > 5:
        print name + " is bad, misassemblies " + str(mis)
        return False
    else:
        print out_dir + "/contigs_reports/nucmer_output/" + name + ".coords.filtered"
        if check_coords:
            print name + " old fraction " + str(fraction)
            fraction = parse_coords_filtered(80, out_dir + "/contigs_reports/minimap_output/" + name_ass + ".coords.filtered", ref_len)
        if (fraction < 90):
            print name + " is bad, updated genome fraction "  + str(fraction)
            return False
        else:
            print name + " is good! updated genome fraction " + str(fraction)
            return True

def process_reference(args, r_dir):
    os.system(args[0])
    ref_dir = os.path.abspath(r_dir)
    out_dir = os.path.abspath(args[0].strip().split()[-1])
    if parse_quast(out_dir, args[1], True):
        shutil.copy2(args[1],ref_dir)

def filter_excessive(dir):
#    refs = os.listdir(dir)
    infos = []
    for line in all_refs:
        if os.path.exists(os.path.join(dir, line[1])):
            print line
            reff = os.path.join(dir,line[1])
            infos.append([reff, int(os.stat(reff).st_size), 1])
    for i in range (0, len (infos)):
        info = infos[i]
        for j in range(0, i):
            other_info = infos[j]
            if other_info[0] == info[0]:
                break
            if other_info[2] != 0 and other_info[1] < 1.05 * info[1] and info[1] < 1.05 * other_info[1]:
                quast_line= "quast.py -a all --ambiguity-usage 0.90 --fast -t 5 --min-contig 0 --min-alignment 0 -R " + other_info[0] + " " + info[0] + " > /dev/null -o tmp_ref_quast"
#                print quast_line
                os.system(quast_line)
                if parse_quast("tmp_ref_quast", info[0]) == True:
                    infos[i][2] = 0
                    break
    for i in range (0, len (infos)):
        if infos[i][2] == 0:
            os.system("rm " + infos[i][0])

def process_mash_output(mash_file, contigs,reference_dir):
    record_number = 0
    contig_file = os.path.abspath(contigs)
    quast_line = []
    id = 0;
    global all_refs
    with open(mash_file) as input:
        for line1 in input:
            path= line1.strip().split()[-1]
            ident = float(line1.strip().split()[0])
            all_refs.append ([ident, os.path.basename(path)])
            quast_line.append(["quast.py --fast -t 1 --min-contig 0 --min-identity 80 --min-alignment 0 -R " + path + " " + contig_file + " > tmp.log -o tmp_mash_quast" + str(id) , path])
            id += 1;
    all_refs.sort(reverse=True)
    Parallel(n_jobs=20)(delayed(process_reference) (args, reference_dir) for args in quast_line)
    os.system("rm -rf tmp_mash_quast*")

def run_mash(reads1, reads2, out):
    os.system(mash_bin + " screen -p 20  -w -i 0.95 " + plasmid_sketch + " " + reads1 + " " + reads2 + ">"+ out)


def main():
    if len(sys.argv) != 4:
        print "Usage: " + sys.argv[0] + " <mash output file> <contigs> <output dir>"
        exit()
    out_dir = os.path.abspath(sys.argv[3])
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)
    process_mash_output(sys.argv[1], sys.argv[2], sys.argv[3])
#    filter_excessive(out_dir)

if __name__ == '__main__':
    main()

