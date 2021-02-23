#!/usr/bin/python3   
import sys
import os
import statistics

def mean_completeness(input_file):

    threads = 20
    filename = os.path.splitext(input_file)[0]
    os.system (f"export CHECKVDB=/Nancy/mrayko/Libs/checkv_db/checkv-db-v0.6")
    os.system (f"/home/mrayko/miniconda3/bin/checkv end_to_end {input_file} checkv_{filename}  -t {threads}")

    completeness = []
    for line in open (f"checkv_{filename}/quality_summary.tsv"):
        if is_number(line.split()[8].split()[0]):
            completeness.append(float(line.split()[8].split()[0]))
    print (statistics.mean(completeness))
    with open(f"checkv_{filename}/mean.txt", 'w') as f:
        f.write('%f' % statistics.mean(completeness))

def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print ("Usage: " + sys.argv[0] + "<fasta file with viral sequences>")
        exit()
    infile = sys.argv[1]
    mean_completeness(infile)
