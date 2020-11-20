#!/usr/bin/python3
import sys
import os
import random
from scipy import signal
import pysam
import numpy as np
from matplotlib import pyplot as plt
'''
import matplotlib.pyplot as plt
import ruptures as rpt
import numpy as np
signal = np.genfromtxt("MK415399.depth", delimiter ="\t", usecols = 2)
signal10 = signal[::10]
signal100 = signal[::100]
algo = rpt.Pelt(model="l1").fit(signal100)
result = algo.predict(pen=200)
rpt.display(signal100, result)
plt.show()

algo = rpt.Window(model="l1", min_size=10, jump=10).fit(signal)
result = algo.predict(n_bkps = 5) # up to 5 breakpoints



all based on https://charles.doffy.net/files/sp-review-2020.pdf 
'''

jump = 1.3
window_size = 100
epsilon = 10
sampling_size = 20
#with -a option!!

def split_multifasta(infile, outdir):
    with open(infile) as input:
        line1 = input.next()
        while len(line1) > 0:
            curname = ""
            if len(line1) > 0 and line1[0] == '>':
                curname = line1.split()[0][1:] + ".fasta"
                out = open(os.path.join(outdir, curname), "w")
                while True:
                    out.write(line1)
                    line1 = input.next()
                    if len(line1) == 0 or line1[0] == '>':
                        break


def extract_suspicious_depth(depth_file):
    depth = []
    genome_len = 0
    for line in open(depth_file, 'r'):
        arr = line.split()
        depth.append(float(arr[2]))
        genome_len += 1
    genome_len = len(depth)
    ratios = []
    reverse_r = []
#Yes, this can be 20x faster)
    for i in range(window_size,  genome_len - window_size):
#1 for zeroes
        left = sum(depth[i - window_size:i]) + 1
        right = sum(depth[i+1: i + window_size + 1]) + 1
        ratios.append (left/right)
        reverse_r.append(right/left)
    peaks = signal.find_peaks(ratios, distance =2*window_size)[0]
    reverse_peaks = signal.find_peaks(reverse_r, distance = 2*window_size)[0]
    good_left = []
    good_right = []
    for i in peaks:
        if (ratios[i] > jump):
           print (f'{i + window_size} {ratios[i]} right')
           good_right.append(i + window_size)
    for i in reverse_peaks:
        if (ratios[i] < 1/jump):
           print (f'{i + window_size} {ratios[i]} left')
           good_left.append(i + window_size)
    if len(good_right) < 1 or len(good_left) < 1:
        print('suspicious repeats not found')
        return [-1, -1, genome_len]
    elif len(good_right) > 1 or len(good_left) > 1:
        print('multiple suspicious repeats found')
        return [-1, -1, genome_len]
    else:
        left = good_left[0]
        right = good_right[0]
        if (right - left < 100) or (right - left > 10000):
            print (f'distance is not similar to a repeat {right - left}')
            return [-1, -1, genome_len]
        else:
            return [left, right, genome_len]


#    print (peaks)
#    print (reverse_peaks)

    return []


def count_traversing_fraction (start_pos, end_pos, bamfile):
    traversing = 0
#TODO: name
    for read in bamfile.fetch('MK415399.1', start_pos, end_pos):
        if start_pos - epsilon > read.reference_start and end_pos + epsilon < read.reference_end :
            traversing += 1


    count = 0
    for read in bamfile.fetch('MK415399.1', (start_pos + end_pos) / 2, (start_pos + end_pos + 2) / 2 ):
        count += 1
    print (f' traversing {traversing} of {count}')
    return [traversing, count]

def check (start_pos, end_pos, bam, genome_length):
    bamfile = pysam.AlignmentFile(bam, "rb")
    if genome_length < (end_pos - start_pos) * 2:
        print('Suspicious region longer than half of genome')
        return False
    [traversing, count] = count_traversing_fraction(start_pos, end_pos, bamfile)
    sample = 0
    random_traversing = 0
    random_count = 0
    print ("checking random samples")
    while sample < sampling_size:
        check_start = random.randint(0, genome_length - end_pos + start_pos)
        if check_start > start_pos - epsilon and check_start < end_pos + epsilon:
            continue
        check_end = check_start + end_pos - start_pos
        ttraverse, tcount = count_traversing_fraction(check_start, check_end, bamfile)
        random_traversing += ttraverse
        random_count += tcount
        sample += 1
#TODO what if repeat is longer than read length
    if traversing * random_count <= count * random_traversing * 10:
        print ("significantly less traversing , repeat confirmed")
        return True
    else:
        return False

#TODO multiple reads files, other options for minimap
def prepare_index_and_depth(circulars, reads, workdir):

    os.makedirs(workdir, exist_ok=True)
    split_multifasta(circulars, workdir)
    for contig in workdir:
        if contig.split('.')[-1] == "fasta":
            full_contig = join(workdir, contig)
            name = contig.split('.')[0]
            bam_file = {join(workdir, name + ".bam")}
            if not os.path.exists(bam_file):
                bam_line = f'minimap2 -x map-pb -a -t 30 --sam-hit-only {full_contig} {reads} | samtools sort -o {bam_file}'
                os.system(bam_line)
            
            depth_file = join(workdir, name +".depth")
            depth_line = f'samtools depth -a {bam_file} >  {depth_file}'
            os.system(depth_line)
            [left, right, genome_len] = extract_suspicious_depth(depth_file)
            if left !=  -1:
                check(left, right, bam_file, genome_len)
#    minimap2 -x map-pb -a -t 30 --sam-hit-only sample.fa reads.fastq.gz | samtools sort -o output.bam


#    samtools view -b -S sample.sam > sample.bam   # convert .sam to .bam
#    samtools sort sample.bam -o sample_sorted.bam # sort .bam file
#    samtools index sample_sorted.bam              # create index (.bam.bai file)



if __name__ == "__main__":
    if len(sys.argv) != 3:
        print ("Usage: " + sys.argv[0] + " <assembly_info.txt file from metaflye assembly> <output_file>")
#        exit()
    depth = sys.argv[1]
    bam = sys.argv[2]
    [left, right] = extract_suspicious_depth(depth)
    if left != -1:
        check(left, right, bam, 94435)