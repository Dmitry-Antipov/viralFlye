#!/usr/bin/python3
import sys
import os
import random
from scipy import signal
import pysam
import numpy as np
from os.path import join
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
    with open(infile, 'r') as input:
        noteof = True
        line1= ""
        try:
            line1 = next(input)
        except:
            noteof = False
            
        while line1 and noteof:
            curname = ""
            if len(line1) > 0 and line1[0] == '>':
                curname = line1.split()[0][1:] + ".fasta"
                out = open(os.path.join(outdir, curname), "w")
                while True:
                    out.write(line1)
                    try:
                        line1 = next(input)
                    except:
                        noteof = False
                        break
                    if (not line1) or len(line1) == 0 or line1[0] == '>':
                        break 
    print("Splitted")


def extract_suspicious_depth(depth_file):
    print (f'Checking contig {os.path.basename(depth_file)} ...')

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
        print('suspicious breakpoints not found')
        return [-1, -1, genome_len]
    elif len(good_right) > 1 or len(good_left) > 1:
        print('multiple suspicious coverage breakpoints found')
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


def count_traversing_fraction (start_pos, end_pos, bamfile, name):
    traversing = 0
    for read in bamfile.fetch(name , start_pos, end_pos):
        if start_pos - epsilon > read.reference_start and end_pos + epsilon < read.reference_end :
            traversing += 1


    count = 0
    for read in bamfile.fetch(name, (start_pos + end_pos) / 2, (start_pos + end_pos + 2) / 2 ):
        count += 1
    print (f' traversing {traversing} of {count}')
    return [traversing, count]

def check (start_pos, end_pos, bam, genome_length, name):
    bamfile = pysam.AlignmentFile(bam, "rb")
    if genome_length < (end_pos - start_pos) * 2:
        print('Suspicious region longer than half of genome')
        return False
    [traversing, count] = count_traversing_fraction(start_pos, end_pos, bamfile, name)
    sample = 0
    random_traversing = 0
    random_count = 0
    print ("checking random samples")
    while sample < sampling_size:
        check_start = random.randint(0, genome_length - end_pos + start_pos)
        if check_start > start_pos - epsilon and check_start < end_pos + epsilon:
            continue
        check_end = check_start + end_pos - start_pos
        ttraverse, tcount = count_traversing_fraction(check_start, check_end, bamfile, name)
        random_traversing += ttraverse
        random_count += tcount
        sample += 1
#TODO what if repeat is longer than read length
    if traversing * random_count <= count * random_traversing * 10:
        print ("significantly less traversing, repeat confirmed")
        return True
    else:
        return False

#TODO multiple reads files, other options for minimap
def prepare_index_and_depth(circulars, assembly, reads, workdir):

    os.makedirs(workdir, exist_ok=True)
    contig_names = []
    for line in open (circulars, 'r'):
        if len(line) > 0 and line[0] == ">":
            contig_names.append(line.strip()[1:])
#    split_multifasta(circulars, workdir)
    res = open(join(workdir, "linears.txt"), 'w')
    bam_file = join(workdir,  "long_reads_realignment.bam")
    bam_line = f'minimap2 -x map-pb -a -t 30 --sam-hit-only --secondary=no {assembly} {reads} | samtools sort -o {bam_file}'
    print(bam_line)
    os.system(bam_line)
    os.system(f'samtools index {bam_file}')
    for contig in contig_names:
        depth_file = join(workdir, f'{contig}.depth')
        bam_contig =  join(workdir, f'{contig}.bam')

        samtools_line = f'samtools view -b {bam_file} {contig} > {bam_contig}; samtools index {bam_contig}'
        depth_line = f'samtools depth -a {bam_contig} > {depth_file}'
#        print(samtools_line)
        os.system(samtools_line)
        os.system(depth_line)
        [left, right, genome_len] = extract_suspicious_depth(depth_file)
        if left != -1:
            try:
                if check(left, right, bam_contig, genome_len, contig):
                    res.write(f'{contig} {left} {right} {right - left} \n')
            except:
                print(f'something wrong in pysam with contig {name}')

'''
    for contig in os.listdir(workdir):
        print (contig)
        if contig.split('.')[-1] == "fasta":
            full_contig = join(workdir, contig)
            name = contig.split('.')[0]
            if not os.path.exists(bam_file):
                print (f'Running minimap for {contig}')
                bam_line = f'minimap2 -x map-pb -a -t 30 --sam-hit-only {full_contig} {reads} | samtools sort -o {bam_file}'
                os.system(bam_line)
                os.system (f'samtools index {bam_file}')            
            depth_file = join(workdir, name +".depth")
            depth_line = f'samtools depth -a {bam_file} >  {depth_file}'
            os.system(depth_line)
            [left, right, genome_len] = extract_suspicious_depth(depth_file)
            if left != -1:
                try:
                    if check(left, right, bam_file, genome_len, name):
                        res.write(f'{name} {left} {right} {right - left} \n')
                except:
                    print(f'something wrong in pysam with contig {name}')
#    minimap2 -x map-pb -a -t 30 --sam-hit-only sample.fa reads.fastq.gz | samtools sort -o output.bam


#    samtools view -b -S sample.sam > sample.bam   # convert .sam to .bam
#    samtools sort sample.bam -o sample_sorted.bam # sort .bam file
#    samtools index sample_sorted.bam              # create index (.bam.bai file)

'''


if __name__ == "__main__":
    if len(sys.argv) != 4:
        print ("Usage: " + sys.argv[0] + " <circular viral contigs> <all contigs> <reads>")
        exit()
    prepare_index_and_depth(sys.argv[1], sys.argv[2], sys.argv[3], os.path.join(os.path.dirname(sys.argv[1]), "linear_check"))
'''
    depth = sys.argv[1]
    bam = sys.argv[2]
    [left, right] = extract_suspicious_depth(depth)
    if left != -1:
        check(left, right, bam, 94435)
'''
