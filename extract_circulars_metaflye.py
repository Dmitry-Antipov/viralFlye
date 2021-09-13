#!/usr/bin/python   
import sys
import os
import random
#minlength = 5000
maxlength = 1000000
mincov = 10


#        >Utg837180 LN:i:9417 RC:i:34 XO:i:1
def extract_circulars_raven(contigs_file, output_file, minlen_limit, rl):
    fo = open(output_file, "w")
    for line in open(contigs_file, 'r'):
        if line[0] == ">":
            arr = line.split()
            length = int(arr[1].split(':')[2])
            reads = int(arr[2].split(':')[2])
            circulars = int(arr[3].split(':')[2])
            cov = reads * rl / length 
            if (circulars == 1 and length > minlen_limit and length < maxlength and cov > mincov):
                fo.write(arr[0][1:] + "\n")
#L       Ctg882442       +       SRR10963010.1217762     +       12046M

def extract_linears_raven(contigs_file, graph_file, output_file, minlen_limit, rl):
    nonisolated = set()
    fo = open(output_file, "w")

    for line in open(graph_file, 'r'):
        if line[0] == "L":
            arr = line.split()
            nonisolated.add(arr[1])
            nonisolated.add(arr[3])
    for line in open(contigs_file, 'r'):
        if line[0] == ">":
            arr = line.split()
            length = int(arr[1].split(':')[2])
            reads = int(arr[2].split(':')[2])
            circulars = int(arr[3].split(':')[2])
            cov = reads * rl / length 

            if (not(arr[0] in nonisolated) and length > minlen_limit and length < maxlength and cov > mincov):
                fo.write(arr[0][1:] + "\n")

def extract_circulars(info_file, output_file, minlen_limit):
    #fasta = sys.argv[2]
    contigs_list = extract_contigs(info_file, minlen_limit)
    fo = open(output_file, "w")
    for line in open(info_file, 'r'):
        arr = line.split()
        if len(arr) < 4:
            continue
        if arr[0] == "#seq_name":
            continue
        length = int(arr[1])
        cov = float(arr[2])
        circ = arr[3]
        if arr[3] == "Y" and length > minlen_limit and length < maxlength and cov >mincov:
            fo.write(arr[0] + "\n")

def extract_linears(info_file, output_file, minlen_limit):
    #fasta = sys.argv[2]
    contigs_list = []
    fo = open(output_file, "w")
    for line in open(info_file, 'r'):
        arr = line.split()
        if len(arr) < 4:
            continue
        if arr[0] == "#seq_name":
            continue
        length = int(arr[1])
        cov = float(arr[2])
        circ = arr[3]
        edge_seq = arr[7]
        edges = edge_seq.split(",")
        if edges[0] == "*" and edges[-1] == "*":
            if length > minlen_limit and length < maxlength and cov >mincov:
                fo.write(arr[0] + "\n")
#                fo.write(edge_seq + "\n")
#for line in goodf:
#    good.add(line.strip())
#res = unique.intersection(good)
#print res
#print len(res)
if __name__ == "__main__":
    if len(sys.argv) != 3:
        print ("Usage: " + sys.argv[0] + " <assembly_info.txt file from metaflye assembly> <output_file>")
        exit()
    info = sys.argv[1]
    out_f = sys.argv[2]
#    extract_circulars(info, out_f)
    extract_linears(info, out_f)
