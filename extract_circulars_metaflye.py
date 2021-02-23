#!/usr/bin/python   
import sys
import os
import random
minlength = 5000
maxlength = 1000000
mincov = 10

def extract_circulars(info_file, output_file):
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
        if arr[3] == "Y" and length > minlength and length < maxlength and cov >mincov:
            fo.write(arr[0] + "\n")

def extract_linears(info_file, output_file):
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
            if length > minlength and length < maxlength and cov >mincov:
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
