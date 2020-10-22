#!/usr/bin/python   
import sys
import os
import random
minlength = 1000
maxlength = 1000000
mincov = 5

def extract_circulars(info_file):
    #fasta = sys.argv[2]
    contigs_list = []
    for line in open(info, 'r'):
        arr = line.split()
        if len(arr) < 4:
            continue
        if arr[0] == "#seq_name":
            continue
        length = int(arr[1])
        cov = float(arr[2])
        circ = arr[3]
        if arr[3] == "Y" and length > minlength and length < maxlength and cov >mincov:
            print (arr[0])
#for line in goodf:
#    good.add(line.strip())
#res = unique.intersection(good)
#print res
#print len(res)

if len(sys.argv) != 2:
    print ("Usage: " + sys.argv[0] + " <assembly_info.txt file from metaflye assembly>")
    exit()
info = sys.argv[1]
extract_circulars(info)
