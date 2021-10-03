#!/usr/bin/python3
############################################################################
# Copyright (c) 2015 Saint Petersburg State University
# Copyright (c) 2011-2015 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

from __future__ import with_statement
import os
import shlex
import shutil
import stat
import subprocess
import sys
import platform
import re
import gzip
import time
import random
from collections import defaultdict
from genericpath import isdir, exists

from os.path import join

#Do not use sra tools from chihua! it's outdated and bugged
sra_tools = "/home/dantipov/other_tools/sratoolkit.2.10.7-ubuntu64/bin/"
seqtk = "/home/dantipov/other_tools/seqtk/seqtk" 

#edges_4_length_18794_total_18794_edge_17+edge_18+edge_7432+edge_19+
def main(components, graph, outdir):
    components = sys.argv[1]
    
    outdir = sys.argv[3]
    
    os.makedirs(outdir, exist_ok="true")
    contigs = join(os.path.dirname(graph), "edges.fasta")
    contigs_f = open(contigs, "w")
    for line in  open (graph, "r"):
        arr = line.split()
        if len(arr) > 0 and arr[0] == "S":
            contigs_f.write (f'>{arr[1]}\n')
            contigs_f.write (f'{arr[2]}\n')
    component = 0
    for line in open (components, "r"):
        if len(line) > 0 and line[0] == ">":
            arr = line.replace('+','_').replace('-', '_').split('_')
            edges = []
            ind = 7
            while ind < len(arr):
                edges.append(f'edge_{arr[ind]}\n')
                ind +=2
            outf_name = join(outdir, f'component_{component}.list')
            outf = open(outf_name, "w")
            for f in  edges:
                outf.write(f)
            outf.close()
            component += 1
            outf_fasta = join(outdir, f'component_{component}.fasta')
            seqtk_line = (f'{seqtk} subseq {contigs} {outf_name} > {outf_fasta}')
            os.system (seqtk_line)

def count_complete_component(file, frac):
    count = 0
    references = []
    ref_names =[]
    for line in open(file, "r"):
        arr = line.split()
        if arr[0] != "Assemblies":
            ind = 0
            for el in arr[1:]:
                ind += 1
                flag = False
                if el != "-" and float(el) > frac * 100.0:
                    if not flag:
                        count +=1
                    flag = True
                    references[ind].append(arr[0])
        else:
            for el in arr[1:]:
                references.append([])
                ref_names.append(el)
                
    print (f'{count} references assembled with genome fraction > {frac} in our components')
    for i in range (0, len(references)):
        if len(references[i]) > 1:
            print (ref_names[i])
            print (references[i])
if __name__ == '__main__':
#components, graph, outdir
#    main(sys.argv[1],sys.argv[2], sys.argv[3])
    count_complete_component(sys.argv[1], 0.8)


