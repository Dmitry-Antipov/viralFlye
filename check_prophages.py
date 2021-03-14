#!/usr/bin/python3
import sys
import os
import random
import logging
from Bio.Seq import Seq
inf = 1e8
import gfapy
#TODO: add logging
# https://github.com/ablab/IsoQuant/blob/master/isoquant.py#L248
# https://docs.python.org/3/howto/logging.html

edges = {}
covs = {}
lngth = {}
rel_cov = 1.3

def check_outside_edge(a):
    return lngth[a] > 2000 and covs[a] > 5

def check_prophage_edge(a):
    return lngth[a] < 50000 

def check_all_covs(start, prophage, fin):
    if covs[start] > rel_cov* covs[fin]:
        return False
    if covs[fin] > rel_cov* covs[start]:
        return False
    if covs[prophage] * rel_cov * 2 > covs[fin] + covs[start]:
        return False
    return True


if __name__ == "__main__":
    gfaf = sys.argv[1]
    graph = gfapy.Gfa.from_file(gfaf)
    orientations = {"-":"+","+":"-"}
    for s in graph.segments:
        for o in orientations.keys():
            edges[s.name + o] = set()
            covs[s.name + o] = s.get("dp")
            lngth [s.name + o] = len(s.sequence) 

    for s in graph.edges:
        from_s = [s.from_segment.name + s.from_orient, s.to_segment.name + orientations[s.to_orient]]
        to_s = [s.to_segment.name + s.to_orient, s.from_segment.name + orientations[s.from_orient]]
        for i in range(0, 2):
            edges[from_s[i]].add(to_s[i])
    to_vv = set()
    for start in edges:
        if len(edges[start]) == 2 and check_outside_edge(start):
            for prophage in edges[start]:
                if len(edges[prophage]) == 1:
                    fin = next(iter(edges[prophage]))
                    if check_prophage_edge(prophage) and check_outside_edge(fin) and fin in edges[start]:
                        if check_all_covs(start, prophage, fin):
                            print (f"found {start} {prophage} {fin}")
                            
                            to_vv.add(prophage[:-1])    
                    
#print (total)
#for f in unique:
#    print (f)

