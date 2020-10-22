#!/usr/bin/python   
import sys
import os
import random
neighbours = {}
global_used = set()
ids = {}

class node_stat:
    def __init__ (self, length, cov):
        self.length = int(length)
        self.cov = float (cov)
def check (name, cutoff):
    arr = name.split(':')
    if len(arr) == 1:
        return False
    cont_name = arr[0].strip()
    contig = cont_name.split('_')
    if int(contig[3]) < cutoff:
        return False
    return True
#params: fastg 
#
#    print (cont_name[-2])
#    if cont_name[-2] == "\'":
#        return False
def get_ids(link_name):
    arr = link_name.split()
    res = [arr[1], arr[3]]    
    return res

def get_small_component(id, low_cutoff, high_cutoff, min_coverage):
    global neighbours
    global global_used
    global ids
    used = set()
    if id in global_used:
        return used
    in_job = [id]
    while in_job:        
        cur_v = in_job.pop()
        if cur_v not in used:
            used.add(cur_v)
            for v in neighbours[cur_v]:
                if v not in used:    
                    in_job.append(v)
#    print (id + " " + str(len(used)))
    global_used.update(used)
    total_len = 0
    total_cov = 0.0
    for f in used:
        total_len += ids[f].length
        total_cov += ids[f].length * ids[f].cov
    total_cov /= total_len
#    if total_len >= low_cutoff and total_len <= high_cutoff and total_cov > min_coverage: 
#    if total_len >= ids[id].length and total_len > 1000 and total_cov > 10 and len(used) < 2:
    if total_len == ids[id].length and total_len > 1000 and total_cov > 5:
# and len(neighbours[id]) > 0:
#    if len (used) > 2 and len (used) < 10:
        return used
    else:
        return set()
#params: fastg file, ids file, max component size
#looks for ids that contains in small connected components of size 1< SIZE <= max_cutoff 

good = set()
low_cutoff = int (sys.argv[2])
high_cutoff = int (sys.argv[3])
min_coverage = float(sys.argv[4])
for line in open(sys.argv[1], 'r'):
    if line[0] == "L":
        arr = get_ids(line)
        if len(arr) == 0:
            print (line)
            exit()
        neighbours[arr[0]].add(arr[1])
        neighbours[arr[1]].add(arr[0])

    elif line[0] == "S":
        arr = line.split()
        length = len(arr[2])
        cov = arr[3].split(':')[2]
        ids[arr[1]] = node_stat(length, cov)
        neighbours[arr[1]] = set()
unique = set()
total = 0
for id in ids:
    s =  get_small_component(id, low_cutoff, high_cutoff, min_coverage)
    if len (s) > 0:
        total += 1
        unique = unique.union(s)
        print (s)
print (total)
#for f in unique:
#    print (f)
good = set()
#goodf = open ("/Iceking/dantipov/metaFlye/Aloha/SRR10378147_haplo_3/checkv_good.names", 'r')
#for line in goodf:
#    good.add(line.strip())
#res = unique.intersection(good)
#print res
#print len(res)
