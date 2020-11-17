#!/usr/bin/python3
import sys
import os
import random
import logging
from Bio.Seq import Seq
neighbours = {}
global_used = set()
segments = {}
links = {}


#for each edge 4 vertices
edges_to_id = {}
min_interesting_size = 1000
max_interesting_size = 500000
min_coverage = 5
class Vertex:
    def __init__ (self, ver_id, edge_id):
#vertices are segments starts & ends, edges are either links(empty seq) or segments (stored in segments ends)
        self.seq = ""
        self.next = []
        self.id = ver_id
        self.edge = edge_id
        self.orientation = "+"
    def __str__(self):
        return (f'edge {self.edge}, id {self.id}, next: {self.next} ')
class node_stat:
    def __init__ (self, length, cov, seq):
        self.length = int(length)
        self.cov = float (cov)
        self.seq = seq

def rc (seq):

    seqS = Seq(seq)

    return seqS.reverse_complement().__str__()
'''
    res = seq
    l = len(seq)
    table = {}
    table['A']= 'G'
    table['G'] = 'A'
    table['C'] = 'T'
    table['T'] = 'C'
    table['N'] = 'N'
    for i in range (0, l):
        if not seq[i] in table:
            print (f'{seq} {i} strange symbol {seq[i]}')
            exit(0)
        res[l - 1 - i] = table[seq[i]]
    return res
'''

def get_ids(link_name):
    arr = link_name.split()
    res = [arr[1], arr[3]]    
    return res

def get_small_component(id, min_size, max_size, min_cov):
    global neighbours
    global global_used
    global segments
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
        total_len += segments[f].length
        total_cov += segments[f].length * segments[f].cov
    total_cov /= total_len
#    if total_len >= low_cutoff and total_len <= high_cutoff and total_cov > min_coverage: 
#    if total_len >= ids[id].length and total_len > 1000 and total_cov > 10 and len(used) < 2:
    if total_len < max_size and total_len > min_size and total_cov > min_cov and total_len > segments[id].length:
        return used
    else:
        return set()
#params: fastg file, ids file, max component size
#looks for ids that contains in small connected components of size 1< SIZE <= max_cutoff 
def construct_graph(edge_component):
    global segments
    global links
    global edges_to_id
    vertices = {}
    id_count = 0
    edge_count = 0
    edges_to_id = {}
    for e in edge_component:
        #starts, end, rc_start, rc_end   0,1,2,3, rc_id = 3 - id
        for i in range(0, 4):
            vertices[id_count] = Vertex(id_count, e)
            id_count += 1
        edges_to_id[e] = edge_count
        vertices[edge_count * 4].next.append(edge_count * 4 + 1)
        vertices[edge_count * 4 + 2].next.append(edge_count * 4 + 3)
        vertices[edge_count * 4 + 1].seq = segments[e].seq
        vertices[edge_count * 4 + 3].seq = rc(segments[e].seq)
        vertices[edge_count * 4 + 3].orientation = "-"

        edge_count += 1

#adding edges from lines
    for e in edge_component:
        for l in links[e]:
            arr = l.split()
    #L	edge_43116	-	edge_6653	+	0M
            edge_start = e
            edge_end = arr[3]
            link_start_shift = 1
            if arr[2] == '-':
                link_start_shift = 3
            link_end_shift = 0
            if arr[4] == '-':
                link_end_shift = 2


            start_id = edges_to_id[edge_start] * 4 + link_start_shift
            end_id = edges_to_id[edge_end]* 4 + link_end_shift
            print (f'from link {l} adding edge {start_id} {end_id} ')
            vertices[start_id].next.append(end_id)

#rc_Link
            start_id = edges_to_id[edge_start] * 4 + 3 - link_start_shift
            end_id = edges_to_id[edge_end]* 4 + 3 - link_end_shift
            print(f'from link {l} adding edge {end_id} {start_id} ')
            vertices[end_id].next.append(start_id)
    return vertices

def get_start_end_vertex(edge_component, vertices):
    global segments
    global links
    global edges_to_id
    max_l = 0
    max_e = ""
    for e in edge_component:
        if segments[e].length > max_l:
            max_l = segments[e].length
            max_e = e
    return [edges_to_id[max_e] * 4 + 1, edges_to_id[max_e] * 4]

def run_dikstra(vertices, source, target):
    q = set()
    inf = 1e8
    dist = {}
    prev = {}
    for v in vertices:
        dist[v] = inf
        prev[v] = -1
        q.add(v)
    dist[source] = 0
    while len(q) > 0:
        min_d = inf
        min_v = -1
        for v in q:
            if dist[v] < min_d:
                min_d = dist[v]
                min_v = v
        u = min_v
        if min_d == inf:
            break
        q.remove(min_v)
        for next_u in vertices[u].next:
            alt = dist[u] + len(vertices[next_u].seq)
            if alt < dist[next_u]:
                dist[next_u] = alt
                prev[next_u] = u
    return dist, prev

if len (sys.argv) != 6:
    print(f'Script for a random traversal of suspicious components(for their further examination on virality)')
    print(f'Usage: {sys.argv[0]} <min_component_length>  <max_component_length> <min_coverage> <output file>')

low_cutoff = int (sys.argv[2])
high_cutoff = int (sys.argv[3])
min_coverage = float(sys.argv[4])
resf = open(sys.argv[5], "w")
good = set()

for line in open(sys.argv[1], 'r'):
    if line[0] == "L":
        arr = get_ids(line)
        if len(arr) == 0:
            print (line)
            exit()
        neighbours[arr[0]].add(arr[1])
        neighbours[arr[1]].add(arr[0])
        if not arr[0] in links:
            links[arr[0]] = []
        if not arr[1] in links:
            links[arr[1]] = []
        links[arr[0]].append(line)

    elif line[0] == "S":
        arr = line.split()
        length = len(arr[2])
        cov = arr[3].split(':')[2]
        segments[arr[1]] = node_stat(length, cov, arr[2])
        neighbours[arr[1]] = set()
unique = set()
total = 0
for seq_id in segments:
    s =  get_small_component(seq_id, low_cutoff, high_cutoff, min_coverage)
    if len (s) > 0:
        print (f'exporing small component of edges {s}')
        vertices = construct_graph(s)

        [source, target] = get_start_end_vertex(s, vertices)
        dist, prev = run_dikstra (vertices,source, target)
#        if "edge_1" in s:
#            for v in vertices:
#                print (vertices[v])
#            for e in s:
#                print (links[e])
#            print (f' source {source} target {target}')
#            print (f' {dist} {prev}')

        path = []
        if prev[target] != -1:
            while target != source:
                path.append(target)
                target = prev[target]
            path.append(source)
            path.reverse()
#            exit()
            sum_len = 0
            for e in s:
                sum_len += segments[e].length
            res =""
            for v in path:
                res += vertices[v].seq
            path_len = len(res)
            header = f'edges_{len(s)}_length_{path_len}_total_{sum_len}_'
            for v in path:
                print (f'{v} {vertices[v].edge}')
                if v%2 == 1:
                    header += vertices[v].edge
                    header += vertices[v].orientation
            if sum_len * 0.7 > path_len:
                print(f'Path is small ({path_len} of {sum_len}). Complex component of {len(s)} edges. Header: {header}')
            else:
                resf.write(">" + header + "\n")
                resf.write(res + "\n")
        else:
            print ("no path from source to sink")
#print (total)
#for f in unique:
#    print (f)
