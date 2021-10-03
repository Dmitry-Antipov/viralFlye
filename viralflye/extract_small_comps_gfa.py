#!/usr/bin/python3
import sys
import os
import random
import logging
from Bio.Seq import Seq
inf = 1e8
#TODO: add logging
# https://github.com/ablab/IsoQuant/blob/master/isoquant.py#L248
# https://docs.python.org/3/howto/logging.html


#for each edge 4 vertices

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
#starts, end, rc_start, rc_end   0,1,2,3, rc_id = 3 - id

    def get_rc_id(self):
        return (self.id // 4) * 4 + 3 - (self.id % 4)
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


def get_small_component(id, min_size, max_size, min_cov, neighbours, global_used, segments):
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


def get_header_id(header):
    pos = header.find("edge_")
    return header[:pos]

def construct_graph(edge_component, segments, links):
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
#            print (f'from link {l} adding edge {start_id} {end_id} ')
            vertices[start_id].next.append(end_id)

#rc_Link
            start_id = edges_to_id[edge_start] * 4 + 3 - link_start_shift
            end_id = edges_to_id[edge_end]* 4 + 3 - link_end_shift
#            print(f'from link {l} adding edge {end_id} {start_id} ')
            vertices[end_id].next.append(start_id)
    return [vertices, edges_to_id]

def get_start_end_vertex(edge_component, segments, edges_to_id):
    max_l = 0
    max_e = ""
    for e in edge_component:
#        print (f'edge {e} length {segments[e].length}')
        if segments[e].length > max_l:
            max_l = segments[e].length
            max_e = e
#    print (max_e)
    return [edges_to_id[max_e] * 4 + 1, edges_to_id[max_e] * 4]

def run_dikstra(vertices, source):
    global inf
    q = set()
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

def get_furtherst_id(dist):
    global inf
    max_d = -1
    max_v = -1
    for v in dist:
        if dist[v] > max_d and dist[v] < inf:
            max_d = dist[v]
            max_v = v
    return max_v

def restore_path(prev, source, target):
    path = []
#    print (prev)
#    print (f'{source} {target} {prev[target]}')

    while target != source:
        if prev[target] == -1:
            print('error from dijkstra')
            exit(0)
        path.append(target)
        target = prev[target]
    path.reverse()
    return path

def extract_paths_in_components(input_file, low_cutoff, high_cutoff, min_coverage, output_file):
    neighbours = {}
    global_used = set()
    segments = {}
    links = {}

    resf = open(output_file, "w")
    stats_f = open(os.path.splitext(output_file)[0] + ".stats", "w" )
    pred = {}
    '''
    pred_file = os.path.dirname(output_file) + "/vv_comp_rerun/components_rerun_result_table.csv" 
    if not os.path.exists(pred_file):
        return

    prediction_csv =  open(pred_file, "r")
    for line in prediction_csv:
        arr = line.split(',')
        pred[get_header_id(arr[0])] = arr[1] '''

    good = set()
    for line in open(input_file, 'r'):
        if line[0] == "L":
            arr = get_ids(line)
            if len(arr) == 0:
#                print (line)
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
        s =  get_small_component(seq_id, low_cutoff, high_cutoff, min_coverage, neighbours, global_used, segments)
        if len (s) > 0:
#            print (f'exporing small component of edges {s}')
            vertices, edges_to_id = construct_graph(s, segments, links)

            [source, target] = get_start_end_vertex(s, segments, edges_to_id)
            dist, prev = run_dikstra (vertices, source)
            if prev[target] != -1:
                path = restore_path(prev, source, target)
                path.append(source)
    # longest path containing longest edge with no cycles:
            else:
                new_end = get_furtherst_id(dist)
#                print (f'new_end {new_end}')
                forward_path = restore_path(prev, source, new_end)
                target_rc = vertices[target].get_rc_id()
#                print (f'target target_rc {target} {target_rc}')
                [new_dist, new_prev] = run_dikstra(vertices, target_rc)
                print (new_dist)
                new_start_rc = get_furtherst_id(new_dist)
#                print(f'new_start_rc {new_start_rc}')
#                print (f'forward_path {forward_path}')
                backward_path = restore_path(new_prev, target_rc, new_start_rc)
#                print(f'backward_path {backward_path}')
                backward_path_rc = []
                for p in backward_path:
                    backward_path_rc.append(vertices[p].get_rc_id())
                backward_path_rc.reverse()
                backward_path_rc.append(source)
                backward_path_rc.extend(forward_path)
                path = backward_path_rc
    #            exit()
            sum_len = 0
            total_cov = 0
            sum_cov = 0
            for e in s:
                sum_len += segments[e].length
                sum_cov += segments[e].cov
                total_cov += segments[e].cov * segments[e].length
            res = ""
            for v in path:
#                print (v)
                res += vertices[v].seq
            path_len = len(res)
#>Utg61358 LN:i:45259 RC:i:17 XO:i:1
            header = f'edges_{len(s)}_length_{path_len}_total_{sum_len}_'
            for v in path:
#                print (f'{v} {vertices[v].edge}')
                if v%2 == 1:
                    header += vertices[v].edge
                    header += vertices[v].orientation
            header += f' LN:i:{sum_len} RC:i:{sum_cov:.0f}'
            if sum_len * 0.3 > path_len:
                print(f'Path is small ({path_len} of {sum_len}). Complex component of {len(s)} edges. Header: {header}')
            else:
                resf.write(">" + header + "\n")
                resf.write(res + "\n")
                hdr = get_header_id(header)
               
                if hdr in pred and pred[hdr] == "Virus":
                    stats_f.write(f'{header}\t{len(s)}\t{sum_len}\t{(total_cov/sum_len):.2f}\n')

if __name__ == "__main__":
    if len (sys.argv) != 6:
        print(f'Script for a random traversal of suspicious components(for their further examination on virality)')
        print(f'Usage: {sys.argv[0]} <input_file> <min_component_length>  <max_component_length> <min_coverage> <output file>')
        exit(0)
    low_cutoff = int (sys.argv[2])
    high_cutoff = int (sys.argv[3])
    min_coverage = float(sys.argv[4])
    extract_paths_in_components(sys.argv[1], low_cutoff, high_cutoff, min_coverage, sys.argv[5])
    
#print (total)
#for f in unique:
#    print (f)
