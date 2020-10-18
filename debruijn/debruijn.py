#!/bin/env python3
# -*- coding: utf-8 -*-
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#    A copy of the GNU General Public License is available at
#    http://www.gnu.org/licenses/gpl-3.0.html

"""Perform assembly based on debruijn graph."""

import argparse
import os
import sys
import networkx as nx
from operator import itemgetter
import random
random.seed(9001)
from random import randint
import statistics
import matplotlib.pyplot as plt
import copy

__author__ = "Your Name"
__copyright__ = "Universite Paris Diderot"
__credits__ = ["Your Name"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Your Name"
__email__ = "your@email.fr"
__status__ = "Developpement"

def isfile(path):
    """Check if path is an existing file.
      :Parameters:
          path: Path to the file
    """
    if not os.path.isfile(path):
        if os.path.isdir(path):
            msg = "{0} is a directory".format(path)
        else:
            msg = "{0} does not exist.".format(path)
        raise argparse.ArgumentTypeError(msg)
    return path


def get_arguments():
    """Retrieves the arguments of the program.
      Returns: An object that contains the arguments
    """
    # Parsing arguments
    parser = argparse.ArgumentParser(description=__doc__, usage=
                                     "{0} -h"
                                     .format(sys.argv[0]))
    parser.add_argument('-i', dest='fastq_file', type=isfile,
                        required=True, help="Fastq file")
    parser.add_argument('-k', dest='kmer_size', type=int,
                        default=21, help="K-mer size (default 21)")
    parser.add_argument('-o', dest='output_file', type=str,
                        default=os.curdir + os.sep + "contigs.fasta",
                        help="Output contigs in fasta file")
    return parser.parse_args()


#==============================================================
# Main program
#==============================================================
def read_fastq(fastq):
    with open(fastq, 'r') as f:
        for l in f:
            yield(next(f).strip())
            next(f)
            next(f)

def cut_kmer(seq, ksize):
	for i in range(len(seq)):
		if (i + ksize) <= len(seq):
			yield(seq[i:i+ksize])

def build_kmer_dict(fastq,ksize):
    occu = {}
    for seq in read_fastq(fastq):
        for kmer in cut_kmer(seq,ksize):
            if kmer not in occu:
                occu[kmer] = 0    
            occu[kmer] += 1
    return(occu)

def build_graph(d_Kmer):
	G = nx.DiGraph()
	for k in list(d_Kmer.keys()):
		n = len(k)
		k1 = k[:n-1]
		k2 = k[1:n]
		G.add_node(k1)
		G.add_node(k2)
		G.add_edge(k1,k2)
		G[k1][k2]['weight'] = d_Kmer[k]
	return(G)

def get_starting_nodes(g):
    start_nodes = []
    for node in g.nodes():
        if not bool(g.pred[node]):
            start_nodes.append(node)
    return start_nodes

def get_sink_nodes(g):
    end_nodes = []
    for node in g.nodes():
        if not bool(g.succ[node]):
            end_nodes.append(node)
    return end_nodes

def get_contigs(g,snodes,fnodes):
    contigs = []
    for sn in snodes:
        for fn in fnodes:
            for path in nx.all_simple_paths(g,sn,fn):
                contig = []
                contig.append(path[0][0:-1]) #départ du path jusqu'a dernière lettre
                for neud in path:
                    contig.append(neud[-1])
                contig = ''.join(contig)
                taille = len(contig)
                contigs.append((contig,taille))
    return contigs

def save_contigs(contigs,filename):
	with open(filename, 'w') as o:
		for i,ele in enumerate(contigs):
			o.write(">contig_{} len={}\n{}\n".format(i,ele[1],fill(ele[0])))

def fill(text, width=80):
    """Split text with a line return to respect fasta format"""
    return os.linesep.join(text[i:i+width] for i in range(0, len(text), width))

def std(liste):
    return(statistics.stdev(liste))

def path_average_weight(g,path):
    w = 0
    for i in range(len(path)-1):
        w += g[path[i]][path[i+1]]["weight"]
    return(w/(i+1))

def remove_paths(g,lpath,delete_entry_node,delete_sink_node):
    for path in lpath:
    	for i,neud in enumerate(path):
    		if (i!=0) and (i!=(len(path)-1)) and neud in g.nodes():
    			g.remove_node(neud)
    	if delete_entry_node:
    		g.remove_node(path[0])
    	if delete_sink_node:
    		g.remove_node(path[len(path)-1])
    return g		

def select_best_path(g,lpath,llenpath,lmeanwei,delete_entry_node=False,delete_sink_node=False):
    cutpath = []
    for pi in range(len(lpath)):
        for pj in range((pi+1),len(lpath)):
            stdl = std([llenpath[pi],llenpath[pj]])
            stdw = std([lmeanwei[pi],lmeanwei[pj]])
            if stdl == 0 and stdw == 0:
                rd = random.choice([pi,pj])
                cutpath.append(lpath[rd])
            elif stdl != 0 and stdw == 0:
                if llenpath[pi] > llenpath[pj]:
                    cutpath.append(lpath[pj])
                else :
                    cutpath.append(lpath[pi])	
            else:
                if lmeanwei[pi] > lmeanwei[pj]:
                    cutpath.append(lpath[pj])
                else :
                    cutpath.append(lpath[pi])	
    print("CUTPATH : {}\n".format(cutpath))	
    g = remove_paths(g,cutpath,delete_entry_node,delete_sink_node)
    return g

def solve_bubble(g,nanc,ndes):
    lpath = []
    llenpath = []
    lmeanwei = []
    for path in nx.all_simple_paths(g,nanc,ndes):
        lpath.append(path)
        llenpath.append(len(path))
        lmeanwei.append(path_average_weight(g,path))
    g = select_best_path(g,lpath,llenpath,lmeanwei)
    return g
"""
def simplify_bubble(g):
    listbubble = []
    in_ = get_starting_nodes(g)
    out_ = get_sink_nodes(g)
    to_del = []
    for ni in in_:
        for nj in out_:
            print("NI {}".format(ni))
            lpath = []
            inNod = set([])
            for path in nx.all_simple_paths(g,ni,nj):
                lpath.append(path)
            for i in range(1,len(lpath)):
                pbuble = set(lpath[i-1])
                diff = pbuble.symmetric_difference(set(lpath[i]))
                inNod = inNod.union(diff)
            print("inod{}".format(inNod))
            rmv = []
            stop = False
            for nbu in list(inNod):
                print(nbu)
                for neud in rmv:
                    if nbu == neud:
                        stop = True
                if stop :
                    stop = False
                    continue
                npred = list(g.predecessors(nbu))[0]
                rmv.append(npred)
                rmv.append(nbu)
                nbsucc = g.successors(npred)
                # print("nbsucc : {} npred : {}".format(nbsucc,npred))
                while len(list(nbsucc))<=1:
                    npred = list(g.predecessors(npred))[0]
                    nbsucc = g.successors(npred)
                    rmv.append(npred)
                    # print("nbsucc : {} npred : {}".format(nbsucc,npred))
                nsuiv = list(g.successors(nbu))[0]
                rmv.append(nsuiv)
                nbpred = g.predecessors(nsuiv)
                while len(list(nbpred))<=1:
                    nsuiv = list(g.successors(nsuiv))[0]
                    nbpred = g.predecessors(nsuiv)
                    rmv.append(nsuiv)
                print("nsuiv : {}\nnpred : {}\nrmv : {}".format(nsuiv,npred,rmv))
                to_del.append([npred,nsuiv])
    print(to_del)
    cpydel = copy.deepcopy(to_del)
    for pi in range(len(to_del)):
        for path in nx.all_simple_paths(g,to_del[pi][0],to_del[pi][1]):
            for pj in range(len(to_del)):
                if pj == pi:
                    continue 
                print("PATH {}".format(path[1:-1]))
                if to_del[pj][0] in path[1:-1] and to_del[pj][1] in path[1:-1] and to_del[pj] in cpydel:
                    print("RMV {}".format(to_del[pj])
                    cpydel.remove(to_del[pj])
    print(cpydel)
    g = solve_bubble(g,npred,nsuiv)
    return g
"""
def test_simplify_bubbles():
    graph_1 = nx.DiGraph()
    graph_1.add_weighted_edges_from([(1, 2, 10), (2, 4, 15), (2, 3, 15),
                                     (2, 5,10), (3, 4,10), (3, 6, 3),
                                     (4, 6, 3), (5, 7, 3), (6, 7, 10),
                                     (7, 8, 10),(7, 9, 10),(8, 10, 1),
                                     (9, 12, 10),(12, 10, 3),(10, 11, 4),
                                     (36, 5, 4)])
    graph_1 = simplify_bubbles(graph_1)

def simplify_bubbles(g):
    start = get_starting_nodes(g)
    sink = get_sink_nodes(g)
    for pi in start:
        for pj in sink:
            lpath = list(nx.all_simple_paths(g,pi,pj))
            while len(lpath) >= 2:
                print(lpath[0],lpath[1])
                i=0
                print("ARBI")
                while lpath[0][i] in lpath[1] and i < len(lpath[0])-1:
                    print(lpath[0][i])
                    i = i + 1
                if i == len(lpath[0]):
                    l_1 = lpath[1]
                    l_0 = lpath[0]
                    while l_1[i] in l_0:
                        i = i + 1
                    bulle = l_1[i-1]
                else:
                    bulle = lpath[0][i-1]
                i=0
                print("BABR")
                while not lpath[0][i+lpath[0].index(bulle)+1] in lpath[1]:
                    print(lpath[0][i+lpath[0].index(bulle)+1])
                    i += 1
                    print(lpath[0][i+lpath[0].index(bulle)+1])
                bulls = lpath[0][i+lpath[0].index(bulle)+1]
                print("bull : {} bulls : {}".format(bulle,bulls)) 
                g = solve_bubble(g,bulle,bulls)
                lpath = list(nx.all_simple_paths(g,pi,pj))
                
    return g

def solve_entry_tips(g,lentry):
    letips = []    
    for enod in lentry:
        npred = 1
        snod = enod
        while npred <= 1:
            pred = snod
            print(pred)
            snod = list(g.successors(pred))[0]
            if not bool(list(g.successors(snod))):
                break
            npred = len(list(g.predecessors(snod)))
            print(npred)
        letips.append([enod,snod])
    lmeanwei = []
    llenpath = []
    lpath = []
    for tips in letips:
        for path in nx.all_simple_paths(g,tips[0],tips[1]):
            lmeanwei.append(path_average_weight(g,path))
            llenpath.append(len(path))
            lpath.append(path)
    g = select_best_path(g,lpath,llenpath,lmeanwei,True)
    return g

def solve_out_tips(g,lout):
    lotips = []
    for onod in lout:
        nsuc = 1
        pnod = onod
        while nsuc <= 1:
            suiv = pnod
            pnod = list(g.predecessors(suiv))[0]
            nsuc = len(list(g.successors(pnod)))
            if not bool(list(g.predecessors(pnod))):
                break 
        lotips.append([pnod,onod])
    lmeanwei = []
    llenpath = []
    lpath = []
    for tips in lotips:
        for path in nx.all_simple_paths(g,tips[0],tips[1]):
            lmeanwei.append(path_average_weight(g,path))
            llenpath.append(len(path))
            lpath.append(path)  
    g = select_best_path(g,lpath,llenpath,lmeanwei,False,True)
    return g

def visuGraph(g):
    nx.draw(g,with_labels = True, font_weight='bold')
    plt.show()

def main():
    """
    Main program function
    """
    # Get arguments
    args = get_arguments()
    d = build_kmer_dict(args.fastq_file, args.kmer_size)
    # print("STEP1")
    g = build_graph(d)
    # print("STEP2")
    # visuGraph(g)
    st_nod = get_starting_nodes(g)
    sk_nod = get_sink_nodes(g)
    print("STEP3")
    g = simplify_bubbles(g)    
    # test_simplify_bubbles()
    print("STEP4")
    g = solve_entry_tips(g,st_nod)
    print("STEP5")
    g = solve_out_tips(g,sk_nod)
    print("STEP6")
    contigs = get_contigs(g,get_starting_nodes(g),get_sink_nodes(g))
    print("STEP7")
    save_contigs(contigs,args.output_file)
    print("STEP8")

if __name__ == '__main__':
    main()
