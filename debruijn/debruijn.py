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
    		if (i!=0) and (i!=(len(path)-1)):
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
	remove_paths(g,cutpath,delete_entry_node,delete_sink_node)
	return g

def solve_bubble():
	pass

def simplify_bubbles():
	pass

def solve_entry_tips():
	pass
	
def solve_out_tips():
	pass

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
    g = build_graph(d)
    # visuGraph(g)
    st_nod = get_starting_nodes(g)
    sk_nod = get_sink_nodes(g)
    contigs = get_contigs(g,st_nod,sk_nod)
    # save_contigs(contigs,args.output_file)
    # path_average_weight(g)

if __name__ == '__main__':
    main()
