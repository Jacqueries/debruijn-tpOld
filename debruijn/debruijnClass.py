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

import random
import statistics
import argparse
import os
import sys
import matplotlib.pyplot as plt
import networkx as nx
random.seed(9001)

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
    """Read fasta file file.
      :Parameters:
          fastq: Path to the file
    """
    with open(fastq, 'r') as file:
        for _ in file:
            yield next(file).strip()
            next(file)
            next(file)

def cut_kmer(seq, ksize):
    """Cut read into kmer.
      :Parameters:
          ksize: size of kmer
          seq: sequence
    """
    for i in range(len(seq)):
        if (i + ksize) <= len(seq):
            yield seq[i:i+ksize]

def build_kmer_dict(fastq,ksize):
    """Cut read into kmer.
      :Parameters:
          ksize: size of kmer
          seq: sequence
    """
    occu = {}
    for seq in read_fastq(fastq):
        for kmer in cut_kmer(seq,ksize):
            if kmer not in occu:
                occu[kmer] = 0
            occu[kmer] += 1
    return occu

def build_graph(dkmer):
    """La fonction ​ build_graph /4 prendra en entrée un dictionnaire de k-mer et créera l’arbre de
k-mers préfixes et suffixes décrit précédemment. Les arcs auront pour paramètre obligatoire un
poids nommé “​ weight​ ”.
    """
    graph = nx.DiGraph()
    for k in list(dkmer.keys()):
        taille = len(k)
        kpre = k[:taille-1]
        ksuf = k[1:taille]
        graph.add_node(kpre)
        graph.add_node(ksuf)
        graph.add_edge(kpre,ksuf)
        graph[kpre][ksuf]['weight'] = dkmer[k]
    return graph

def get_starting_nodes(graph):
    """get_starting_nodes ​ /1 prend en entrée un graphe et retourne une liste de noeuds
d’entrée
    """
    start_nodes = []
    for node in graph.nodes():
        if not bool(graph.pred[node]):
            start_nodes.append(node)
    return start_nodes

def get_sink_nodes(graph):
    """get_sink_nodes​ /1 prend en entrée un graphe et retourne une liste de noeuds de sortie
    """
    end_nodes = []
    for node in graph.nodes():
        if not bool(graph.succ[node]):
            end_nodes.append(node)
    return end_nodes

def get_contigs(graph,snodes,fnodes):
    """get_contigs / 3 prend un graphe, une liste de noeuds d’entrée et une liste de sortie et
retourne une liste de tuple(contig, taille du contig)
    """
    contigs = []
    for snod in snodes:
        for fnod in fnodes:
            for path in nx.all_simple_paths(graph,snod,fnod):
                contig = []
                contig.append(path[0][0:-1]) #départ du path jusqu'a dernière lettre
                for neud in path:
                    contig.append(neud[-1])
                contig = ''.join(contig)
                taille = len(contig)
                contigs.append((contig,taille))
    return contigs

def save_contigs(contigs,filename):
    """save_contigs ​ /2 ​ qui prend une liste de tuple (contig, taille du contig) et un nom de
fichier de sortie et écrit un fichier de sortie contenant les contigs selon le format fasta
(retour chariot tous les 80 charactères) à l’aide de la fonction fill
    """
    with open(filename, 'w') as output:
        for i,ele in enumerate(contigs):
            output.write(">contig_{} len={}\n{}\n".format(i,ele[1],fill(ele[0])))

def fill(text, width=80):
    """Split text with a line return to respect fasta format"""
    return os.linesep.join(text[i:i+width] for i in range(0, len(text), width))

def std(liste):
    """std ​ qui prend une liste de valeur, qui retourne l’écart type.
    """
    return statistics.stdev(liste)

def path_average_weight(graph,path):
    """path_average_weight /1 qui prend un graphe et un chemin et qui retourne un poids
moyen
    """
    weight = 0
    for i in range(len(path)-1):
        weight += graph[path[i]][path[i+1]]["weight"]
    return weight/(i+1)

def remove_paths(graph,lpath,delete_entry_node,delete_sink_node):
    """remove_paths /3 qui prend un graphe et une liste de chemin, la variable booléenne
delete_entry_node pour indiquer si les noeuds d’entrée seront supprimés et la variable
booléenne ​ delete_sink_node pour indiquer si les noeuds de sortie seront supprimés et
retourne un graphe nettoyé des chemins indésirables.
    """
    for path in lpath:
        for i,neud in enumerate(path):
            if (i!=0) and (i!=(len(path)-1)) and neud in graph.nodes():
                graph.remove_node(neud)
        if delete_entry_node:
            graph.remove_node(path[0])
        if delete_sink_node:
            graph.remove_node(path[len(path)-1])
    return graph

def select_best_path(graph,lpath,llenpath,lmeanwei,delete_entry_node=False,delete_sink_node=False):
    """select_best_path /4 qui prend un graphe, une liste de chemin, une liste donnant la
longueur de chaque chemin, une liste donnant le poids moyen de chaque chemin,
delete_entry_node pour indiquer si les noeuds d’entrée seront supprimés et
delete_sink_node pour indiquer si les noeuds de sortie seront supprimés et retourne un
graphe nettoyé des chemins indésirables. Par défaut, ​ delete_entry_node et
delete_sink_node​ seront ici à ​ False​ .
    """
    cutpath = []
    for pathi in range(len(lpath)):
        for pathj in range((pathi+1),len(lpath)):
            stdl = std([llenpath[pathi],llenpath[pathj]])
            stdw = std([lmeanwei[pathi],lmeanwei[pathj]])
            if stdl == 0 and stdw == 0:
                rdm = random.choice([pathi,pathj])
                cutpath.append(lpath[rdm])
            elif stdl != 0 and stdw == 0:
                if llenpath[pathi] > llenpath[pathj]:
                    cutpath.append(lpath[pathj])
                else :
                    cutpath.append(lpath[pathi])
            else:
                if lmeanwei[pathi] > lmeanwei[pathj]:
                    cutpath.append(lpath[pathj])
                else :
                    cutpath.append(lpath[pathi])
    print("CUTPATH : {}\n".format(cutpath))
    graph = remove_paths(graph,cutpath,delete_entry_node,delete_sink_node)
    return graph

def solve_bubble(graph,nanc,ndes):
    """solve_bubble /2 qui prend un graphe, un noeud ancêtre, un noeud descendant et
retourne un graph nettoyé de la bulle se trouvant entre ces deux noeuds en utilisant les
fonctions précédemment développée
    """
    lpath = []
    llenpath = []
    lmeanwei = []
    for path in nx.all_simple_paths(graph,nanc,ndes):
        lpath.append(path)
        llenpath.append(len(path))
        lmeanwei.append(path_average_weight(graph,path))
    graph = select_best_path(graph,lpath,llenpath,lmeanwei)
    return graph

def simplify_bubbles(graph):
    """simplify_bubbles /1 ​ qui prend un graphe et retourne un graphe sans bulle. Il faut ici
identifier le noeud ancêtre et le noeud descendant que l’on va indiquer à ​ solve_bubble​
    """
    start = get_starting_nodes(graph)
    sink = get_sink_nodes(graph)
    for pathi in start:
        for pathj in sink:
            lpath = list(nx.all_simple_paths(graph,pathi,pathj))
            while len(lpath) >= 2:
                i=0
                while lpath[0][i] in lpath[1] and i < len(lpath[0])-1:
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
                while not lpath[0][i+lpath[0].index(bulle)+1] in lpath[1]:
                    i += 1
                bulls = lpath[0][i+lpath[0].index(bulle)+1]
                graph = solve_bubble(graph,bulle,bulls)
                lpath = list(nx.all_simple_paths(graph,pathi,pathj))
    return graph

def solve_entry_tips(graph,lentry):
    """solve_entry_tips /4 qui prend un graphe et une liste de noeuds d’entrée et retourne
graphe sans chemin d’entrée indésirable
    """
    letips = []
    for enod in lentry:
        npred = 1
        snod = enod
        while npred <= 1:
            pred = snod
            snod = list(graph.successors(pred))[0]
            if not bool(list(graph.successors(snod))):
                break
            npred = len(list(graph.predecessors(snod)))
        letips.append([enod,snod])
    lmeanwei = []
    llenpath = []
    lpath = []
    for tips in letips:
        for path in nx.all_simple_paths(graph,tips[0],tips[1]):
            lmeanwei.append(path_average_weight(graph,path))
            llenpath.append(len(path))
            lpath.append(path)
    graph = select_best_path(graph,lpath,llenpath,lmeanwei,True)
    return graph

def solve_out_tips(graph,lout):
    """solve_out_tips /4 qui prend un graphe et une liste de noeuds de sortie et retourne
graphe sans chemin de sortie indésirable
    """
    lotips = []
    for onod in lout:
        nsuc = 1
        pnod = onod
        while nsuc <= 1:
            suiv = pnod
            pnod = list(graph.predecessors(suiv))[0]
            nsuc = len(list(graph.successors(pnod)))
            if not bool(list(graph.predecessors(pnod))):
                break
        lotips.append([pnod,onod])
    lmeanwei = []
    llenpath = []
    lpath = []
    for tips in lotips:
        for path in nx.all_simple_paths(graph,tips[0],tips[1]):
            lmeanwei.append(path_average_weight(graph,path))
            llenpath.append(len(path))
            lpath.append(path)
    graph = select_best_path(graph,lpath,llenpath,lmeanwei,False,True)
    return graph

def main():
    """
    Main program function
    """
    # Get arguments
    args = get_arguments()
    dicc = build_kmer_dict(args.fastq_file, args.kmer_size)
    print("STEP1")
    graph = build_graph(dicc)
    print("STEP2")
    st_nod = get_starting_nodes(graph)
    sk_nod = get_sink_nodes(graph)
    print("STEP3")
    graph = simplify_bubbles(graph)
    print("STEP4")
    graph = solve_entry_tips(graph,st_nod)
    print("STEP5")
    graph = solve_out_tips(graph,sk_nod)
    print("STEP6")
    contigs = get_contigs(graph,get_starting_nodes(graph),get_sink_nodes(graph))
    print("STEP7")
    save_contigs(contigs,args.output_file)
    print("STEP8")

if __name__ == '__main__':
    main()
