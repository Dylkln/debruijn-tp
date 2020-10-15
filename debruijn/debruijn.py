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
import random
import statistics
import networkx as nx
import matplotlib.pyplot as plt
from operator import itemgetter
from random import randint
random.seed(9001)


__author__ = "Dylan KLEIN"
__copyright__ = "Universite Paris Diderot"
__credits__ = ["Dylan KLEIN"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Dylan KLEIN"
__email__ = "klein.dylan@outlook.com"
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


def read_fastq(fastq_file : str):

    with open(fastq_file, "r") as filin:
        for _ in filin:
            yield next(filin).strip()
            next(filin)
            next(filin)

def cut_kmer(sequence : str, k_mer : int):
    
    for i in range(len(sequence) - k_mer + 1):
        yield sequence[i : i + k_mer]


def build_kmer_dict(fastq_file : str, k_mer : int):

    k_mer_dict = {}
    
    for sequence in read_fastq(fastq_file):
        for kmer in cut_kmer(sequence, k_mer):
            if kmer not in k_mer_dict.keys():
                k_mer_dict[kmer] = 0
            k_mer_dict[kmer] += 1


    return k_mer_dict


def build_graph(k_mer_dict : dict):
    
    graph = nx.DiGraph()
    
    for k_mer, weight in k_mer_dict.items():
        graph.add_edge(k_mer[:-1], k_mer[1:], weight = weight)

    plt.subplot(222)

    nx.draw(graph, with_labels = True)
    plt.savefig("graph")

    return graph

def get_starting_nodes(graph):
    pass


def get_sink_nodes(graph):
    pass


def get_contigs(graph, nodes_in : list, nodes_out : list):
    pass


def save_contigs(tuple_list : list, file_name : str):
    pass


def fill(text, width=80):
"""Split text with a line return to respect fasta format"""

return os.linesep.join(text[i:i+width] for i in range(0, len(text), width))


def std(valeur : int):
    pass


def path_average_weight(graph, path):
    pass


def remove_paths(graph, path_list : list):
    pass


def select_best_path(graph, path_list : list, path_length : list,
                        medium_path_weight : list):
    pass


def solve_bubble(graph, ancestral_node, descendant_node):
    pass


def simplify_bubbles(graph):
    pass


def solve_entry_tips(graph, in_node_list : list):
    pass


def solve_out_tips(graph, out_node_list : list):
    pass    


#==============================================================
# Main program
#==============================================================
def main():
    """
    Main program function
    """
    # Get arguments

    args = get_arguments()
    k_mer_dict = build_kmer_dict(args.fastq_file, args.kmer_size)
    graph = build_graph(k_mer_dict)


if __name__ == '__main__':
    main()
