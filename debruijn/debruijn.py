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
from random import randint
import networkx as nx
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
    """
    Check if path is an existing file.
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
    """
    Retrieves the arguments of the program.
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
    """
    Takes file and return sequences yield
    """
    with open(fastq_file, "r") as filin:
        for _ in filin:
            yield next(filin).strip()
            next(filin)
            next(filin)

def cut_kmer(sequence : str, k_mer : int):
    """
    takes a sequence and a k_mer size and return a k_mer yield
    """
    for i in range(len(sequence) - k_mer + 1):
        yield sequence[i : i + k_mer]


def build_kmer_dict(fastq_file : str, k_mer : int):
    """
    takes a fastq file and a k_mer size and return a k_mer dictionnary
    """
    k_mer_dict = {}
    for sequence in read_fastq(fastq_file):
        for kmer in cut_kmer(sequence, k_mer):
            if kmer not in k_mer_dict:
                k_mer_dict[kmer] = 0
            k_mer_dict[kmer] += 1


    return k_mer_dict


def build_graph(k_mer_dict : dict):
    """
    Build De Bruijn Graph based on k_mer dictionnary
    """
    bruijn_graph = nx.DiGraph()

    for k_mer, weight in k_mer_dict.items():
        bruijn_graph.add_edge(k_mer[:-1], k_mer[1:], weight = weight)

    return bruijn_graph

def get_starting_nodes(bruijn_graph : nx.DiGraph):
    """
    Take a graph and return a list of starting nodes
    """
    nodes = list(bruijn_graph.nodes())
    nodes_in = []

    for node in nodes:
        pred_len = len(list(bruijn_graph.predecessors(node)))
        if pred_len == 0:
            nodes_in.append(node)

    return nodes_in


def get_sink_nodes(bruijn_graph : nx.DiGraph):
    """
    Take a graph and return a list of sinking nodes
    """
    nodes = list(bruijn_graph.nodes())
    nodes_out = []

    for node in nodes:
        succ_len = len(list(bruijn_graph.successors(node)))
        if succ_len == 0:
            nodes_out.append(node)

    return nodes_out


def get_contigs(bruijn_graph : nx.DiGraph, nodes_in : list, nodes_out : list):
    """
    Get contigs from the files and the starting and sinkng nodes
    """
    contigs = []
    for node_in in nodes_in:
        for node_out in nodes_out:
            if list(nx.all_simple_paths(bruijn_graph, node_in, node_out)):
                liste = list(nx.all_simple_paths(bruijn_graph, node_in, node_out))
                contig = liste[0][0]
                for kmer in range(1, len(liste[0])):
                    contig += liste[0][kmer][-1]
                contigs.append((contig, len(contigs)))

    return contigs

def save_contigs(tuple_list : list, file_name : str):
    """
    Save the contigs in a fasta file
    """
    with open(file_name, "w") as filout:
        for index, contig in enumerate(tuple_list):
            filout.write("> contig_" + str(index + 1) + " len = " + str(contig[1]) + "\n")
            filout.write(fill(contig[0]))
            filout.write("\n")


def fill(text : str, width = 80):
    """
    Split text with a line return to respect fasta format
    """

    return os.linesep.join(text[i:i+width] for i in range(0, len(text), width))


def std(value_list : list):
    """
    Return the standard deviation of an input list
    """
    standard_deviation = statistics.stdev(value_list)
    return standard_deviation

def path_average_weight(bruijn_graph : nx.DiGraph, path : str):
    """
    Return the average weight of a path
    """
    tot_weight = 0
    path_length = len(path) - 1

    for index in path:

        tot_weight += bruijn_graph.out_degree(index, weight = "weight")

    av_weight = tot_weight / path_length

    return av_weight



def remove_paths(bruijn_graph : nx.DiGraph, path_list : list,
    delete_entry_node : bool,
    delete_sink_node : bool):
    """
    Remove paths from the graph.
    Delete entry nodes and/or sink nodes when set to True
    """
    for path in path_list:
        bruijn_graph.remove_nodes_from(path[1:-1])

        if delete_sink_node:
            bruijn_graph.remove_node(path[-1])

        if delete_entry_node:
            bruijn_graph.remove_node(path[0])

    return bruijn_graph

def select_best_path(bruijn_graph : nx.DiGraph, path_list : list,
    path_length : list, medium_path_weight : list, delete_entry_node : bool,
    delete_sink_node : bool):
    """
    Select the best path from the graph
    """
    weight_max = max(medium_path_weight)
    length_max = max(path_length)

    best_path_weight_w = []
    best_path_length_w = []
    best_path_candidates_w = []

    for index, value in enumerate(path_list):

        if value == weight_max:

            best_path_weight_w.append(medium_path_weight[index])
            best_path_length_w.append(path_length[index])
            best_path_candidates_w.append(path_list[index])

    best_path_candidates = []

    for index, value in enumerate(best_path_candidates_w):

        if value == length_max:

            best_path_candidates += best_path_candidates_w[index]


    if len(best_path_candidates) >= 1:

        best_path_index = randint(0, len(best_path_candidates))


    to_remove_paths = []


    for path in best_path_candidates:
        if path is not best_path_candidates[best_path_index]:
            to_remove_paths.append(path)


    bruijn_graph = remove_paths(bruijn_graph, to_remove_paths,
        delete_entry_node = False, delete_sink_node = False)


    return bruijn_graph


def solve_bubble(bruijn_graph : nx.DiGraph, ancestral_node : str, descendant_node : str):
    """
    Remove a bubble from a graph
    """
    bubble_path = []
    bubble_length = []
    bubble_weight = []

    for path in nx.all_simple_paths(bruijn_graph, ancestral_node, descendant_node):
        bubble_path.append(path)
        bubble_length.append(len(path))
        bubble_weight.append(path_average_weight(bruijn_graph, path))

    bruijn_graph = select_best_path(bruijn_graph, bubble_path, bubble_length, bubble_weight,
        delete_entry_node = False, delete_sink_node = False)

    return bruijn_graph


def simplify_bubbles(bruijn_graph : nx.DiGraph):
    """
    Remove all bubbles from a graph
    """
    nodes_to_delete = []

    for descendant_node in bruijn_graph:
        anc_node = list(bruijn_graph.predecessors(descendant_node))

        if len(anc_node) >= 2:
            ancestor = nx.lowest_common_ancestor(bruijn_graph, anc_node[0], anc_node[1])
            nodes_to_delete.append([ancestor, descendant_node])

    for couple_anc_des in nodes_to_delete:
        bruijn_graph = solve_bubble(bruijn_graph, couple_anc_des[0], couple_anc_des[1])

    return bruijn_graph


def solve_entry_tips(bruijn_graph : nx.DiGraph, in_node_list : list):
    """
    Remove all unwanted starting nodes
    """
    path_list = []
    path_length = []
    path_weight = []

    for node in in_node_list:

        des_node = list(bruijn_graph.successors(node))[0]
        anc_node = list(bruijn_graph.predecessors(des_node))
        length_succ = len(list(bruijn_graph.successors(des_node)))
        while len(anc_node) == 1 and length_succ:

            des_node = list(bruijn_graph.successors(des_node))[0]
            anc_node = list(bruijn_graph.predecessors(des_node))

        if len(anc_node) > 1:

            for path in nx.all_simple_paths(bruijn_graph, node, des_node):
                path_list.append(path)
                path_length.append(path)
                path_weight.append(path_average_weight(bruijn_graph, path))

    bruijn_graph = select_best_path(bruijn_graph, path_list, path_length, path_weight,
        delete_entry_node = True, delete_sink_node = False)

    return bruijn_graph

def solve_out_tips(bruijn_graph : nx.DiGraph, out_node_list : list):
    """
    Remove all unwanted sinking nodes
    """
    path_list = []
    path_length = []
    path_weight = []

    for node in out_node_list:

        anc_node = list(bruijn_graph.predecessors(node))[0]
        des_node = list(bruijn_graph.successors(anc_node))
        length_anc = len(list(bruijn_graph.predecessors(anc_node)))
        while len(des_node) == 1 and length_anc:

            anc_node = list(bruijn_graph.predecessors(anc_node))[0]
            des_node = list(bruijn_graph.successors(anc_node))

        if len(des_node) > 1:

            for path in nx.all_simple_paths(bruijn_graph, anc_node, node):
                path_list.append(path)
                path_length.append(path)
                path_weight.append(path_average_weight(bruijn_graph, path))

    bruijn_graph = select_best_path(bruijn_graph, path_list, path_length, path_weight,
        delete_entry_node = False, delete_sink_node = True)

    return bruijn_graph

#==============================================================
# Main program
#==============================================================

def main():
    """
    Main program function
    """
    # Get arguments
    args = get_arguments()

    # Build k_mer_dict
    k_mer_dict = build_kmer_dict(args.fastq_file, args.kmer_size)

    # Build associated graph
    bruijn_graph = build_graph(k_mer_dict)

    # Remove all bubbles
    bruijn_graph = simplify_bubbles(bruijn_graph)

    # Remove non opti in and out tips
    bruijn_graph = solve_entry_tips(bruijn_graph, get_starting_nodes(bruijn_graph))
    bruijn_graph = solve_out_tips(bruijn_graph, get_sink_nodes(bruijn_graph))

    # Get contigs, and save it in fasta file
    contigs = get_contigs(bruijn_graph, get_starting_nodes(bruijn_graph),
        get_sink_nodes(bruijn_graph))
    save_contigs(contigs, args.output_file)


if __name__ == '__main__':
    main()
