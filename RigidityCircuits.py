#!/usr/bin/env python

"""
Generate a list of all Assur Graphs.

Author:                     Ciaran Mc Glue
Date of creation start:     25/03/2020

The purpose of this script is to create a comprehensive list of all
possible combinations of Assur Graphs starting with the smallest
'basic' Assur Graphs up to a graph size of N. This program will be
using the definitions and methods set forth in the paper
"Combinatorial Characterization of the Assur Graphs from Engineering"
by Servatius et al.

Graphs of a vertex size N are stored in text files labelled by the
vertex count (cumulative count of pinned and unpinned vertices) in
the format of:

[[(n,n+1), ...], [(n, {'pinned':BOOLEAN}), (n+1, {'pinned':BOOLEAN}), ...]]

This is two lists, the first is the edge set, the second is the vertex
set. The edge set is listed as the pair of connected vertices while the
vertex set maintains information about whether a vertex is pinned or
not using a boolean dictionary key value pair.

This minimises impact to storage space while maintaining a somewhat
human-readable format.
"""

"""  ------------------- Imports --------------------  """

import argparse
import ast
import os

import networkx as nx

from os import listdir
from os.path import isfile, join


"""  ---------------- Global Constants ----------------  """

DIR = "AssurGraphFiles"

CIRCUIT_DIR = DIR + "/RigidityCircuits"
CIRCUIT_FILENAME = "/RigidityCircuitsOf_"


"""  --------------- Graphing Functions ----------------  """


# edge_split:
# Splits an edge, adds a new vertex, and connects it to the two vertices
# that the edge was removed from, as well as a third pre-existing vertex.
# Returns a new graph created by this operation
def edge_split(graph, a, b, c):
    # Arrays start at zero, so the size will return the next integer to be used as a node.
    enum = len(graph)

    # maintains the original used graph integrity
    h = graph.copy()

    h.remove_edge(a, b)  # Remove the edge from the nodes supplied
    h.add_node(enum)  # add new node to the graph

    # add edges from this new node to the previous two nodes plus the third supplied one
    h.add_edge(a, enum)
    h.add_edge(b, enum)
    h.add_edge(c, enum)

    return h


# two_sum
# Takes in an argument of two graphs and the edges of each where the graphs will be
# joined together and then returns the newly created graph. The size of the graph returned
# will be the total sum of the vertices between the two minus 2: (|G(V)| + |H(V)|) - 2.
# Graphs will have their nodes relabelled to make it clear which nodes belong to which graph.
# Then depending on the edges where the sum is going to be the respective nodes are labelled
# the same so that a graph union can be applied to cause the graphs to be "glued" along
# those respective edges. nx.compose does not relabel all edges after composing so they will
# need to be renumbered so graph storage format is maintained.
def two_sum(grapha, graphb, a, b):
    print("two_sum...")
    mapa = []
    mapb = []

    # relabel tuple edge information to relabel nodes to be unionised
    a = ("G1_" + str(a[0]), "G1_" + str(a[1]))
    b = ("G2_" + str(b[0]), "G2_" + str(b[1]))

    # relabel the nodes in each graph to keep each set unique
    # creates a mapping for each graph and then applies it so each
    # vertex has a G#_ as a prefix
    for num1 in range(0, len(grapha)):
        mapa.append("G1_" + str(num1))

    for num2 in range(0, len(graphb)):
        mapb.append("G2_" + str(num2))

    g1_mapping = dict(zip(grapha, mapa))
    g1 = nx.relabel_nodes(grapha, g1_mapping)

    # These two lines label the nodes to be joined at 'a' and 'b' so a compose of the two
    # graphs will result in a two sum of the graphs
    g1_union_nodes = {a[0]: 'a', a[1]: 'b'}
    g1 = nx.relabel_nodes(g1, g1_union_nodes, copy=False)

    g2_mapping = dict(zip(graphb, mapb))
    g2 = nx.relabel_nodes(graphb, g2_mapping)

    # These two lines label the nodes to be joined at 'a' and 'b' so a compose of the two
    # graphs will result in a two sum of the graphs
    g2_union_nodes = {b[0]: 'a', b[1]: 'b'}
    g2 = nx.relabel_nodes(g2, g2_union_nodes, copy=False)

    # compose the two graphs so they will be joined at the common nodes and edge
    u = nx.compose(g1, g2)

    # relabel the graphs from 0 - N so they can be stored correctly
    u = nx.convert_node_labels_to_integers(u, first_label=0, ordering="default")

    return u


# basic_graphs
# Creates the rigidity circuit K_4 and adds it to the file of rigidity circuits on 4 vertices
def basic_graphs():
    print("basic_graphs... \n")

    # Create the complete graph on 4 vertices and write it to
    # the rigidity circuit files to create more rigidity circuits
    k_4 = nx.complete_graph(4)
    if not duplicate_or_isomorphic(k_4, CIRCUIT_DIR + CIRCUIT_FILENAME):
        write_graph_to_file(k_4, CIRCUIT_DIR + CIRCUIT_FILENAME)


# create_rigidity_circuits
# This function creates all possible rigidity circuits of the integer input
# argument to the script using two_sum and edge_splits.
def create_rigidity_circuits(max_size):
    print("creating rigidity circuits...")

    # To create graphs of size N, the edge split has to be performed on graphs of
    # size N-1 since it only ever adds one more vertex. N has to be 5 at the
    # absolute minimum for this section.
    edge_split_file = CIRCUIT_DIR + CIRCUIT_FILENAME + str(max_size - 1) + ".txt"

    with open(edge_split_file, "r") as graphs_to_split:

        for graph in graphs_to_split:
            edge_list = ast.literal_eval(graph)
            g = nx.Graph()
            g.add_edges_from(edge_list[0])
            nodes = list(nx.nodes(g))

            for edge in edge_list[0]:
                # copy list so items can be removed without affecting the source list
                temp_nodes = nodes.copy()

                a, b = edge
                temp_nodes.remove(a)
                temp_nodes.remove(b)

                for node in temp_nodes:
                    new_temp_graph = edge_split(g, a, b, node)

                    if not duplicate_or_isomorphic(new_temp_graph, CIRCUIT_DIR + CIRCUIT_FILENAME):
                        write_graph_to_file(new_temp_graph, CIRCUIT_DIR + CIRCUIT_FILENAME)

    # To achieve a correct two sum for the desired graph size, the two circuits being summed
    # together cannot have a total vertex count of size larger than the desired size plus 2.
    # Exceptions can be made since graphs smaller than 4 do not exist for creating Assur Graphs
    # outside of the AG_3 exception, which cannot be made from smaller graphs anyway.
    two_sum_addends = get_set_of_number_combinations(max_size)

    for addends in two_sum_addends:
        num1, num2 = addends
        graph_file_1 = CIRCUIT_DIR + CIRCUIT_FILENAME + str(num1) + ".txt"
        graph_file_2 = CIRCUIT_DIR + CIRCUIT_FILENAME + str(num2) + ".txt"

        with open(graph_file_1, "r") as file_1:
            with open(graph_file_2, "r") as file_2:

                # These nested loops will create 2-sums for all graphs contained in the files
                # For each graph in the first file it will perform a 2 sum operation for each
                # graph in the second. A 2-sum will be performed for each edge of each graph
                # so this method is a brute-force approach and could become time-expensive
                # for large sets of graphs. Only unique graphs that no isomorphisms already
                # exist for will be added to the file of all graphs of that length
                for lines_1 in file_1:
                    edge_list_1 = ast.literal_eval(lines_1)
                    h = nx.Graph()
                    h.add_edges_from(edge_list_1[0])

                    for lines_2 in file_2:
                        edge_list_2 = ast.literal_eval(lines_2)
                        i = nx.Graph()
                        i.add_edges_from(edge_list_2[0])

                        for edge_h in edge_list_1[0]:
                            for edge_i in edge_list_2[0]:
                                temp_graph_2sum_result = two_sum(h, i, edge_h, edge_i)

                                if not duplicate_or_isomorphic(temp_graph_2sum_result,
                                                               CIRCUIT_DIR + CIRCUIT_FILENAME):
                                    write_graph_to_file(temp_graph_2sum_result,
                                                        CIRCUIT_DIR + CIRCUIT_FILENAME)


# duplicate_or_isomorphic
# this accepts a graph as an argument and returns a boolean result based on
# whether that graph has already been created and stored. It will check for
# a duplicate graph first before checking for isomorphism between graphs.
def duplicate_or_isomorphic(graph_to_check, graph_type_filename):
    filename = graph_type_filename + str(len(graph_to_check)) + ".txt"

    if not os.path.exists(filename):
        with open(filename, "w"): pass

    with open(filename, "r+") as temp_file:
        for line_number, line in enumerate(temp_file):

            try:
                temp_graph = nx.Graph()
                temp_list = ast.literal_eval(line)

                temp_graph.add_edges_from(temp_list[0])
                temp_graph.add_nodes_from(temp_list[1])

                # if a graph is isomorphic then it is a duplicate entry
                if nx.is_isomorphic(temp_graph, graph_to_check):
                    return True

            except Exception as error:
                print("Error with graph line: " + str(error))
                print("Error line string is: \n" + line)
                # it would be useful here to delete the error line and
                # then just continue the parsing of the file but it
                # is both difficult and possibly dangerous to edit a
                # file while parsing it. Easier to remove the error line
                # manually after being alerted by an error handler.

    return False




"""  --------------- Utility Functions ----------------  """


# write_list_to_file
# Creates a list format of the set of edges and node information of the graph.
# this is useful because we can then reconstruct the graphs from this format,
# and it is human readable from the text files they are stored in.
def write_graph_to_file(graph_to_write, graph_type_filename):
    filename = graph_type_filename + str(len(graph_to_write)) + ".txt"

    if not os.path.exists(filename):
        with open(filenamen, 'w'): pass

    edge_list = list(graph_to_write.edges)
    node_info_list = list(graph_to_write.nodes.data())

    graph = [edge_list, node_info_list]

    with open(filename, "a") as temp_file:
        temp_file.write(str(graph) + "\n")


# get_files
# returns a list of the files in the supplied directory, list will
# need to be cast and parsed somewhere else
def get_files(directory):
    files_info = [os.path.abspath(directory + f) for f in os.listdir(directory)]

    return files_info


# validate_files
#
def validate_files():
    # Check to see if the directory for the graph files exists
    try:
        # Create target Directory
        os.mkdir(DIR)
        print("Directory " + DIR + " Created.\n")
    except Exception:
        print("Directory " + DIR + " already exists \n")

    try:
        # Create target Directory
        os.mkdir(CIRCUIT_DIR)
        print("Directory " + CIRCUIT_DIR + " Created.\n")
    except Exception:
        print("Directory " + CIRCUIT_DIR + " already exists \n")


# get_set_of_number_combinations
# This function takes in a number and finds the list of all possible two
# number combinations that can add up to a certain number. This is useful
# for determining what graphs to 2-sum together to find a resulting graph
# of a desired size. For example, getting 2-summed graphs of total vertex
# count 10 requires the list of all two summed numbers that result in 12,
# since two vertices are "lost" in the 2-sum operation.
def get_set_of_number_combinations(max):
    print("getting number combinations")
    number_list = []

    # when performing a 2-sum 2 vertices will be lost from the total as they
    # are "glued" together in a sense, so the max number has an offset of 2
    offset_max = max + 2

    # it starts at 4 because that is the smallest rigidity circuit available
    for num in range(4, offset_max):
        result = offset_max - num

        if result >= 4 and num >= 4:
            combo = num, result
            number_list.append(combo)

    return number_list


# apply_arg_parser
# setup to make sure the use supplies the program with an integer argument
# when calling the script Also has a -h flag to tell a user what the arguments
# are for and their format.
def apply_arg_parser():
    parser = argparse.ArgumentParser(description='Will calculate the Assur Graphs for all inner '
                                                 'vertex counts from 3 up until the supplied '
                                                 'integer argument n.')

    parser.add_argument('integer',
                        metavar='n',
                        type=int,
                        nargs=1,
                        help='The largest Assur Graph of inner vertex count n to be generated. '
                             'A value less than 3 will not generate any Assur graphs, and will '
                             'result in no output.')

    return parser


# handle_args
# Used to handle any arguments supplied by a user on the command line.
# Some conditions like when the integer input is less than 1 will be
# accounted for here by just abandoning the executing of the script to
# avoid crashes or wasting time.
def handle_args(arg_parser):
    arg, unknown = arg_parser.parse_known_args()
    unknowns = ""

    for unknown_arg in unknown:
        unknowns = unknowns + unknown_arg + " "

    if unknowns is not "":
        print("Ignoring unrecognised arguments: " + unknowns)
        print("Using first argument: " + str(arg.integer[0]))

    if arg.integer[0] <= 3:
        print("Invalid input of " + str(arg.integer[0]) + ". Use an input of n > 3")
        exit()

    return arg.integer[0]


"""  --------------- Main Call ----------------  """

# main
# This is the top most function being called. It functions as the control of
# all other functions. It is called below the Main Call line.
def main():
	# Get the input arguments supplied to the program on the command line
    arguments = apply_arg_parser()
    vertex_count = handle_args(arguments)

    validate_files()
	basic_graphs()
	create_rigidity_circuits(vertex_count)

try:
    main()
except Exception as mainExecuteError:
    print(mainExecuteError)
