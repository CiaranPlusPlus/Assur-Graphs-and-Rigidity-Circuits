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

Assur Graphs will be created using a calculated number of rigidity
graphs that can then be modified to have a set of pinned vertices.

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
import itertools

import networkx as nx
import networkx.algorithms.isomorphism as iso

from os import listdir
from os.path import isfile, join


"""  ---------------- Global Constants ----------------  """

DIR = "AssurGraphFiles/"

ASSUR_DIR = DIR + "AssurGraphs/"
CIRCUIT_DIR = DIR + "RigidityCircuits/"

ASSUR_FILENAME = "AssurGraphsOf_"
CIRCUIT_FILENAME = "RigidityCircuitsOf_"


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
def two_sum(ga, gb, a, b):
    # copying graphs avoids risk of modifying the original object
    grapha = ga.copy()
    graphb = gb.copy()

    # create lists for vertex label mapping
    mapa = []
    mapb = []

    # relabel tuple edge information to relabel nodes to be unionised
    a = ("G1_" + str(a[0]), "G1_" + str(a[1]))
    b = ("G2_" + str(b[0]), "G2_" + str(b[1]))

    # relabel the nodes in each graph to keep each set unique
    # creates a mapping for each graph and then applies it so each
    # vertex has a G#_ as a prefix to the original vertex numbers
    for num1 in list(grapha.nodes):
        mapa.append("G1_" + str(num1))

    for num2 in list(graphb.nodes):
        mapb.append("G2_" + str(num2))

    g1_mapping = dict(zip(grapha, mapa))
    g1 = nx.relabel_nodes(grapha, g1_mapping)

    # Composing two graphs does not remove edges so we need to do it manually
    g1.remove_edge(a[0], a[1])

    # These two lines label the nodes to be joined at 'a' and 'b' so a compose of the two
    # graphs will result in a two sum of the graphs
    g1_union_nodes = {a[0]: 'a', a[1]: 'b'}
    g1 = nx.relabel_nodes(g1, g1_union_nodes, copy=False)

    g2_mapping = dict(zip(graphb, mapb))
    g2 = nx.relabel_nodes(graphb, g2_mapping)

    # Composing two graphs does not remove edges so we need to do it manually
    g2.remove_edge(b[0], b[1])

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
# There are a few basic known Assur Graphs. The AG on three vertices is the
# smallest AG possible. The complete graph of 4 vertices k_4 is used to create
# all possible rigidity circuits which are in turn used to create assur graphs.
# This function is run only  to replace these graphs in in the list, and to
# create them on the first time running.
def basic_graphs():
    # Create the most basic AG of three vertices where 2 are pinned
    # This hardcoded as it is the exception that cannot be generated from a
    # rigidity circuit.
    ag_3 = nx.Graph()
    ag_3.add_nodes_from([(0, {'pinned': True}), (1, {'pinned': False}), (2, {'pinned': True})])
    ag_3.add_edges_from([(0, 1), (1, 2)])

    if not duplicate_or_isomorphic(ag_3, ASSUR_DIR + ASSUR_FILENAME, "assur"):
        write_graph_to_file(ag_3, ASSUR_DIR + ASSUR_FILENAME)

    # Create the complete graph on 4 vertices and write it to
    # the rigidity circuit files to create more rigidity circuits
    k_4 = nx.complete_graph(4)
    if not duplicate_or_isomorphic(k_4, CIRCUIT_DIR + CIRCUIT_FILENAME, "rigidity"):
        write_graph_to_file(k_4, CIRCUIT_DIR + CIRCUIT_FILENAME)


# create_rigidity_circuits
# This function creates all possible rigidity circuits of the integer input
# argument to the script using two_sum and edge_splits.
def create_rigidity_circuits(rc_size):
    # To create graphs of size N, the edge split has to be performed on graphs of
    # size N-1 since it only ever adds one more vertex. N has to be 5 at the
    # absolute minimum for this section.
    edge_split_file = CIRCUIT_DIR + CIRCUIT_FILENAME + str(rc_size - 1) + ".txt"

    with open(edge_split_file, "r") as graphs_to_edge_split:

        for graph in graphs_to_edge_split:
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

                    if not duplicate_or_isomorphic(new_temp_graph,
                                                   CIRCUIT_DIR + CIRCUIT_FILENAME,
                                                   "rigidity"):
                        write_graph_to_file(new_temp_graph, CIRCUIT_DIR + CIRCUIT_FILENAME)

    # To achieve a correct two sum for the desired graph size, the two circuits being summed
    # together cannot have a total vertex count of size larger than the desired size plus 2.
    # Exceptions can be made since graphs smaller than 4 do not exist for creating Assur Graphs
    # outside of the AG_3 exception, which cannot be made from smaller graphs anyway.
    two_sum_addends = get_set_of_number_combinations(rc_size)

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
                                                               CIRCUIT_DIR + CIRCUIT_FILENAME,
                                                               "rigidity"):
                                    write_graph_to_file(temp_graph_2sum_result,
                                                        CIRCUIT_DIR + CIRCUIT_FILENAME)

                                # Now check the 2-sum with one graph flipped
                                edge_i_flipped = tuple(reversed(edge_i))
                                temp_graph_2sum_result_flip = two_sum(h, i, edge_h, edge_i_flipped)

                                if not duplicate_or_isomorphic(temp_graph_2sum_result_flip,
                                                               CIRCUIT_DIR + CIRCUIT_FILENAME,
                                                               "rigidity"):
                                    write_graph_to_file(temp_graph_2sum_result,
                                                        CIRCUIT_DIR + CIRCUIT_FILENAME)


# create_assur_graphs
# This function takes in the desired Assur Graph size and will generate all Assur Graphs
# for that size. It does this by using the rigidity circuits that are of size less than
# this. For each graph, each node of that graph will be 'split' such that the resulting
# graph is of the desired size. The function does this by checking for how many nodes
# it needs to add to the initial rigidity circuit and then creates resulting pinned
# graphs for every possible combination of the edges that connected to the node that
# was split.
def create_assur_graphs(ag_size):
    files = []

    # To create graphs of size N, we need to look at all rigidity circuits that are
    # smaller than the desired assur graph size so that we can split out nodes
    # and pin them.
    for x in range(1, ag_size):
        temp_name = CIRCUIT_FILENAME + str(x) + ".txt"

        # If that generated file name exists, it will be added to the list of files to be parsed
        if temp_name in os.listdir(CIRCUIT_DIR):
            files.append(os.path.abspath(CIRCUIT_DIR + temp_name))

    # this loop parses though the rigidity circuits in each of the files and creates
    # Assur Graphs from them, then it it check the destination files for isomorphic
    # graphs and write them if the graph is unique
    for file in files:
        with open(file) as graph_file:
            for graph in graph_file:
                edge_list = ast.literal_eval(graph)
                g = nx.Graph()  # graph object to hold the temporary graph
                g.add_edges_from(edge_list[0])
                nodes = list(nx.nodes(g))

                # The number of nodes to add to get the desired size
                nodes_to_add = ag_size - len(g)

                for node in nodes:
                    # this creates a list of all possible Assur Graphs for this rigidity circuit
                    assur_graphs = split_vertex_combinations(g, node, nodes_to_add)

                    # check for isometric graphs in the file before writing
                    for ag in assur_graphs:
                        if not duplicate_or_isomorphic(ag, ASSUR_DIR + ASSUR_FILENAME, "assur"):
                            write_graph_to_file(ag, ASSUR_DIR + ASSUR_FILENAME)


# split_vertex_combinations
# This takes in a graph, a node on that graph, and the number of vertices to add. Using
# this information, it will split the node into the set number of vertices and then for
# each combination of distributed edges it will create a graph with a set of pinned
# and unpinned edges and return the list of all of the graphs created this way. This
# should result in a list of Assur Graphs with some possible isomorphic entries.
def split_vertex_combinations(h, n, m):
    a_graphs = []  # list of graphs being created by this function
    pinned_node_list = [n]  # this is a list of the split nodes
    i = h.copy()  # Copy the input graph to avoid original object changes

    # this loop adds in new nodes to the graph so it meets the size requirements
    for x in range(0, m):
        # The next node will be one integer higher than what is already there, nodes start at '0'
        enum = len(i)

        # add a new node labelled sequentially higher than the ones already in the graph
        i.add_node(enum)
        pinned_node_list.append(enum)

    # this sets a bound for how many combinations of edges there can be
    valence = i.degree[n]

    # If the valence of the node being 'split' is equal to the number of
    # nodes being added then there is not enough edges for all of these
    # nodes and the graph cannot be split as there cannot be isolates.
    if valence <= m:
        a_graphs = []
        return a_graphs

    # This sets each node to have a true or false pinned status
    attributes_dict = {}
    for node in i.nodes:
        if node in pinned_node_list:
            attributes_dict[node] = {'pinned': True}
        else:
            attributes_dict[node] = {'pinned': False}

    nx.set_node_attributes(i, attributes_dict)  # apply the node attributes

    # get all combinations of edges distributed among the pinned vertices
    set_of_edge_combinations = edge_combinations(list(i.edges(n)), pinned_node_list)

    # remove the edges that will be replaced from the original graph
    edges_to_remove = list(i.edges(n))
    i.remove_edges_from(edges_to_remove)

    for combo in set_of_edge_combinations:
        j = i.copy()  # copy i so the loop does not affect the base graph object
        j.add_edges_from(combo)  # add the set of edges for this iteration
        a_graphs.append(j)

    return a_graphs


# edge_combinations
# This takes in the list of edges adjacent to the vertex being split, it then creates
# a list of all possible combinations of edge connections between the set of pinned
# vertices, leaving no pinned vertex unconnected to avoid isolates. This functions works
# by getting the list of edges connected to the node of interest, and then creates a full
# list of all combinations, replacing the number that represents the node of interest in
# each edge connection data structure with one of the nodes that is pinned. It returns a
# list of lists that contains the pinned nodes edges that can just be applied to a graph
# data structure where the old edges have been removed. This gets all possible Assur Graphs
# for that graph configuration, even isometric entries.
def edge_combinations(list_of_edges, list_of_nodes):
    number_combinations = []
    combinations = []
    edge_count = len(list_of_edges)

    # To get all useful combinations without leaving an isolate, we get the product
    # of the list repeated for the number of edges that need to be distributed. We
    # then check for isolated vertices and ignore these combinations as they will not
    # result in valid Assur Graphs
    for n in itertools.product(list_of_nodes, repeat=edge_count):
        result = set(list_of_nodes).issubset(n)  # This makes sure there are no isolates

        if result:
            number_combinations.append(n)

    # This nested loop will create the set of edge connection combinations
    for combo in number_combinations:
        temp_list = []

        for e, edge in enumerate(list_of_edges):
            new_edge = combo[e], edge[1]  # create the new edge combination
            temp_list.append(new_edge)

        combinations.append(temp_list)  # add this set of combinations to the list to return
    return combinations


# duplicate_or_isomorphic
# this accepts a graph as an argument and returns a boolean result based on
# whether that graph has already been created and stored. It will check for
# a duplicate graph first before checking for isomorphism between graphs.
def duplicate_or_isomorphic(graph_to_check, graph_type_filename, mode):
    filename = graph_type_filename + str(len(graph_to_check)) + ".txt"

    # create the file if it does not already exist
    if not os.path.exists(filename):
        with open(filename, "w"):
            pass

    with open(filename, "r+") as temp_file:
        for line_number, line in enumerate(temp_file):

            try:
                temp_graph = nx.Graph()
                temp_list = ast.literal_eval(line)

                temp_graph.add_edges_from(temp_list[0])
                temp_graph.add_nodes_from(temp_list[1])

                # Since graphs are stored slightly differently the Assur Graphs themselves
                # require a node match check for their pinned vertices. While rigidity
                # circuits do not use any status for the nodes as it is not needed.
                if mode == "assur":
                    nm = iso.categorical_node_match('pinned', [True, False])
                    # if a graph is isomorphic then it is a duplicate entry
                    if nx.is_isomorphic(temp_graph, graph_to_check, node_match=nm):
                        return True

                elif mode == "rigidity":
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


# write_graph_to_file
# Creates a list format of the set of edges and node information of the graph.
# this is useful because we can then reconstruct the graphs from this format,
# and it is human readable from the text files they are stored in.
def write_graph_to_file(graph_to_write, graph_type_filename):
    filename = graph_type_filename + str(len(graph_to_write)) + ".txt"

    # create the file if it does not already exist
    if not os.path.exists(filename):
        with open(filenamen, 'w'):
            pass

    edge_list = list(graph_to_write.edges)
    node_info_list = list(graph_to_write.nodes.data())

    graph = [edge_list, node_info_list]

    with open(filename, "a") as temp_file:
        temp_file.write(str(graph) + "\n")


# validate_files
# This function checks to see if all directories and sub-directories
# exist for containing the files that host the set of graphs for each
# vertex count. It uses the program input to check if all files
# below that desired number exist and replace them if they do not.
def validate_files(input_num):
    # Check to see if the directory for the graph files exists
    try:
        # Create target file Directory
        os.mkdir(DIR)
        print("Directory " + DIR + " Created.\n")
    except Exception:
        print("Directory " + DIR + " already exists \n")

    try:
        # Create target rigidity circuit Directory
        os.mkdir(CIRCUIT_DIR)
        print("Directory " + CIRCUIT_DIR + " Created.\n")
    except Exception:
        print("Directory " + CIRCUIT_DIR + " already exists \n")

    try:
        # Create target assur graph Directory
        os.mkdir(ASSUR_DIR)
        print("Directory " + ASSUR_DIR + " Created.\n")
    except Exception:
        print("Directory " + ASSUR_DIR + " already exists \n")

    # Check for files for the rigidity circuits. These files need to exist
    # to create the circuits the Assur Graphs are generated from.
    for rc_num in range(4, input_num):
        # if the input num is 3 or 4 just make sure the basic graph set is done
        if rc_num <= 4:
            basic_graphs()
            continue

        file_name_rc_num = CIRCUIT_DIR + CIRCUIT_FILENAME + str(rc_num) + ".txt"
        if not os.path.exists(file_name_rc_num):
            # if the files are not there, create them
            print("Creating sub-files for rigidity circuits of size " + str(rc_num) + "...\n")
            create_rigidity_circuits(rc_num)


# get_set_of_number_combinations
# This function takes in a number and finds the list of all possible two
# number combinations that can add up to a certain number. This is useful
# for determining what graphs to 2-sum together to find a resulting graph
# of a desired size. For example, getting 2-summed graphs of total vertex
# count 10 requires the list of all two summed numbers that result in 12,
# since two vertices are "lost" in the 2-sum operation.
def get_set_of_number_combinations(max_num):
    number_list = []

    # when performing a 2-sum 2 vertices will be lost from the total as they
    # are "glued" together in a sense, so the max number has an offset of 2
    offset_max = max_num + 2
    cap = round(offset_max/2)+1

    # it starts at 4 because that is the smallest rigidity circuit available.
    for num in range(4, cap):
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
                        help='The largest Assur Graph of total vertex count n to be generated. '
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

    if arg.integer[0] <= 2:
        print("Invalid input of " + str(arg.integer[0]) + ". Use an input of n > 2")
        exit()

    if arg.integer[0] == 0:
        print("No Assur Graphs on 4 vertices exist.")
        exit()

    return arg.integer[0]


"""  --------------- Main Call ----------------  """


# main
# This is the top most function being called. It functions as the flow control
# of all other functions. It is called below the Main Call line.
def main():
    # Get the input arguments supplied to the program on the command line
    arguments = apply_arg_parser()
    vertex_count = handle_args(arguments)

    # if he input was larger than the current set of files that exist
    # this method will create the smaller rigidity circuit files
    validate_files(vertex_count)

    # This checks the size of the input. At less than 3 the arg_parser will fail
    # At 3, the basic set of graphs will be created. At 5+ it will create those
    # sets of Assur Graphs. Also
    if vertex_count >= 4:
        create_assur_graphs(vertex_count)
    else:
        basic_graphs()


# General execution catch error handling.
try:
    main()
except Exception as mainExecuteError:
    print(mainExecuteError)
