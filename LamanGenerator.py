# Imports
import argparse
import ast
import os

import networkx as nx


# Constants
DIR = "LamanGraphFiles"
FILENAME = DIR + "/LamanGraphsOf"


# Definitions

# hen_1:
# Adds a new vertex and connects it to two previously existing vertices.
# Returns a new graph created by this operation
def hen_1(g, edge):
    a, b = edge
    h = g.copy()  # maintains the original used graph integrity
    enum = len(g)  # Arrays start at zero, so the size will return the next integer to be used as a node.

    # create an edge between the new node and the two nodes the edge connects to
    h.add_edge(a, enum)
    h.add_edge(b, enum)
    return h


# hen_2:
# Splits an edge, adds a new vertex, and connects it to the two vertices
# that the edge was removed from, as well as a third pre-existing vertex.
# Returns a new graph created by this operation
def hen_2(g, a, b, c):
    enum = len(g)  # Arrays start at zero, so the size will return the next integer to be used as a node.
    h = g.copy()  # maintains the original used graph integrity

    h.remove_edge(a, b)  # Remove the edge from the nodes supplied
    h.add_node(enum)  # add new node to the graph
    # add edges from this new node to the previous two nodes plus the third supplied one
    h.add_edge(a, enum)
    h.add_edge(b, enum)
    h.add_edge(c, enum)

    return h


# do_files_exist:
# Checks for the existence of the Laman Graph files that will hold the set
# of all Laman graphs that exist for a certain number of vertices.
# Returns a boolean state and a list of the files missing for graph generation
def do_files_exist(n):
    missing_files = []
    x = 2

    # check to see if previous graphs exist so the new ones can be created
    while x < n:
        if not os.path.isfile(FILENAME + str(x) + ".txt"):
            missing_files.append(x)

        x += 1

    if len(missing_files) > 0 and x >= 2:  # if there are missing files
        state = False
    else:
        state = True

    return state, missing_files


# create_laman_graph_set
# Creates a set of laman graphs for input integer n. If n is 2 it will create
# the smallest Laman Graph of 2 vertices and a single edge. Because of the use
# of previously generated graphs to construct higher order graphs, this is an
# edge case and has it's own functionality for when n = 2.
# Returns the file name of the set of graphs and the set of Laman graphs
def create_laman_graph_set(n):
    print("def create_laman_graph_set for n = " + str(n))
    file_name = FILENAME + str(n) + ".txt"
    set_of_graphs = []

    if n == 2:  # edge case condition of the first Laman graph
        set_of_graphs.append(nx.complete_graph(2))
    else:
        print('create graphs of size ' + str(n) + ' using lower graphs')  # Debug text
        previous_graphs = build_graphs_from_file(FILENAME + str(n-1) + ".txt")  # Load lower order graphs

        for graph in previous_graphs:
            temp_graph = list(nx.edges(graph, nbunch=None))  # Store the edges as a list to be used

            for edge in temp_graph:
                temp_nodes = list(nx.nodes(graph))  # Temporary list of nodes to be used in parsing the graph
                next_laman_hen_1 = hen_1(graph, edge)  # Use the first henneberg operation to create graphs

                if not graph_is_isomorphic(set_of_graphs, next_laman_hen_1):
                    set_of_graphs.append(next_laman_hen_1)

                if len(graph) > 2:  # hen_2 cannot be performed on graphs with < 3 vertices. Arrays start at 0
                    a, b = edge  # hen_2 needs to know what edge to remove and what nodes to reconnect
                    temp_nodes.remove(a)
                    temp_nodes.remove(b)

                    # This will add an edge between the new node and a node not connected to the removed edge
                    # This will be done for each node not connected to the removed edge so all permutations
                    # are accounted for
                    for node in temp_nodes:
                        next_laman_hen_2 = hen_2(graph, a, b, node)  # Use the 2nd henneberg operation to create graphs
                        if not graph_is_isomorphic(set_of_graphs, next_laman_hen_2):
                            set_of_graphs.append(next_laman_hen_2)

    return file_name, set_of_graphs


# graph_is_isomorphic
# Checks in the supplied graphs is an isomorphic of any of the graphs contained
# within the supplied set of graphs
# Returns a True if the graph is an isomorphic, false otherwise
def graph_is_isomorphic(graph_set, graph_to_test):
    if not graph_set:
        return False
    else:
        for finished_graph in graph_set:  # For each graph in this list
            if nx.is_isomorphic(finished_graph, graph_to_test):
                return True


# build_graphs_from_file:
# restores a list of existing graphs based on the file being called.
# Returns the list of graphs created from the file.
def build_graphs_from_file(called_file):
    graphs = []  # will hold the list of graphs from the file

    with open(called_file) as f:
        lines = [line.rstrip("\n") for line in f]

    # Creates a list of NetworkX graphs objects instead of a list of edge connection information
    for L in lines:
        g = nx.Graph()
        temp = list(ast.literal_eval(L))
        g.add_edges_from(temp)
        graphs.append(g)
    return graphs


# recreate_graphs
# If graphs are missing in the list of files for creating higher order
# graphs, then this will attempt to restore those missing files.
def recreate_graphs(missing_graphs):
    for n in missing_graphs:
        created_file, set_of_laman_graphs = create_laman_graph_set(n)
        write_list_to_file(created_file, set_of_laman_graphs)


# write_list_to_file
# Creates a list format of the set of edges since it also shows the node information of the graph
# this is useful because we can then reconstruct the graphs from this format, and it is human
# readable from the text files they are stored in.
def write_list_to_file(filename, list_to_write):

    with open(filename, "w+") as temp_file:

        for element in list_to_write:
            temp = nx.edges(element, nbunch=None)  # Stores the graphs in a list of edge connection information
            temp_file.write(str(temp) + "\n")


# apply_arg_parser
# setup to make sure the use supplies the program with an integer argument when calling the script
# Also has a -h flag to tell a user what the arguments are for and their format
def apply_arg_parser():
    parser = argparse.ArgumentParser(description='Will calculate the Laman Graphs for all vertex counts'
                                                 ' from 2 up until the supplied integer argument n.')

    parser.add_argument('integer', metavar='n', type=int, nargs=1, help='The largest LamanGraph of vertex count'
                                                                        ' n to be generated. A value of 1 will not'
                                                                        ' generate any Laman graphs, and will result in' 
                                                                        ' no output.')

    return parser


# handle_args
# Used to handle any arguments supplied by a user on the command line. Some conditions like
# when the integer input is less than 1 will be accounted for here by just abandoning the
# executing of the script to avoid crashes or wasting time.
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

    return arg.integer[0]


#  --------------- Main ----------------  #

arguments = apply_arg_parser()
vertex_count = handle_args(arguments)

try:
    # Create target Directory
    os.mkdir(DIR)
    print("Directory " + DIR + " Created ")
except FileExistsError:
    print("Directory " + DIR + " already exists")

files_exist, missing_graphs_by_index_count = do_files_exist(int(vertex_count))

if not files_exist:
    for e in missing_graphs_by_index_count:
        print("Missing LamanGraphFiles/LamanGraphOf" + str(e))
    recreate_graphs(missing_graphs_by_index_count)

file, set_of_lgraphs = create_laman_graph_set(vertex_count)
write_list_to_file(file, set_of_lgraphs)
