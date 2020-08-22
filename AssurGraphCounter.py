#!/usr/bin/env python

"""
Count the number of graphs in each file for both the Assur Graphs
and the rigidity circuits..

Author:                     Ciaran Mc Glue
Date of creation start:     19/04/2020

"""

"""  ------------------- Imports --------------------  """

import os, re

import networkx as nx

from os import listdir
from os.path import isfile, join


"""  ---------------- Global Constants ----------------  """


DIR = "AssurGraphFiles/"

ASSUR_DIR = DIR + "AssurGraphs/"
CIRCUIT_DIR = DIR + "RigidityCircuits/"

ASSUR_FILENAME = "AssurGraphsOf_"
CIRCUIT_FILENAME = "RigidityCircuitsOf_"


"""  --------------- Utility Functions ----------------  """


# get_files
# returns a list of the files in the supplied directory, list will
# need to be cast and parsed somewhere else
def get_files(directory):
    files_info = [os.path.abspath(directory + f) for f in os.listdir(directory)]

    return files_info


# unit_testing
# Place in method calls here with instanced variables and then call
# this method from main to test partitioned behaviour
def unit_testing():
    print("testing... \n")


"""  --------------- Main Call ----------------  """


# main
# This is the top most function being called. It functions as the control of
# all other functions. It is called below the Main Call line.
def main():
    num = None

    rigidty_files = get_files(CIRCUIT_DIR)
    assur_files = get_files(ASSUR_DIR)

    print("Rigidity Circuits:")
    print("Vertex # : Number of Graphs")

    for r_file in rigidty_files:
        r_counter = 0
        file_num = re.search('_(\\d).txt', r_file)

        if file_num:
            num = file_num.group(1)

        with open(r_file, "r") as file_to_count:
            for line in file_to_count:
                r_counter += 1

        print(str(num) + " : " + str(r_counter))

    print("\nAssur Graphs:")
    print("Vertex # : Number of Graphs")

    for a_file in assur_files:
        a_counter = 0
        file_num = re.search('_(\\d).txt', a_file)

        if file_num:
            num = file_num.group(1)

        with open(a_file, "r") as file_to_count:
            for line in file_to_count:
                a_counter += 1

        print(str(num) + " : " + str(a_counter))


try:
    main()
except Exception as mainExecuteError:
    print(mainExecuteError)
