# Sanath kahagalage
# UNSW Canberra
# 2021/03/12
#______________________________________________________________________________
from numpy import array, mean, abs, linalg
import pprint
p = pprint.pprint
from enum import IntEnum
import networkx as nx
from itertools import combinations
from collections import defaultdict
import numpy
from scipy.sparse import csr_matrix


#______________________________________________________________________________
def get_raw_particle_particle_data():
    filename = 'contactData/DTLZ7_8D.txt'
    with open(filename) as infile:
        #next(infile) # skip header
        particle_particle_data = infile.readlines()
    return particle_particle_data

#______________________________________________________________________________
def get_contact_network():
    contact_network = nx.Graph()

    for line in get_raw_particle_particle_data():
        tokens = line.split()
        pid1, pid2 = int(tokens[0]), int(tokens[1])
        contact_network.add_edge(pid1, pid2)

    return contact_network

#________________________________________________________________________________________
def get_contacts():
    contacts = dict()
    # These are particle-particle contacts
    for line in get_raw_particle_particle_data():
        tokens = line.split()
        pid1, pid2 = int(tokens[0]), int(tokens[1])
        dis = float(tokens[2])
        contacts[pid1, pid2] = (
            dis,
            
        )
        contacts[pid2, pid1] = (
            dis,
            
        )


    return contacts

 
#_______________________________________________ get minimum cut tree (Gomory-Hu tree)

def get_cut_tree():

    network = get_contact_network()
    contacts = get_contacts()
        
    # Construct a weighted flow network
    print('Constructing weighted flow network')
    flow_network = nx.Graph()
    for node in network.nodes():
        neighbours = network.neighbors(node)        
        neighbours_always = set(neighbours) 
        for neighbour in neighbours_always:
            euc_dis = contacts[node,neighbour][0]
            if euc_dis == 0.0:
                cap = 1/(1e-50)**2
            else:
                cap = 1/(euc_dis)**2
        
            edge_weight             = epl           


            flow_network.add_edge(node, neighbour, capacity= edge_weight)

    # Calculate max flow and corresponding minimum cut
    print('Starting calculation')
    T   = nx.gomory_hu_tree(flow_network, capacity='capacity', flow_func=None)

    return T

T = get_cut_tree()
#print T.nodes()
#print T.edges()
print (len(T.edges()))
print (len(T.nodes()))

def minimum_edge_weight_in_shortest_path(T, u, v):
    path = nx.shortest_path(T, u, v, weight="weight")
    return min((T[u][v]["weight"], (u, v)) for (u, v) in zip(path, path[1:]))

edgeAll = []
cap = []

for edge in T.edges():
    u = edge[0]
    v = edge[1]
    cut_value, edge = minimum_edge_weight_in_shortest_path(T, u, v)
    edgeAll.append(list(edge))
    cap.append(cut_value)

with open('results/gomory-hu-tree_DTLZ7_8D.txt', 'w') as outfile:
    for edge in edgeAll:
        print('%d %d' %(edge[0],edge[1]),file=outfile)


with open('results/capacity_DTLZ7_8D.txt', 'w') as outfile:
    numpy.savetxt(outfile, cap)
