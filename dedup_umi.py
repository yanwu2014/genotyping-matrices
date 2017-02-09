'''
DedupUMI.py - Deduplicate reads that are coded with a UMI
=========================================================

:Author: Ian Sudbery
:Release: $Id$
:Date: |today|
:Tags: Python

'''

### Slightly modified version of Ian Sudbery's code ###

import sys
import pysam
import random
import collections
import itertools
import pandas as pd
import numpy as np
import cPickle as cp
from editdistance import eval
from distance import hamming
from collections import defaultdict

def breadth_first_search(node, adj_list):
    searched = set()
    found = set()
    queue = set()
    queue.update((node,))
    found.update((node,))

    while len(queue) > 0:
        node = (list(queue))[0]
        found.update(adj_list[node])
        queue.update(adj_list[node])
        searched.update((node,))
        queue.difference_update(searched)

    return found

def edit_distance(a,b):
    return eval(a,b)

def get_umi(read):
    return read.qname.split("_")[-1]

def get_average_umi_distance(umis):
    if len(umis) == 1:
        return -1
    dists = [edit_distance(*pair) for pair in itertools.combinations(umis, 2)]
    return float(sum(dists))/(len(dists))


def remove_umis(adj_list, cluster, nodes):
    '''removes the specified nodes from the cluster and returns
    the remaining nodes '''

    # list incomprehension: for x in nodes: for node in adj_list[x]: yield node
    nodes_to_remove = set([node
                           for x in nodes
                           for node in adj_list[x]] + nodes)

    return cluster - nodes_to_remove


class ClusterAndReducer:
    '''A functor that clusters a bundle of reads,
    indentifies the parent UMIs and returns the selected reads, umis and counts

    The initiation of the functor defines the methods:

      ** get_adj_list ** - returns the edges connecting the UMIs

      ** connected_components ** - returns clusters of connected components
                                   using the edges in the adjacency list

      ** get_best ** - returns the parent UMI(s) in the connected_components

      ** reduce_clusters ** - loops through the connected components in a
                              cluster and returns the unique reads. Optionally
                              returns lists of umis and counts per umi also

    Note: The get_adj_list and connected_components methods are not required by
    all custering methods. Where there are not required, the methods return
    None or the input parameters.

    '''

    ######## "get_best" methods ##########

    def _get_best_min_account(self, cluster, adj_list, counts):
        ''' return the min UMI(s) need to account for cluster'''
        if len(cluster) == 1:
            return list(cluster)

        sorted_nodes = sorted(cluster, key=lambda x: counts[x],
                              reverse=True)

        for i in range(len(sorted_nodes) - 1):
            if len(remove_umis(adj_list, cluster, sorted_nodes[:i+1])) == 0:
                return sorted_nodes[:i+1]

    def _get_best_higher_counts(self, cluster, counts):
        ''' return the UMI with the highest counts'''
        if len(cluster) == 1:
            return list(cluster)[0]
        else:
            sorted_nodes = sorted(cluster, key=lambda x: counts[x],
                                  reverse=True)
            return sorted_nodes[0]

    def _get_best_percentile(self, cluster, counts):
        ''' return all UMIs with counts >1% of the
        median counts in the cluster '''

        if len(cluster) == 1:
            return list(cluster)
        else:
            threshold = np.median(counts.values())/100
            return [read for read in cluster if counts[read] > threshold]

    def _get_best_null(self, cluster, counts):
        ''' return all UMIs in the cluster'''

        return list(cluster)


    ######## "get_adj_list" methods ##########

    def _get_adj_list_adjacency(self, umis, counts, threshold):
        ''' identify all umis within hamming distance threshold'''

        return {umi: [umi2 for umi2 in umis if
                      edit_distance(umi, umi2) <= threshold]
                for umi in umis}

    def _get_adj_list_directional_adjacency(self, umis, counts, threshold,
                                            use_hamming=False, countRatio = 1.5):
        ''' identify all umis within the hamming distance threshold
        and where the counts of the first umi is > (1.5 * second umi counts)-1'''
        if use_hamming:
            return {umi: [umi2 for umi2 in umis if 
                          hamming(umi, umi2) <= threshold 
                          and counts[umi] >= (counts[umi2]*countRatio)-1] 
                          for umi in umis}
        else:
            return {umi: [umi2 for umi2 in umis if 
                          edit_distance(umi, umi2) <= threshold 
                          and counts[umi] >= (counts[umi2]*countRatio)-1] 
                          for umi in umis}


    ######## "get_connected_components" methods ##########

    def _get_connected_components_adjacency(self, umis, graph, counts):
        ''' find the connected UMIs within an adjacency dictionary'''

        found = list()
        components = list()

        for node in sorted(graph, key=lambda x: counts[x], reverse=True):
            if node not in found:
                component = breadth_first_search(node, graph)
                found.extend(component)
                components.append(component)

        return components


    ######## "reduce_clusters" methods ##########

    def _reduce_clusters_multiple(self, clusters,
                                  adj_list, counts, stats=False):
        ''' collapse clusters down to the UMI(s) which account for the cluster
        using the adjacency dictionary and return the list of final UMIs'''

        # TS - the "adjacency" variant of this function requires an adjacency
        # list to identify the best umi, whereas the other variants don't
        # As temporary solution, pass adj_list to all variants

        final_umis = []
        umi_counts = []

        for cluster in clusters:
            parent_umis = self.get_best(cluster, adj_list, counts)

            if stats:
                final_umis.extend(parent_umis)
                umi_counts.extend([counts[umi] for umi in parent_umis])

        return final_umis, umi_counts

    def _reduce_clusters_single(self, clusters,
                                adj_list, counts, stats=True):
        ''' collapse clusters down to the UMI which accounts for the cluster
        using the adjacency dictionary and return the list of final UMIs'''

        final_umis = {}
        umi_counts = {}

        for cluster in clusters:
            parent_umi = self.get_best(cluster, counts)

            if stats:
                final_umis[parent_umi] = cluster 
                umi_counts[parent_umi] = sum([counts[x] for x in cluster])

        return final_umis, umi_counts

    
    
    def __init__(self, cluster_method="directional-adjacency"):
        ''' select the required class methods for the cluster_method'''

        if cluster_method == "adjacency":
            self.get_adj_list = self._get_adj_list_adjacency
            self.get_connected_components = self._get_connected_components_adjacency
            self.get_best = self._get_best_min_account
            self.reduce_clusters = self._reduce_clusters_multiple

        elif cluster_method == "directional-adjacency":
            self.get_adj_list = self._get_adj_list_directional_adjacency
            self.get_connected_components = self._get_connected_components_adjacency
            self.get_best = self._get_best_higher_counts
            self.reduce_clusters = self._reduce_clusters_single

        elif cluster_method == "cluster":
            self.get_adj_list = self._get_adj_list_adjacency
            self.get_connected_components = self._get_connected_components_adjacency
            self.get_best = self._get_best_higher_counts
            self.reduce_clusters = self._reduce_clusters_single

        elif cluster_method == "percentile":
            self.get_adj_list = self._get_adj_list_null
            self.get_connected_components = self._get_connected_components_null
            self.get_best = self._get_best_percentile
            self.reduce_clusters = self._reduce_clusters_no_network

        if cluster_method == "unique":
            self.get_adj_list = self._get_adj_list_null
            self.get_connected_components = self._get_connected_components_null
            self.get_best = self._get_best_null
            self.reduce_clusters = self._reduce_clusters_no_network

    

def main(argv=None):
    ClusterAndReducer()


if __name__ == "__main__": main()
