import pandas as pd
import leidenalg as la
import igraph as ig
import re

'''
Takes in a igraph object and returns the result partition object
Utilizes best experimental settings, subject to change
'''
def generate_leiden_optimized(graph):
    partitions = la.find_partition(graph, la.RBConfigurationVertexPartition, resolution_parameter=0.8, n_iterations=100, weights=graph.es['weight'],
                               initial_membership=la.find_partition(graph, la.ModularityVertexPartition, n_iterations=100, weights=graph.es['weight']).membership)
    return partitions

'''
The following function takes in a list of node names(community) and selects a adjacency matrix subset
'''
def make_adj_matrices_from_community(partition, master_adj):
    node_names = master_adj.columns
    # make adjacency matrices out of the communities
    communities = [[] for i in range(partition._len)]
    for idx, elem in enumerate(partition.membership):
        communities[elem].append(node_names[idx])
    adj_communities = []
    for idx in range(len(communities)):
        adj_communities.append(master_adj.loc[communities[idx]][communities[idx]])
    return adj_communities
'''
Breaks chrX:A-B into {"chr":"x","start":"A","stop":"B"}
'''
def breakup_name_val(name):
    parts = re.split('[:-]', name)
    return {"chr":parts[0],"start":parts[1],"stop":parts[2]}

'''
The following function takes a list of node names(community) and generates a string of a track file
also_swap_loci duplicate each track line with loci2 loci1 along with the standard loci1 loci2
'''
def convert_adj_mat_to_track_file(adj_mat, also_swap_loci = True):
    final_string = ""
    # Iterate through every row
    for idx, row in adj_mat.iterrows():
        locus1 = breakup_name_val(idx)
        # Iterate through interactions
        for name in adj_mat.columns:
            score = row[name]
            if score > 0:
                locus2 = breakup_name_val(name)
                final_string += locus1["chr"]+"\t"+locus1["start"]+"\t"+locus1["stop"]+"\t"+name+","+str(score)+"\n"
                # Might need a second one according to epigenome browser documentation but give user the option
                # True by default
                if also_swap_loci:
                    final_string += locus2["chr"]+"\t"+locus2["start"]+"\t"+locus2["stop"]+"\t"+idx+","+str(score)+"\n"
    return final_string  

'''
Take in a list of adjacency matrices for communities and save them as longrange .bedgraph files
The adjacency matrix weight for each contact is used for the score parameter
format in final file:
chrX	A	B	chrX:C-D,WEIGHT
chrX	C	D	chrX:A-B,WEIGHT
'''
def save_tracks_for_communities(adj_mats, path):
    for idx in range(len(adj_mats)):
        track_file = convert_adj_mat_to_track_file(adj_mats[idx], True)
        text_file = open(path + "/" + "community_"+str(idx)+".bedgraph", "w")
        text_file.write(track_file)
        text_file.close()

'''
Example run of the code using a hypothetical sum.csv file
'''
def example_code():
    sums = pd.read_csv("sum.csv", index_col="name")
    node_names = sums.index.values.tolist()
    # Get the values as np.array, it's more convenenient.
    sum_vals = sums.values
    # Create graph, A.astype(bool).tolist() or (A / A).tolist() can also be used.
    g = ig.Graph.Adjacency((sum_vals > 0).tolist())
    # Add edge weights and node labels.
    g.es['weight'] = sum_vals[sum_vals.nonzero()]
    g.vs['label'] = node_names  # or a.index/a.columns
    best_partition = generate_leiden_optimized(g)
    adj_mats = make_adj_matrices_from_community(best_partition, sums)
    save_tracks_for_communities(adj_mats, "./output_folder")
