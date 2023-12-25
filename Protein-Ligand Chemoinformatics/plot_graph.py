import pandas as pd
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt

import sys

def main():
	"""
	Function is called when the file is called in the terminal
	"""
	# check that the file is being properly used
	if (len(sys.argv) !=4):
		print("Please specify an input file and an output file as args.")
		return

	# Get variables from shell
	edge_file = sys.argv[1]
	nodes_file = sys.argv[2]
	output_file = sys.argv[3]

	# Read in dataframes
	edges_df = pd.read_csv(edge_file, delimiter=' ', names = ["from", "to"])
	nodes_df = pd.read_csv(nodes_file, delimiter=',')

	# Get dictionary mapping 'uniprot_accession' => 'uniprot_id'
	acc_to_id = nodes_df.set_index('uniprot_accession')['uniprot_id'].to_dict()

	# Replace all 'uniprot_accession' by 'uniprot_id' in edges dataframe
	edges_df = edges_df.replace(acc_to_id)

	# Define dictionary mapping 'indications' => colors
	colors_dict = {"bp":"red", "bp;cholesterol":"green", "bp;cholesterol;diabetes":"blue","bp;diabetes":"purple"}

	# Replace 'indications' with 'colors' from colors_dict and drop 'uniprot_accession' column
	nodes_df = nodes_df.replace({"indications": colors_dict})
	nodes_df = nodes_df.drop('uniprot_accession', axis=1)

	# Reference for plotting: https://www.python-graph-gallery.com/324-map-a-color-to-network-nodes
	G = nx.from_pandas_edgelist(edges_df, 'from', 'to', create_using=nx.Graph())
	nodes_df= nodes_df.set_index('uniprot_id')
	nodes_df=nodes_df.reindex(G.nodes())

	# Customize the nodes and draw
	plt.figure(figsize=(8,8), dpi=150) 
	nx.draw_shell(G, with_labels=True, node_color=nodes_df['indications'], node_size=150,font_size=8, font_color="#000000",edge_color="grey")
	# plt.show()
	plt.savefig(output_file)

if __name__=="__main__":
    main()
