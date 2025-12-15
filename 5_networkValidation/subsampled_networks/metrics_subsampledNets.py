# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import networkx as nx

from networkx.readwrite import json_graph
import json

import copy
import random
import csv

import seaborn as sns

import os
import re
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor
import glob
from functools import partial

from joypy import joyplot
import tempfile
from pathlib import Path

# === HELPER FUNCTIONS ===
def get_connected_graph(g):
    if nx.is_connected(g):
        return g
    else:
        largest_cc = max(nx.connected_components(g), key=len)
        return g.subgraph(largest_cc).copy()

def add_weights(net, weights):
  for _, row in weights.iterrows():
      v1, v2, w = row["OTU_1"], row["OTU_2"], row["weight"]
      if net.has_edge(v1, v2):
          net[v1][v2]["weight"] = w
          net[v1][v2]["invWeight"] = 1- w
      else:
          print(f"edge {v1},{v2} not found")


def create_graph(base_dir, prefix):

    match = re.match(r"([A-Za-z]+)_(\d+)", prefix)
    if match:
        MetObesity = match.group(1)
        number_network = int(match.group(2))
    else:
        raise ValueError(f"error extracting info: {prefix}")

    graphml_path = os.path.join(base_dir, f"{prefix}.graphml")

    G = nx.read_graphml(graphml_path)
    connectedG = get_connected_graph(G)

    connectedG.graph["MetObesity"] = MetObesity
    connectedG.graph["number_network"] = number_network
    connectedG.graph["number_ccs"] = nx.number_connected_components(G)

    mapping = {node: data['name'] for node, data in connectedG.nodes(data=True)}
    nx.relabel_nodes(connectedG, mapping, copy=False)

    weights_path = graphml_path.replace("_ig_mb.graphml", "_edge_weights.csv")

    if os.path.exists(weights_path):
            weights = pd.read_csv(weights_path)
            add_weights(connectedG, weights)
    else:
        print(f"#####  weights not found: {weights_path}")

    return connectedG

def add_abundance(net, metadata):
    abundance_dict = pd.Series(metadata.Abundance.values, index=metadata.Label).to_dict()
    nx.set_node_attributes(net, abundance_dict, 'Abundance')

def get_pos_edges(net):
  pos_edges_net, neg_edges_net = [], []
  for edge in net.edges(data = True):
    if edge[-1]['weight'] > 0:
      pos_edges_net.append(edge[-1]['weight'])

  return(100 * len(pos_edges_net)/net.size())


def compute_graph_properties(G, file_path):
    clsness = nx.closeness_centrality(G, distance="invWeight")
    bet = nx.betweenness_centrality(G, weight="invWeight")
    degree_sequence = sorted((d for n, d in G.degree()), reverse=True)

    props = {
        "file": file_path,
        "MetObesity": G.graph.get("MetObesity"),
        "number_network": G.graph.get("number_network"),
        "number_ccs" : G.graph.get("number_ccs"),
        "num_nodes": G.number_of_nodes(),
        "num_edges": G.number_of_edges(),
        "pos_edges": get_pos_edges(G),
        "density": nx.density(G),
        "average_clustering": nx.average_clustering(G, weight="invWeight"),
        "average_shortest_path_length": nx.average_shortest_path_length(G, weight="invWeight"),
        "mean_betweenness": sum(bet.values()) / len(bet),
        "mean_closeness": sum(clsness.values()) / len(clsness),
        "mean_degree": sum(degree_sequence) / len(degree_sequence),
        "max_degree": max(degree_sequence),
        "max_KCore": max((nx.core_number(G)).values())
    }
    return props


def process_network(prefix, base_dir):

    try:
        print(f"########## processing {prefix} ##########")
        G = create_graph(base_dir, prefix+'_ig_mb')
        props = compute_graph_properties(G, os.path.join(base_dir, prefix))
        return props
    except Exception as e:
        print(f"### error processing {prefix}: {e}")
        return None

# === GET FILES ===
folder_path = '.'
graphml_files = glob.glob(os.path.join(folder_path, '*.graphml'))
len(graphml_files) # bien.


all_files = os.listdir(folder_path)
prefixes = sorted({
    re.match(r"([A-Za-z]+_\d+)\_ig_mb.graphml", f).group(1)
    for f in all_files
    if f.endswith(".graphml") and re.match(r"([A-Za-z]+_\d+)\_ig_mb.graphml", f)
})

results = []
process_func = partial(process_network, base_dir=folder_path)

with ProcessPoolExecutor() as executor:  
    results = list(executor.map(process_func, prefixes))

results = [r for r in results if r is not None]


results_df = pd.DataFrame(results)
print(results_df.head())

results_df.to_csv('subsampledNetsData.csv', index=False)

# === GET KEYSTONE TAXA ===
def find_keystone_taxa(G, name_converter=None, quantile_cutoff=0.9):

    bet = nx.betweenness_centrality(G, weight="invWeight")
    deg = dict(G.degree())

    df = pd.DataFrame({
        'OTU': list(bet.keys()),
        'bet': list(bet.values()),
        'deg': [deg[n] for n in bet.keys()]
    })
    
    if name_converter is not None:
        df = df.merge(name_converter, on='OTU', how='left')

    # Cuantiles
    bet_q = df['bet'].quantile(quantile_cutoff)
    deg_q = df['deg'].quantile(quantile_cutoff)

    df['bottleneck'] = df['bet'] > bet_q
    df['hub'] = df['deg'] > deg_q
    df['keystone'] = df['bottleneck'] & df['hub']
    df['score'] = df[['bottleneck', 'hub', 'keystone']].sum(axis=1)

    return df


def keystone_network(prefix, base_dir, quantile_cutoff=0.9):
    try:
      print(f"########## PROCESSING {prefix} ##########")

      G = create_graph(base_dir, prefix+'_ig_mb')

      meta_path = os.path.join(base_dir, f"{prefix}_metadata.csv")
      name_converter = pd.read_csv(meta_path) if os.path.exists(meta_path) else None
      name_converter = name_converter.rename(columns={'Label': 'OTU'})
      df = find_keystone_taxa(G, name_converter=name_converter, quantile_cutoff=quantile_cutoff)

      df["Network"] = prefix
      df["Group"] = prefix.split("_")[0]

      return {prefix: df}

    except Exception as e:
        print(f"### error processing {prefix}: {e}")
        return None


# === Run (parallel execution) ===
def run_keystone_parallel(prefixes, base_dir, quantile_cutoff=0.9):
    process_func = partial(keystone_network, base_dir=base_dir, quantile_cutoff=quantile_cutoff)

    with ProcessPoolExecutor() as executor:
        results = list(executor.map(process_func, prefixes))

    results = [r for r in results if r is not None]

    keystone_dict = {k: v for d in results for k, v in d.items()}

    df_all = pd.concat(keystone_dict.values(), ignore_index=True)

    return keystone_dict, df_all



keystone_dict, df_keystone = run_keystone_parallel(prefixes, folder_path)


output_path = 'keystone_data_subsampledNets.csv'
df_keystone.to_csv(output_path, index=False)