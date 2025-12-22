# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import networkx as nx

from networkx.readwrite import json_graph
import json

import pandas

import copy
import random
import csv

import seaborn as sns


# === LOAD NETWORKS ===
MHNO = nx.read_graphml("mhno.graphml")
MHO = nx.read_graphml("mho.graphml")
MUNO = nx.read_graphml("muno.graphml")
MUO = nx.read_graphml("muo.graphml")

for net in (MHO, MHNO, MUO, MUNO):
  print(net.order(), net.size(),
        nx.number_connected_components(net))

netdict = {
    "MHNO": MHNO,
    "MHO": MHO,
    "MUNO": MUNO,
    "MUO": MUO
}

# fix node name mapping
for net in netdict.keys():
  mapping = {node: data['name'] for node, data in netdict[net].nodes(data=True)}
  nx.relabel_nodes(netdict[net], mapping, copy=False)


def get_connected_graph(net):
  # Get a list of connected components
  components = list(nx.connected_components(net))
  # Find the largest connected component
  largest_component = max(components, key=len)
  # Create a new graph containing only the largest component
  net_largest = net.subgraph(largest_component)
  print('Connected graph:', nx.is_connected(net_largest))
  return(net_largest)


MHO = get_connected_graph(MHO)
MHNO = get_connected_graph(MHNO)
MUO = get_connected_graph(MUO)
MUNO = get_connected_graph(MUNO)

# add abundances -------------------------------
# MHO
metadata_MHO  = pd.read_csv("mho_metadata.csv")
metadata_MHO.rename(columns=lambda x: x.replace('.', '|'), inplace=True)

# MHNO
metadata_MHNO  = pd.read_csv("mhno_metadata.csv")
metadata_MHNO.rename(columns=lambda x: x.replace('.', '|'), inplace=True)

# MUO
metadata_MUO  = pd.read_csv("muo_metadata.csv")
metadata_MUO.rename(columns=lambda x: x.replace('.', '|'), inplace=True)

# MUNO
metadata_MUNO  = pd.read_csv("muno_metadata.csv")
metadata_MUNO.rename(columns=lambda x: x.replace('.', '|'), inplace=True)

def add_abundance(net, metadata):
  abundance_dict = pd.Series(metadata.Abundance.values, index=metadata.Label).to_dict()
  nx.set_node_attributes(net, abundance_dict, 'Abundance')

add_abundance(MHO, metadata_MHO)
add_abundance(MHNO, metadata_MHNO)
add_abundance(MUO, metadata_MUO)
add_abundance(MUNO, metadata_MUNO)


# add edge weights -----------------------------------
# MHNO
weights_MHNO  = pd.read_csv("edgeweights_mhno.csv")
weights_MHNO.rename(columns=lambda x: x.replace('.', '|'), inplace=True)

# MHO
weights_MHO  = pd.read_csv("/edgeweights_mho.csv")
weights_MHO.rename(columns=lambda x: x.replace('.', '|'), inplace=True)

# MUNO
weights_MUNO  = pd.read_csv("/edgeweights_muno.csv")
weights_MUNO.rename(columns=lambda x: x.replace('.', '|'), inplace=True)

# MUO
weights_MUO  = pd.read_csv("/edgeweights_muo.csv")
weights_MUO.rename(columns=lambda x: x.replace('.', '|'), inplace=True)

def add_weights(net, weights):
  for _, row in weights.iterrows():
      v1, v2, w = row["v1"], row["v2"], row["edgeweight"]
      if net.has_edge(v1, v2):
          net[v1][v2]["weight"] = w
          net[v1][v2]["invWeight"] = 1 - w
      else:
          print(f"edge {v1},{v2} not found")

add_weights(MHNO, weights_MHNO)
add_weights(MHO, weights_MHO)
add_weights(MUNO, weights_MUNO)
add_weights(MUO, weights_MUO)

# === Calculate global properties ===
def get_pos_neg_edges(net):
  pos_edges_net, neg_edges_net = [], []
  for edge in net.edges(data = True):
    if edge[-1]['weight'] > 0:
      pos_edges_net.append(edge[-1]['weight'])
    else:
      neg_edges_net.append(edge[-1]['weight'])
  print('% positive edge: {0}\n% negative edges: {1}'.
      format(100 * len(pos_edges_net)/net.size(),
             100 * len(neg_edges_net)/net.size()))

print("######### MHO #########")
get_pos_neg_edges(MHO)
print("######### MHNO #########")
get_pos_neg_edges(MHNO)
print("######### MUO #########")
get_pos_neg_edges(MUO)
print("######### MUNO #########")
get_pos_neg_edges(MUNO)

# === Get Local Properties ===
# Degree
def get_basic_properties(net):
  degree_sequence = sorted((d for n, d in net.degree()), reverse=True)
  print('Number of nodes: {0} \nNumber of edges: {1} \nMax degree: {2} \nDensity: {3}'
  .format(net.order(), net.size(), max(degree_sequence), nx.density(net)))

print("MHO #############"); get_basic_properties(MHO)
print("MHNO ############"); get_basic_properties(MHNO)
print("MUO #############"); get_basic_properties(MUO)
print("MUNO ############"); get_basic_properties(MUNO)

MHO_names = list(MHO.nodes())
MHNO_names = list(MHNO.nodes())
MUO_names = list(MUO.nodes())
MUNO_names = list(MUNO.nodes())

meandegMHO = sum([i[1] for i in MHO.degree])/MHO.order()
meandegMHNO = sum([i[1] for i in MHNO.degree])/MHNO.order()
meandegMUO = sum([i[1] for i in MUO.degree])/MUO.order()
meandegMUNO = sum([i[1] for i in MUNO.degree])/MUNO.order()
print(meandegMHO, meandegMHNO, meandegMUO, meandegMUNO)


# Betweenness
def get_betw(G):
  dc = nx.degree_centrality(G)
  clsness = nx.closeness_centrality(G, distance = "invWeight")
  btwness_G = nx.betweenness_centrality(G, weight="invWeight")

  print('Centrality (mean):', sum(dc.values())/len(dc))
  print('Closeness (mean):', sum(clsness.values())/len(clsness))
  print('Betweenness (mean):', sum(btwness_G.values())/len(btwness_G))
  return(dc, clsness, btwness_G)

print("MHNO ###########################")
netdict["MHNO"]["centrality"], netdict["MHNO"]["closeness"], netdict["MHNO"]["betweenness"] = get_betw(MHNO)
print("MHO ###########################")
netdict["MHO"]["centrality"], netdict["MHO"]["closeness"], netdict["MHO"]["betweenness"] = get_betw(MHO)
print("MUNO ###########################")
netdict["MUNO"]["centrality"], netdict["MUNO"]["closeness"], netdict["MUNO"]["betweenness"] = get_betw(MUNO)
print("MUO ###########################")
netdict["MUO"]["centrality"], netdict["MUO"]["closeness"], netdict["MUO"]["betweenness"] = get_betw(MUO)

names = ["MHO", "MHNO", "MUO", "MUNO"]
for (i,net) in enumerate([MHO, MHNO, MUO, MUNO]):
  print('#### ', names[i], '###################')
  print('Average shortest path length:', nx.average_shortest_path_length(net, weight = "invWeight"))

# === Export Dataframes ===

netdict = {
    "MHNO": MHNO,
    "MHO": MHO,
    "MUNO": MUNO,
    "MUO": MUO
}

all_nodes = set().union(*[G.nodes() for G in netdict.values()])
print("Total number of nodes: ", len(all_nodes))

metrics = {}

for net, G in netdict.items():
    metrics[net] = {
        "degree": dict(G.degree()),
        "betweenness": nx.betweenness_centrality(G, weight = "invWeight"),
        "closeness": nx.closeness_centrality(G)
        }

data = []

for node in all_nodes:
    for net in netdict.keys():
        degree = metrics[net]["grado"].get(node, None)
        betweenness = metrics[net]["betweenness"].get(node, None)
        closeness = metrics[net]["closeness"].get(node, None)

        data.append({
            "nodo": node,
            "red": net,
            "grado": degree,
            "betweenness": betweenness,
            "closeness": closeness
        })

df = pd.DataFrame(data)
df.to_csv("./network_metrics.csv")


shortest_path_lengths = {}

for net, G in netdict.items():
    path_lengths = []
    for nodo_origen, lengths in dict(nx.shortest_path_length(G, weight = "invWeight")).items():
        path_lengths.extend(lengths.values()) 

    shortest_path_lengths[net] = path_lengths

with open("./shortest_paths.json", "w") as f:
    json.dump(shortest_path_lengths, f)

# === K-Core distributions ===

cores_MHO = nx.core_number(MHO)
cores_MHNO = nx.core_number(MHNO)
cores_MUO = nx.core_number(MUO)
cores_MUNO = nx.core_number(MUNO)

from collections import Counter
mho_col = "#2a9d8f"
mhno_col = "#264653"
muo_col = "#edafb8"
muno_col = "#703d67"

fig, axes = plt.subplots(1, 4, figsize=(8, 2.5), dpi=1200)
xlim = (0, 6)

freqs = Counter(cores_MHO.values())
vals, counts = zip(*sorted(freqs.items(), key=lambda x: x[1], reverse=True))
axes[0].bar(vals, counts, color=mho_col)
axes[0].set_title("MHO", fontsize=12)
axes[0].set_ylabel("Frequency", fontsize=12)
axes[0].set_xlabel("K-Cores", fontsize=12)
axes[0].set_xticks(range(6))
axes[0].set_xlim(xlim)
axes[0].tick_params(axis="both", which="major", labelsize=10)
axes[0].set_title("MHO", fontsize=12, fontweight='bold')

freqs = Counter(cores_MHNO.values())
vals, counts = zip(*sorted(freqs.items(), key=lambda x: x[1], reverse=True))
axes[1].bar(vals, counts, color=mhno_col)
axes[1].set_title("MHNO", fontsize=12)
axes[1].set_xlabel("K-Cores", fontsize=12)
axes[1].set_xticks(range(6))
axes[1].set_xlim(xlim)
axes[1].tick_params(axis="both", which="major", labelsize=10)
axes[1].set_title("MHNO", fontsize=12, fontweight='bold')

freqs = Counter(cores_MUO.values())
vals, counts = zip(*sorted(freqs.items(), key=lambda x: x[1], reverse=True))
axes[2].bar(vals, counts, color=muo_col)
axes[2].set_title("MUO", fontsize=12)
axes[2].set_xlabel("K-Cores", fontsize=12)
axes[2].set_xticks(range(6))
axes[2].set_xlim(xlim)
axes[2].tick_params(axis="both", which="major", labelsize=10)
axes[2].set_title("MUO", fontsize=12, fontweight='bold')

freqs = Counter(cores_MUNO.values())
vals, counts = zip(*sorted(freqs.items(), key=lambda x: x[1], reverse=True))
axes[3].bar(vals, counts, color=muno_col)
axes[3].set_title("MUNO", fontsize=12)
axes[3].set_xlabel("K-Cores", fontsize=12)
axes[3].set_xticks(range(6))
axes[3].set_xlim(xlim)
axes[3].tick_params(axis="both", which="major", labelsize=10)
axes[3].set_title("MUNO", fontsize=12, fontweight='bold')

plt.tight_layout()

plt.savefig("./kcores.png",
            format='png', dpi=1200,
            bbox_inches='tight')



df = pd.DataFrame(list(cores_MHO.items()), columns=['Node', 'k-core'])
df.to_csv('./MHO_kcores.csv', index=False)

df = pd.DataFrame(list(cores_MHNO.items()), columns=['Node', 'k-core'])
df.to_csv('./MHNO_kcores.csv', index=False)

df = pd.DataFrame(list(cores_MUO.items()), columns=['Node', 'k-core'])
df.to_csv('./MUO_kcores.csv', index=False)

df = pd.DataFrame(list(cores_MUNO.items()), columns=['Node', 'k-core'])
df.to_csv('./MUNO_kcores.csv', index=False)

# === Keystone Taxa Calculation ===

def find_keystone_taxa(G, name_converter=None, network_name=None, quantile_cutoff=0.9):

    bet = nx.betweenness_centrality(G, weight = "invWeight")
    deg = dict(G.degree())

    df = pd.DataFrame({
        'OTU': list(bet.keys()),
        'bet': list(bet.values()),
        'deg': [deg[n] for n in bet.keys()]
    })

    if name_converter is not None:
        df = df.merge(name_converter, on='OTU', how='left')

    bet_q = df['bet'].quantile(quantile_cutoff)
    deg_q = df['deg'].quantile(quantile_cutoff)

    df['bottleneck'] = df['bet'] > bet_q
    df['hub'] = df['deg'] > deg_q
    df['keystone'] = df['bottleneck'] & df['hub']

    if network_name is not None:
        df['Network'] = network_name

    return df, bet_q, deg_q


df_MHNO, _, _ = find_keystone_taxa(MHNO, metadata_MHNO.rename(columns={'Label': 'OTU'}), network_name='MHNO')
df_MHO,  _, _ = find_keystone_taxa(MHO,  metadata_MHO.rename(columns={'Label': 'OTU'}), network_name='MHO')
df_MUNO, _, _ = find_keystone_taxa(MUNO, metadata_MUNO.rename(columns={'Label': 'OTU'}), network_name='MUNO')
df_MUO,  _, _ = find_keystone_taxa(MUO,  metadata_MUO.rename(columns={'Label': 'OTU'}), network_name='MUO')


df_all = pd.concat([df_MHNO, df_MHO, df_MUNO, df_MUO], ignore_index=True)


df_all.to_csv("./keystoneTaxa.csv",
              index = False)

df_all.head()