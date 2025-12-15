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

# === LOAD NETWORKS ===
MHNO = nx.read_graphml("mhno.graphml")
MHO = nx.read_graphml("mho.graphml")
MUNO = nx.read_graphml("muno.graphml")
MUO = nx.read_graphml("muo.graphml")

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

# add metadata:
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

# add edge weights:
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

# === Add weights & edit networks ===
def add_weights(net, weights):
  for _, row in weights.iterrows():
      v1, v2, w = row["v1"], row["v2"], row["edgeweight"]
      if net.has_edge(v1, v2):
          net[v1][v2]["weight"] = w
          net[v1][v2]["invWeight"] = 1- w
      else:
          print(f"edge {v1},{v2} not found")

add_weights(MHNO, weights_MHNO)
add_weights(MHO, weights_MHO)
add_weights(MUNO, weights_MUNO)
add_weights(MUO, weights_MUO)

# keep only positive edges:
def filter_positive_edges(net):
    
    G_filtered = net.__class__()

    for n, attrs in net.nodes(data=True):
        G_filtered.add_node(n, **attrs)

    for u, v, data in net.edges(data=True):
        w = data.get("weight", None)
        if w is not None and w > 0:
            G_filtered.add_edge(u, v, **data)

    return G_filtered

MHNO = filter_positive_edges(MHNO)
MHO = filter_positive_edges(MHO)
MUNO = filter_positive_edges(MUNO)
MUO = filter_positive_edges(MUO)

# Keep only largest connected component
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

# Include node abundances

def add_abundance(net, metadata):
  # Crear un diccionario donde la clave sea el nombre del nodo y el valor sea el valor de 'Abundance'
  abundance_dict = pd.Series(metadata.Abundance.values, index=metadata.Label).to_dict()
  # Asignar el atributo 'Abundance' a los nodos
  nx.set_node_attributes(net, abundance_dict, 'Abundance')

add_abundance(MHO, metadata_MHO)
add_abundance(MHNO, metadata_MHNO)
add_abundance(MUO, metadata_MUO)
add_abundance(MUNO, metadata_MUNO)

for n, d in list(MHNO.nodes(data=True))[:5]:
            print(n, d)

# === NETWORK ATTACK FRAMEWORK ===

mho_col = "#2a9d8f"
mhno_col = "#264653"
muo_col = "#edafb8"
muno_col = "#703d67"

x_vals_MHO = [ 100 *i/MHO.order() for i in range(MHO.order())]
x_vals_MHNO = [ 100 *i/MHNO.order() for i in range(MHNO.order())]
x_vals_MUO = [ 100 *i/MUO.order() for i in range(MUO.order())]
x_vals_MUNO = [ 100 *i/MUNO.order() for i in range(MUNO.order())]

## Helper functions -------------------------------
### Calculate percolation threshold & NR50

def calculate_critical_threshold(x_vals_fraction, lcc_sizes_normalized_fraction):

    if len(lcc_sizes_normalized_fraction) <= 1:
        return 0.0 # No se puede calcular la derivada con menos de 2 puntos

    # Calcular la diferencia (aproximación de la pendiente)
    slopes = np.diff(lcc_sizes_normalized_fraction)

    # Encontrar el índice donde la pendiente es más negativa (mayor caída)
    # Excluir el último punto de slopes si corresponde al final del ataque (LCC=0)
    # y queremos el punto de inflexión ANTES de que se vacíe.
    # Si la curva es plana o casi plana, argmin podría dar un valor cercano a 0 o 1.
    if len(slopes) == 0:
        return 0.0

    critical_index = np.argmin(slopes)

    # El umbral crítico es la fracción de nodos removidos en este índice.
    # `x_vals_fraction` tiene un punto más que `slopes`,
    # por lo que el `slopes[i]` corresponde al cambio entre `x_vals_fraction[i]` y `x_vals_fraction[i+1]`.
    # Por lo tanto, el punto de la caída más pronunciada se considera en `x_vals_fraction[critical_index + 1]`.
    f_c = x_vals_fraction[critical_index + 1]

    return f_c

# nr50:
def find_i50(x_vals_fraction, lcc_sizes_normalized_fraction):
  i50 = []
  for i,j in zip(x_vals_fraction, lcc_sizes_normalized_fraction):
    if j>=0.50:
      i50.append(i)
  return(i50[-1])

### Perform targeted attacks:
def attack_degree(g_original, n):

  g = g_original.copy() # Work on a copy to preserve the original graph
  lcc_sizes = []

  # Initial LCC size before any removal
  try:
    initial_lcc_size = len(max(nx.connected_components(g), key=len))
  except ValueError: # Handle empty graph or no connected components
    initial_lcc_size = 0

  lcc_sizes.append(initial_lcc_size) # Add initial state

  for i in range(n):
    if not g.nodes(): # If graph becomes empty
        lcc_sizes.append(0)
        continue

    # Get node with max degree; handle cases where degree might be zero for all remaining nodes
    try:
        node = max(g.degree, key=lambda x: x[1])[0]
    except ValueError: # No nodes left with non-zero degree (e.g., isolated nodes)
        lcc_sizes.append(0)
        g.clear() # Clear the graph to avoid re-evaluating
        continue

    g.remove_node(node)

    try:
      # Find LCC in the remaining graph
      if g.nodes(): # Check if graph still has nodes
        lcc_size = len(max(nx.connected_components(g), key=len))
      else:
        lcc_size = 0 # No nodes left
    except ValueError: # No connected components (e.g., isolated nodes only)
      lcc_size = 0

    lcc_sizes.append(lcc_size)

  return np.array(lcc_sizes)


def attack_betweenness(g_original, n):
 
  g = g_original.copy() 
  lcc_sizes = []

  try:
    initial_lcc_size = len(max(nx.connected_components(g), key=len))
  except ValueError: 
    initial_lcc_size = 0

  lcc_sizes.append(initial_lcc_size) 

  for i in range(n):
    if not g.nodes(): 
        lcc_sizes.append(0)
        continue

    betweenness_dict = nx.betweenness_centrality(g, weight = "invWeight")

    try:
        node_to_remove = max(betweenness_dict, key=betweenness_dict.get)
    except ValueError: 
        lcc_sizes.append(0)
        g.clear()
        continue

    if betweenness_dict[node_to_remove] == 0 and g.order() > 0:
        pass


    g.remove_node(node_to_remove)

    try:
      if g.nodes(): 
        lcc_size = len(max(nx.connected_components(g), key=len))
      else:
        lcc_size = 0 
    except ValueError:
      lcc_size = 0

    lcc_sizes.append(lcc_size)

  return np.array(lcc_sizes)


def attack_abundance(g_original, n):
  g = g_original.copy() 
  lcc_sizes = []


  try:
    initial_lcc_size = len(max(nx.connected_components(g), key=len))
  except ValueError: 
    initial_lcc_size = 0

  lcc_sizes.append(initial_lcc_size) 

  for i in range(n):
    if not g.nodes(): 
        lcc_sizes.append(0)
        continue

    nodes_with_abundance = [node for node in g.nodes() if 'Abundance' in g.nodes[node]]

    if not nodes_with_abundance: 
        lcc_sizes.append(0)
        g.clear() 
        continue

    node_to_remove = sorted(nodes_with_abundance, key=lambda node: g.nodes[node]['Abundance'])[-1]

    g.remove_node(node_to_remove)

    try:
      if g.nodes(): 
        lcc_size = len(max(nx.connected_components(g), key=len))
      else:
        lcc_size = 0 
    except ValueError:
      lcc_size = 0

    lcc_sizes.append(lcc_size)

  return np.array(lcc_sizes)

def attack_inverse_abundance(g_original, n):
  g = g_original.copy() 
  lcc_sizes = []

  try:
    initial_lcc_size = len(max(nx.connected_components(g), key=len))
  except ValueError: 
    initial_lcc_size = 0

  lcc_sizes.append(initial_lcc_size) 
  for i in range(n):
    if not g.nodes():
        lcc_sizes.append(0)
        continue

    nodes_with_abundance = [node for node in g.nodes() if 'Abundance' in g.nodes[node]]

    if not nodes_with_abundance:
        lcc_sizes.append(0)
        g.clear()
        continue

    node_to_remove = sorted(nodes_with_abundance, key=lambda node: g.nodes[node]['Abundance'])[0]

    g.remove_node(node_to_remove)

    try:
      if g.nodes():
        lcc_size = len(max(nx.connected_components(g), key=len))
      else:
        lcc_size = 0
    except ValueError:
      lcc_size = 0

    lcc_sizes.append(lcc_size)

  return np.array(lcc_sizes)


fc_values_inv_abundance = []

# === Attack simulation function ===



def simulate_attack(g_original, attack_type):
  
  g = g_original.copy() 
  n = g.order()

  if attack_type == 'degree':
    lcc_sizes = attack_degree(g, n)
  elif attack_type == 'betweenness':
    lcc_sizes = attack_betweenness(g, n)
  elif attack_type == 'abundance':
    lcc_sizes = attack_abundance(g, n)
  elif attack_type == 'inverse_abundance':
    lcc_sizes = attack_inverse_abundance(g, n)
  else:
    raise ValueError(f"Unknown attack type: {attack_type}")

  initial_lcc_size = lcc_sizes[0] if lcc_sizes[0] > 0 else 1

  lcc_sizes_normalized_fraction = lcc_sizes / initial_lcc_size

  x_vals_fraction = np.array([i / n for i in range(n + 1)])

  # get percolation threshold & nr50 values
  fc_value = calculate_critical_threshold(x_vals_fraction, lcc_sizes_normalized_fraction)
  i50_value = find_i50(x_vals_fraction, lcc_sizes_normalized_fraction)


  return {
      'x_values': x_vals_fraction,
      'y_values': lcc_sizes_normalized_fraction,
      'fc_value': fc_value,
      'i50_value': i50_value
  }


# Attack by Degree
print("Attacks by Degree -------")
print("\tMHO")
mho_degree_results = simulate_attack(MHO, 'degree')
print("\tMHNO")
mhno_degree_results = simulate_attack(MHNO, 'degree')
print("\tMUO")
muo_degree_results = simulate_attack(MUO, 'degree')
print("\tMUNO")
muno_degree_results = simulate_attack(MUNO, 'degree')

# Attack by Betweenness
print("Attacks by Betweenness -------")
print("\tMHO")
mho_betweenness_results = simulate_attack(MHO, 'betweenness')
print("\tMHNO")
mhno_betweenness_results = simulate_attack(MHNO, 'betweenness')
print("\tMUO")
muo_betweenness_results = simulate_attack(MUO, 'betweenness')
print("\tMUNO")
muno_betweenness_results = simulate_attack(MUNO, 'betweenness')

# Attack by Abundance
print("Attacks by Descending Abundance -------")
print("\tMHO")
mho_abundance_results = simulate_attack(MHO, 'abundance')
print("\tMHNO")
mhno_abundance_results = simulate_attack(MHNO, 'abundance')
print("\tMUO")
muo_abundance_results = simulate_attack(MUO, 'abundance')
print("\tMUNO")
muno_abundance_results = simulate_attack(MUNO, 'abundance')

# Attack by Inverse Abundance
print("Attacks by Increasing Abundance -------")
print("\tMHO")
mho_inv_abundance_results = simulate_attack(MHO, 'inverse_abundance')
print("\tMHNO")
mhno_inv_abundance_results = simulate_attack(MHNO, 'inverse_abundance')
print("\tMUO")
muo_inv_abundance_results = simulate_attack(MUO, 'inverse_abundance')
print("\tMUNO")
muno_inv_abundance_results = simulate_attack(MUNO, 'inverse_abundance')


# Define attack results for Betweenness
betweenness_attack_results = {
    "MHO": mho_betweenness_results,
    "MHNO": mhno_betweenness_results,
    "MUO": muo_betweenness_results,
    "MUNO": muno_betweenness_results
}

# Define attack results for Degree
degree_attack_results = {
    "MHO": mho_degree_results,
    "MHNO": mhno_degree_results,
    "MUO": muo_degree_results,
    "MUNO": muno_degree_results
}

# Define attack results for Abundance
abundance_attack_results = {
    "MHO": mho_abundance_results,
    "MHNO": mhno_abundance_results,
    "MUO": muo_abundance_results,
    "MUNO": muno_abundance_results
}

# Define attack results for Inverse Abundance
inv_abundance_attack_results = {
    "MHO": mho_inv_abundance_results,
    "MHNO": mhno_inv_abundance_results,
    "MUO": muo_inv_abundance_results,
    "MUNO": muno_inv_abundance_results
}

# === Plot ===
colors = {
    "MHO": mho_col,
    "MHNO": mhno_col,
    "MUO": muo_col,
    "MUNO": muno_col
}
# Create a figure with subplots side by side for Betweenness and Degree attacks
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8, 3.25), dpi=300)
fig.subplots_adjust(wspace=0.2)

# --- Panel A: Targeted Attack (Degree Removal) ---


for network_name, results in degree_attack_results.items():
    x_vals_fraction = results['x_values']
    y_vals_normalized = results['y_values']
    fc_value = results['fc_value']
    i50_value = results['i50_value']
    color = colors[network_name]

    # Plot the decay curve
    ax1.plot(x_vals_fraction * 100, y_vals_normalized * 100,
             color=color, linewidth=2, label=network_name)

    # Add marker for critical threshold (f_c)
    y_fc = np.interp(fc_value * 100, x_vals_fraction * 100, y_vals_normalized * 100)
    ax1.plot(fc_value * 100, y_fc, marker=10, markersize=6,
        markerfacecolor='none',
        markeredgecolor=color,
        markeredgewidth=2,
        linestyle='None',
             #label=f"{network_name} (f_c) = {fc_value*100:.2f}%"
             )

    # Add marker for NR50 value
    y_i50 = np.interp(i50_value * 100, x_vals_fraction * 100, y_vals_normalized * 100)
    ax1.plot(i50_value * 100, y_i50, marker='o', markersize=6,
        markerfacecolor='none',
        markeredgecolor=color,
        markeredgewidth=2,
        linestyle='None',
             #label=f"{network_name} (NR50) = {i50_value*100:.2f}%"
             )

ax1.set_xlabel("% removed nodes", fontsize=12)
ax1.set_ylabel("% remaining nodes (LCC)", fontsize=12)
ax1.tick_params(axis='both', which='major', labelsize=10)
ax1.grid(True, linestyle='--', alpha=0.7)
legend1 = ax1.legend(loc="best")
legend1._legend_box.align = "left"
# ax1.legend().remove()


# --- Panel B: Targeted Attack (Betweenness Removal) ---


for network_name, results in betweenness_attack_results.items():
    x_vals_fraction = results['x_values']
    y_vals_normalized = results['y_values']
    fc_value = results['fc_value']
    i50_value = results['i50_value']
    color = colors[network_name]

    # Plot the decay curve
    ax2.plot(x_vals_fraction * 100, y_vals_normalized * 100,
             color=color, linewidth=2, label=network_name)

    # Add marker for critical threshold (f_c)
    y_fc = np.interp(fc_value * 100, x_vals_fraction * 100, y_vals_normalized * 100)
    ax2.plot(fc_value * 100, y_fc,
             marker=10, markersize=6,
        markerfacecolor='none',
        markeredgecolor=color,
        markeredgewidth=2,
        linestyle='None',
        # zorder=5
             #label=f"{network_name} (f_c) = {fc_value*100:.2f}%"
             )

    # Add marker for NR50 value
    y_i50 = np.interp(i50_value * 100, x_vals_fraction * 100, y_vals_normalized * 100)
    ax2.plot(i50_value * 100, y_i50,
             marker='o',
        markersize=6,
        markerfacecolor='none',
        markeredgecolor=color,
        markeredgewidth=2,
        linestyle='None',
        # zorder=5
            #  marker='o', markersize=8, color=color, fillstyle='none', linestyle='None',
             #label=f"{network_name} (NR50) = {i50_value*100:.2f}%"
             )


ax2.set_xlabel("% removed nodes", fontsize=12)
# ax2.set_ylabel("% remaining nodes (LCC)", fontsize=12)
ax2.tick_params(axis='both', which='major', labelsize=10)
ax2.grid(True, linestyle='--', alpha=0.7)
# ax2.legend(loc='upper center', bbox_to_anchor=(0.5, -0.2), ncol=4)
legend2 = ax2.legend(loc="best")
legend2._legend_box.align = "left"


plt.tight_layout()

plt.savefig("./betweenness_and_degree.png",
            format='png', dpi=300,
            bbox_inches='tight')



# Create a figure with subplots side by side for Abundance and Inverse Abundance attacks
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8, 3.25), dpi=300)
fig.subplots_adjust(wspace=0.2)

# --- Panel A: Targeted Attack (Abundance Removal) ---
for network_name, results in abundance_attack_results.items():
    x_vals_fraction = results['x_values']
    y_vals_normalized = results['y_values']
    fc_value = results['fc_value']
    i50_value = results['i50_value']
    color = colors[network_name]

    # Plot the decay curve
    ax1.plot(x_vals_fraction * 100, y_vals_normalized * 100,
             color=color, linewidth=2, label=network_name)

    # Add marker for critical threshold (f_c)
    y_fc = np.interp(fc_value * 100, x_vals_fraction * 100, y_vals_normalized * 100)
    ax1.plot(fc_value * 100, y_fc, marker=10, markersize=6,
        markerfacecolor='none',
        markeredgecolor=color,
        markeredgewidth=2,
        linestyle='None',
             zorder = 5
             #label=f"{network_name} (f_c) = {fc_value*100:.2f}%"
             )

    # Add marker for NR50 value
    y_i50 = np.interp(i50_value * 100, x_vals_fraction * 100, y_vals_normalized * 100)
    ax1.plot(i50_value * 100, y_i50, marker='o', markersize=6,
        markerfacecolor='none',
        markeredgecolor=color,
        markeredgewidth=2,
        linestyle='None',
             zorder = 5
             #label=f"{network_name} (NR50) = {i50_value*100:.2f}%"
             )


ax1.set_xlabel("% removed nodes", fontsize=12)
ax1.set_ylabel("% remaining nodes (LCC)", fontsize=12)
ax1.tick_params(axis='both', which='major', labelsize=10)
ax1.grid(True, linestyle='--', alpha=0.7)
# ax1.legend(loc='upper center', bbox_to_anchor=(0.5, -0.2), ncol=4)
legend1 = ax1.legend(loc="best")
legend1._legend_box.align = "left"

# --- Panel B: Targeted Attack (Inverse Abundance Removal) ---
for network_name, results in inv_abundance_attack_results.items():
    x_vals_fraction = results['x_values']
    y_vals_normalized = results['y_values']
    fc_value = results['fc_value']
    i50_value = results['i50_value']
    color = colors[network_name]

    # Plot the decay curve
    ax2.plot(x_vals_fraction * 100, y_vals_normalized * 100,
             color=color, linewidth=2, label=network_name)

    # Add marker for critical threshold (f_c)
    y_fc = np.interp(fc_value * 100, x_vals_fraction * 100, y_vals_normalized * 100)
    ax2.plot(fc_value * 100, y_fc,
             marker=10,
             markersize=6,
        markerfacecolor='none',
        markeredgecolor=color,
        markeredgewidth=2,
        linestyle='None',
             #label=f"{network_name} (f_c) = {fc_value*100:.2f}%"
             )

    # Add marker for NR50 value
    y_i50 = np.interp(i50_value * 100, x_vals_fraction * 100, y_vals_normalized * 100)
    ax2.plot(i50_value * 100, y_i50,
            marker = 'o',
             markersize=6,
        markerfacecolor='none',
        markeredgecolor=color,
        markeredgewidth=2,
        linestyle='None',
             #label=f"{network_name} (NR50) = {i50_value*100:.2f}%"
             )


ax2.set_xlabel("% removed nodes", fontsize=12)
# ax2.set_ylabel("% remaining nodes (LCC)", fontsize=12)
ax2.tick_params(axis='both', which='major', labelsize=10)
ax2.grid(True, linestyle='--', alpha=0.7)
legend2 = ax2.legend(loc="best")
legend2._legend_box.align = "left"
# ax2.legend().remove()

plt.tight_layout()

plt.savefig("./abundance.png",
            format='png', dpi=300,
            bbox_inches='tight')

plt.show()

# Create a dictionary to store the results in the desired format
results_data_reshaped = {}

attack_types = {
    "Betweenness Attack": betweenness_attack_results,
    "Degree Attack": degree_attack_results,
    "Abundance Attack": abundance_attack_results,
    "Inverse Abundance Attack": inv_abundance_attack_results,
}

for attack_type, results_dict in attack_types.items():
    for network_name, results in results_dict.items():
        row_name = f"{attack_type.split(' ')[0]} - {network_name}" # e.g., "Degree - MHNO"
        results_data_reshaped[row_name] = {
            "i50_value": results['i50_value'] * 100,  # Convert to percentage
            "fc_value": results['fc_value'] * 100  # Convert to percentage

        }

# Create the DataFrame from the reshaped data
results_df_reshaped = pd.DataFrame.from_dict(results_data_reshaped, orient='index')

results_df_reshaped.to_csv('/nr50_fc.csv', index=True)

# === RANDOM ATTACK FRAMEWORK ===
def remove_random_node(g, n, threshold = .5):
  lcc_sizes = []
  number_of_attacks = 0
  gOrder = g.order()
  lcc_size = g.order()
  for i in range(n):
    number_of_attacks += 1
    node = random.choice(list(g.nodes()))
    g.remove_node(node)
    try:
      lcc_size = len(max(nx.connected_components(g), key=len))
    except:
      lcc_size = 0
    lcc_sizes.append(lcc_size)

  return(lcc_sizes)

nreps = 1000
random.seed(505)

print("MHO")
MHO_lcc_sizes_random = []
while len(MHO_lcc_sizes_random) < nreps:
  MHO_kk = nx.Graph(MHO)
  MHO_lcc_sizes_random.append(remove_random_node(MHO_kk, MHO_kk.order()))

print("MHNO")
MHNO_lcc_sizes_random = []
while len(MHNO_lcc_sizes_random) < nreps:
  MHNO_kk = nx.Graph(MHNO)
  MHNO_lcc_sizes_random.append(remove_random_node(MHNO_kk, MHNO_kk.order()))

print("MUO")
MUO_lcc_sizes_random = []
while len(MUO_lcc_sizes_random) < nreps:
  MUO_kk = nx.Graph(MUO)
  MUO_lcc_sizes_random.append(remove_random_node(MUO_kk, MUO_kk.order()))

print("MUNO")
MUNO_lcc_sizes_random = []
while len(MUNO_lcc_sizes_random) < nreps:
  MUNO_kk = nx.Graph(MUNO)
  MUNO_lcc_sizes_random.append(remove_random_node(MUNO_kk, MUNO_kk.order()))

arr_MHO = np.array(MHO_lcc_sizes_random)
arr_MHNO = np.array(MHNO_lcc_sizes_random)
arr_MUO = np.array(MUO_lcc_sizes_random)
arr_MUNO = np.array(MUNO_lcc_sizes_random)

plot_MHO = np.mean(arr_MHO, axis = 0)
plot_MHNO = np.mean(arr_MHNO, axis = 0)
plot_MUO = np.mean(arr_MUO, axis = 0)
plot_MUNO = np.mean(arr_MUNO, axis = 0)

sd_MHO = np.std(arr_MHO, axis = 0)
sd_MHNO = np.std(arr_MHNO, axis = 0)
sd_MUO = np.std(arr_MUO, axis = 0)
sd_MUNO = np.std(arr_MUNO, axis = 0)

x_vals_MHO = [100 * i/MHO.order() for i in range(MHO.order())]
x_vals_MHNO = [100 * i/MHNO.order() for i in range(MHNO.order())]
x_vals_MUO = [100 * i/MUO.order() for i in range(MUO.order())]
x_vals_MUNO = [100 * i/MUNO.order() for i in range(MUNO.order())]

MHO_rand_to_plot = np.array([100 * i/MHO.order() for i in plot_MHO])
MHNO_rand_to_plot = np.array([100 * i/MHNO.order() for i in plot_MHNO])
MUO_rand_to_plot = np.array([100 * i/MUO.order() for i in plot_MUO])
MUNO_rand_to_plot = np.array([100 * i/MUNO.order() for i in plot_MUNO])

MHO_rand_to_plot_sd = np.array([100 * i/MHO.order() for i in sd_MHO])
MHNO_rand_to_plot_sd = np.array([100 * i/MHNO.order() for i in sd_MHNO])
MUO_rand_to_plot_sd =  np.array([100 * i/MUO.order() for i in sd_MUO])
MUNO_rand_to_plot_sd =  np.array([100 * i/MUNO.order() for i in sd_MUNO])

def find_i50(x_vals_fraction, lcc_sizes_normalized_fraction):
  i50 = []
  for i,j in zip(x_vals_fraction, lcc_sizes_normalized_fraction):
    if j>=50:
      i50.append(i)
  return(i50[-1] if i50 else 0)

def calculate_critical_threshold(x_vals_fraction, lcc_sizes_normalized_fraction):
    if len(lcc_sizes_normalized_fraction) <= 1:
        return 0.0 
    
    slopes = np.diff(lcc_sizes_normalized_fraction)

    if len(slopes) == 0:
        return 0.0

    critical_index = np.argmin(slopes)

    f_c = x_vals_fraction[critical_index + 1]

    return f_c


results = {
    'MHO': {'fc': [], 'nr50': []},
    'MHNO': {'fc': [], 'nr50': []},
    'MUO': {'fc': [], 'nr50': []},
    'MUNO': {'fc': [], 'nr50': []},
}

networks_data = {
    'MHO': {'lcc_sizes': MHO_lcc_sizes_random, 'order': MHO.order()},
    'MHNO': {'lcc_sizes': MHNO_lcc_sizes_random, 'order': MHNO.order()},
    'MUO': {'lcc_sizes': MUO_lcc_sizes_random, 'order': MUO.order()},
    'MUNO': {'lcc_sizes': MUNO_lcc_sizes_random, 'order': MUNO.order()},
}

for net_name, data in networks_data.items():
    print(f"Calculating percolation threshold & nr50 for: {net_name}")
    for lcc_sizes in data['lcc_sizes']:

        x_vals_fraction = np.arange(1, len(lcc_sizes) + 1) / data['order']
        lcc_sizes_normalized = np.array(lcc_sizes) / data['order']
        lcc_sizes_normalized_fraction = lcc_sizes_normalized * 100

        fc = calculate_critical_threshold(x_vals_fraction, lcc_sizes_normalized)
        results[net_name]['fc'].append(fc)

        nr50 = find_i50(x_vals_fraction, lcc_sizes_normalized_fraction)
        results[net_name]['nr50'].append(nr50)


all_fc = {
    'MHNO': results['MHNO']['fc'],
    'MHO': results['MHO']['fc'],
    'MUNO': results['MUNO']['fc'],
    'MUO': results['MUO']['fc'],
}
df_fc = pd.DataFrame(all_fc)
df_fc.to_csv('./random_fc.csv', index=False)
df_fc.head()

all_nr50 = {
    'MHNO': results['MHNO']['nr50'],
    'MHO': results['MHO']['nr50'],
    'MUNO': results['MUNO']['nr50'],
    'MUO': results['MUO']['nr50'],
}
df_nr50 = pd.DataFrame(all_nr50)
df_nr50.to_csv('./random_nr50.csv', index=False)
df_nr50.head()

import seaborn as sns
sns.set(style="whitegrid")
my_pal = {"MHO" : mho_col,
          "MHNO" : mhno_col,
          "MUO" : muo_col,
          "MUNO" : muno_col}

plt.style.use("default")


fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8, 3.25), dpi=300) 
fig.subplots_adjust(wspace=0.2)

ax1.plot(x_vals_MHNO, MHNO_rand_to_plot, linewidth=2, label="MHNO", color=mhno_col)
ax1.fill_between(x_vals_MHNO,
                 MHNO_rand_to_plot - MHNO_rand_to_plot_sd,
                 MHNO_rand_to_plot + MHNO_rand_to_plot_sd,
                 edgecolor="none", color=mhno_col, alpha=0.2)

ax1.plot(x_vals_MHO, MHO_rand_to_plot, linewidth=2, label="MHO", color=mho_col)
ax1.fill_between(x_vals_MHO,
                 MHO_rand_to_plot - MHO_rand_to_plot_sd,
                 MHO_rand_to_plot + MHO_rand_to_plot_sd,
                 edgecolor="none", color=mho_col, alpha=0.2)

ax1.plot(x_vals_MUNO, MUNO_rand_to_plot, linewidth=2, label="MUNO", color=muno_col)
ax1.fill_between(x_vals_MUNO,
                 MUNO_rand_to_plot - MUNO_rand_to_plot_sd,
                 MUNO_rand_to_plot + MUNO_rand_to_plot_sd,
                 edgecolor="none", color=muno_col, alpha=0.2)

ax1.plot(x_vals_MUO, MUO_rand_to_plot, linewidth=2, label="MUO", color=muo_col)
ax1.fill_between(x_vals_MUO,
                 MUO_rand_to_plot - MUO_rand_to_plot_sd,
                 MUO_rand_to_plot + MUO_rand_to_plot_sd,
                 edgecolor="none", color=muo_col, alpha=0.2)


ax1.set_xlabel("% removed nodes", fontsize=12)
ax1.set_ylabel("% remaining nodes (LCC)", fontsize=12)
ax1.tick_params(axis='both', which='major', labelsize=10)
ax1.legend(loc="best")
ax1.grid(True, linestyle='--', alpha=0.7)



df_nr50_percent = df_nr50 * 100

sns.set_style("whitegrid")
sns.stripplot(
    data=df_nr50_percent,
    palette=my_pal,
    alpha=0.5,
    size=3,
    jitter=0.1,
    ax=ax2  
)
sns.boxplot(
    data=df_nr50_percent,
    palette=my_pal,
    showfliers=False,
    width=0.5,
    boxprops={'alpha': 0.5},
    # whis=[0, 100],
    ax=ax2 
)

ax2.grid(True, linestyle='--', alpha=0.7)


ax2.set_xlabel("Phenotype", fontsize=12)
ax2.set_ylabel("$NR_{50}$ values", fontsize=12)

ax2.tick_params(axis='both', which='major', labelsize=10)
y_ticks = [30, 35, 40, 45, 50]
ax2.set_yticks(y_ticks)

ax2.set_ylim(27, 52)


plt.tight_layout()


plt.savefig("/random.png",
            format='png', dpi=1200,
            bbox_inches='tight')  
plt.show();
