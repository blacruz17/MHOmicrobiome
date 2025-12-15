# -*- coding: utf-8 -*-
import os, glob, math
import pandas as pd
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import kruskal


base_path = '.'
output_path = os.path.join(base_path, "results")
os.makedirs(output_path, exist_ok=True)

# === GROUPS AND ATTACKS ===
groups = ["MHNO", "MUO", "MUNO"]
attack_types = ["degree", "betweenness", "abundance", "inverse_abundance", "random"]

# === COLORS ===
colors = {
    "MHNO": "#264653",
    "MUO": "#edafb8",
    "MUNO": "#703d67"
}

# === HELPER FUNCTIONS ===

def get_connected_graph(g):
    if nx.is_connected(g):
        return g
    else:
        largest_cc = max(nx.connected_components(g), key=len)
        return g.subgraph(largest_cc).copy()

def add_abundance(net, metadata):
          abundance_dict = pd.Series(metadata.Abundance.values, index=metadata.Label).to_dict()
          nx.set_node_attributes(net, abundance_dict, 'Abundance')

def add_weights(net, weights):
  for _, row in weights.iterrows():
      v1, v2, w = row["OTU_1"], row["OTU_2"], row["weight"]
      if net.has_edge(v1, v2):
          net[v1][v2]["weight"] = w
          net[v1][v2]["invWeight"] = 1- w
      else:
          print(f"edge {v1},{v2} not found")


def filter_positive_edges(net):
    G_filtered = net.__class__()

    for n, attrs in net.nodes(data=True):
        G_filtered.add_node(n, **attrs)

    for u, v, data in net.edges(data=True):
        w = data.get("weight", None)
        if w is not None and w > 0:
            G_filtered.add_edge(u, v, **data)

    return G_filtered

def get_pos_edges(net):
  pos_edges_net, neg_edges_net = [], []
  for edge in net.edges(data = True):
    if edge[-1]['weight'] > 0:
      pos_edges_net.append(edge[-1]['weight'])

  return(100 * len(pos_edges_net)/net.size())

def calculate_critical_threshold(x_vals_fraction, lcc_sizes_normalized_fraction):

    if len(lcc_sizes_normalized_fraction) <= 1:
        return 0.0 

    slopes = np.diff(lcc_sizes_normalized_fraction)

    if len(slopes) == 0:
        return 0.0

    critical_index = np.argmin(slopes)
    
    f_c = x_vals_fraction[critical_index + 1]

    return f_c


def find_i50(x_vals_fraction, lcc_sizes_normalized_fraction):
  i50 = []
  for i,j in zip(x_vals_fraction, lcc_sizes_normalized_fraction):
    if j>=0.50:
      i50.append(i)
  return(i50[-1])

def attack_degree(g_original, n):

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

    try:
        node = max(g.degree, key=lambda x: x[1])[0]
    except ValueError: 
        lcc_sizes.append(0)
        g.clear() 
        continue

    g.remove_node(node)

    try:
      if g.nodes(): 
        lcc_size = len(max(nx.connected_components(g), key=len))
      else:
        lcc_size = 0 
    except ValueError: 
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

    betweenness_dict = nx.betweenness_centrality(g, weight="invWeight")

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

  lcc_sizes.append(initial_lcc_size) # Añadir estado inicial

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

def simulate_random_attack(g, n_reps=10):
    n_nodes = len(g)
    y_values_all = []
    i50_list = []
    fc_list = []

    for _ in range(n_reps):
        g_temp = g.copy()
        nodes = list(g_temp.nodes())
        np.random.shuffle(nodes)

        y_values = []
        for node in nodes:
            g_temp.remove_node(node)
            if len(g_temp) == 0:
                y_values.append(0)
            else:
                largest_cc = len(max(nx.connected_components(g_temp), key=len))
                y_values.append(largest_cc / n_nodes)

        y_values = np.array(y_values)
        y_values_all.append(y_values)

        x_vals_fraction = np.arange(1, n_nodes + 1) / n_nodes

        try:
            i50 = find_i50(x_vals_fraction, y_values)
        except IndexError:
            i50 = 1.0  

        fc = calculate_critical_threshold(x_vals_fraction, y_values)

        i50_list.append(i50 * 100)
        fc_list.append(fc * 100)

    return {
        "y_values_mean": np.mean(y_values_all, axis=0),
        "i50_mean": np.mean(i50_list),
        "fc_mean": np.mean(fc_list)
    }

def simulate_attack(g_original, attack_type):

  g = g_original.copy() # Work on a copy to preserve the original graph
  n = g.order() # Number of nodes to remove

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

  # Calculate initial LCC size from the first element of lcc_sizes
  initial_lcc_size = lcc_sizes[0] if lcc_sizes[0] > 0 else 1

  # Normalize LCC sizes by the initial size
  lcc_sizes_normalized_fraction = lcc_sizes / initial_lcc_size

  # Calculate x-values (fraction of nodes removed)
  x_vals_fraction = np.array([i / n for i in range(n + 1)])

  # Calculate critical threshold (f_c) and i50
  fc_value = calculate_critical_threshold(x_vals_fraction, lcc_sizes_normalized_fraction)
  i50_value = find_i50(x_vals_fraction, lcc_sizes_normalized_fraction)


  return {
      'x_values': x_vals_fraction,
      'y_values': lcc_sizes_normalized_fraction,
      'fc_value': fc_value,
      'i50_value': i50_value
  }


# === PART 1: LOAD NETWORKS ===

networks_ready = []  


all_files = sorted(glob.glob(os.path.join(base_path, "*.graphml")))

for path in all_files:
    try:
        file_name = os.path.basename(path).replace(".graphml", "")

        group = next((g for g in groups if g in file_name.upper()), None)
        if group is None:
            print(f"group for {file_name} not found, continue.")
            continue

        g = nx.read_graphml(path)
        g = get_connected_graph(g)
        mapping = {node: data['name'] for node, data in g.nodes(data=True)}

        nx.relabel_nodes(g, mapping, copy=False)

        metadata_path = path.replace("_ig_mb.graphml", "_metadata.csv")
        if os.path.exists(metadata_path):
            metadata = pd.read_csv(metadata_path)
            add_abundance(g, metadata)
        else:
            print(f"metadata not found {file_name}: {metadata_path}")


        weights_path = path.replace("_ig_mb.graphml", "_edge_weights.csv")
        if os.path.exists(weights_path):
            weights = pd.read_csv(weights_path)
            add_weights(g, weights)
            g = filter_positive_edges(g)
        else:
            print(f"weights not found {file_name}: {weights_path}")

        networks_ready.append({
            "file": file_name,
            "group": group,
            "graph": g
        })

    except Exception as e:
        print(f"error loading {path}: {e}")

print(f"\n{len(networks_ready)} networks ready for attack")

for item in networks_ready[:2]:
    print(item["file"], list(item["graph"].nodes(data=True))[:3])

# === PART 2: ATTACKS SIMULATIONS ===

results = {g: {a: {"y_arrays": [], "nr50": [], "fc": []} for a in attack_types} for g in groups}
metric_rows = []

for item in networks_ready:
    file_name = item["file"]
    group = item["group"]
    g = item["graph"]

    print(f"\n########################\nperforming attacks: {file_name} ({group})")

    try:
        # --- Targeted attacks ---
        for ataque in ["degree", "betweenness", "abundance", "inverse_abundance"]:
            print(f" → Ataque: {ataque}")
            res = simulate_attack(g, ataque)
            results[group][ataque]["y_arrays"].append(res["y_values"])
            results[group][ataque]["nr50"].append(res["i50_value"] * 100)
            results[group][ataque]["fc"].append(res["fc_value"] * 100)
            metric_rows.append({
                "file": file_name,
                "MetObesity": group,
                "attack_type": ataque,
                "NR50": res["i50_value"] * 100,
                "f_c": res["fc_value"] * 100
            })

        # --- Random attacks (10 reps) ---
        res_rand = simulate_random_attack(g, n_reps=10)
        results[group]["random"]["y_arrays"].append(res_rand["y_values_mean"])
        results[group]["random"]["nr50"].append(res_rand["i50_mean"])
        results[group]["random"]["fc"].append(res_rand["fc_mean"])
        metric_rows.append({
            "file": file_name,
            "MetObesity": group,
            "attack_type": "random",
            "NR50": res_rand["i50_mean"],
            "f_c": res_rand["fc_mean"]
        })

    except Exception as e:
        print(f"error {file_name}: {e}")


# === EXPORT results ===
results_df = pd.DataFrame(metric_rows)
csv_path = os.path.join(output_path, "attack_metrics_subsampledNets.csv")
results_df.to_csv(csv_path, index=False)


# == Plots ===
### Calculate mean $NR_{50}$ and $f_c$ values for the MHNO group to add them to the plot:

mhno_df = results_df[results_df["MetObesity"] == "MHNO"]

mean_lines = (
    mhno_df.groupby("attack_type")[["NR50", "f_c"]]
           .mean()
           .reset_index()
)
mean_dict = mean_lines.set_index("attack_type").to_dict(orient="index")

mean_lines

## Plot all curves in a grid:
n_rows = len(groups)
n_cols = len(attack_types)

fig_width = 3.8 * n_cols
fig_height = 3.0 * n_rows
fig = plt.figure(figsize=(fig_width, fig_height), dpi=300)

axes = fig.subplots(
    n_rows, n_cols,
    sharex=True,
    sharey=True
)

axes = np.atleast_2d(axes)

for i, group in enumerate(["MHNO", "MUNO", "MUO"]):
    for j, ataque in enumerate(attack_types):
        ax = axes[i, j]

        y_arrays = results[group][ataque]["y_arrays"]
        nr50_vals = results[group][ataque]["nr50"]
        fc_vals = results[group][ataque]["fc"]

        if len(y_arrays) == 0:
            ax.set_visible(False)
            continue

        min_len = min(len(y) for y in y_arrays)
        y_trim = [y[:min_len] for y in y_arrays]
        x_vals = np.linspace(0, 100, min_len)

        for idx, y in enumerate(y_trim):
            ax.plot(
                x_vals, y * 100,
                color=colors[group],
                linewidth=1.0,
                alpha=0.65
            )

            if idx < len(nr50_vals) and idx < len(fc_vals):
                nr50 = nr50_vals[idx]
                fc = fc_vals[idx]

                ax.scatter(nr50, np.interp(nr50, x_vals, y * 100),
                           color=colors[group], s=14, marker='o', zorder=3,
                           edgecolor='black', linewidth=0.4)

                ax.scatter(fc, np.interp(fc, x_vals, y * 100),
                           color=colors[group], s=14, marker='x', zorder=3)

        # add vertical lines:
        if ataque in mean_dict:
          nr50_ref = mean_dict[ataque]["NR50"]
          fc_ref   = mean_dict[ataque]["f_c"]

          ax.axvline(nr50_ref, color="black", linestyle="--", linewidth=1.2)
          ax.axvline(fc_ref, color="black", linestyle=":", linewidth=1.2)

        if i == n_rows - 1:
            ax.set_xlabel("% removed nodes", fontsize=14)

        if j == 0:
            ax.set_ylabel("% remaining nodes (LCC)", fontsize=14)

        ax.set_xlim(0, 100)
        ax.set_ylim(0, 100)
        ax.grid(True, linestyle="--", alpha=0.4)

        ax.tick_params(axis='both', labelsize=12)
        if i == 0:
            titulo = ataque.capitalize().replace("_abundance", " abundance")
            ax.set_title(titulo, fontsize=16, pad=10)

plt.tight_layout(h_pad=2.0, w_pad=1.0)
save_path = os.path.join(output_path, "subsampledNetworks_attacks.png")
plt.savefig(save_path, dpi=300, bbox_inches='tight')
plt.show()

palette = {
    "MHO": "#2a9d8f",
    "MHNO": "#264653",
    "MUO": "#edafb8",
    "MUNO": "#703d67"
}
