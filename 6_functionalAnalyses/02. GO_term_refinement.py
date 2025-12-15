#!/usr/bin/env python3
"""
GO term post-processing for HUMAnN outputs.

"""
# import standard libraries
import json
from collections import defaultdict, deque

import pandas as pd



# 1) Load HUMAnN GO-mapped gene families table and exclude noise rows

df = pd.read_csv("Data/genefamilies_merged_unstrat_to_go.tsv", sep="\t")
col_gf = "# Gene Family"

df_go = df[~df[col_gf].str.upper().isin({"UNMAPPED", "UNGROUPED"})]
go_ids = df_go[col_gf].tolist() # List of GO IDs in the data


# 2) Load GO graph (go-basic.json) and define helpers for node metadata

with open("Data/go-basic.json", "r", encoding="utf-8") as f:
    gob = json.load(f)

graph = gob["graphs"][0]  # The file contains a single graph with all GO nodes.

IAO_REPLACED_BY = "http://purl.obolibrary.org/obo/IAO_0100001"
OBO_CONSIDER = "http://www.geneontology.org/formats/oboInOwl#consider"


def get_aspect_from_node(node: dict):
    """Return the GO aspect (namespace) for the given node, if available."""
    meta = node.get("meta", {})
    for bpv in meta.get("basicPropertyValues", []):
        pred = bpv.get("pred") or bpv.get("property")
        if pred and "hasOBONamespace" in pred:
            val = bpv.get("val") or bpv.get("value")
            if val in ("biological_process", "molecular_function", "cellular_component"):
                return val
    return None


def is_deprecated(node: dict) -> bool:
    """Whether the GO node is marked as deprecated in go-basic.json."""
    meta = node.get("meta", {})
    return bool(meta.get("deprecated", False))


def iri_to_go(iri):
    """Convert GO IRIs to GO: identifiers; preserve GO: identifiers; otherwise return None."""
    if isinstance(iri, str) and iri.startswith("http://purl.obolibrary.org/obo/GO_"):
        return iri.replace("http://purl.obolibrary.org/obo/GO_", "GO:")
    return iri if isinstance(iri, str) and iri.startswith("GO:") else None


def basic_prop_values(node: dict):
    """Return the list of basicPropertyValues for a GO node."""
    return node.get("meta", {}).get("basicPropertyValues", []) or []


def get_first_basic(node: dict, pred_iri: str):
    """Return the first basicPropertyValues 'val' matching a predicate IRI."""
    for bp in basic_prop_values(node):
        if bp.get("pred") == pred_iri:
            return bp.get("val")
    return None


def get_all_basic(node: dict, pred_iri: str):
    """Return all basicPropertyValues 'val' entries matching a predicate IRI."""
    return [bp.get("val") for bp in basic_prop_values(node) if bp.get("pred") == pred_iri]


# ------------------------------------------------------------
# 3) Build a GO node index and a resolver for obsolete GO IDs
# ------------------------------------------------------------
id_index = {}
for n in graph.get("nodes", []):
    nid = n.get("id", "")
    if isinstance(nid, str) and nid.startswith("http://purl.obolibrary.org/obo/GO_"):
        gid = iri_to_go(nid)
        if gid:
            id_index[gid] = n


def resolve_go(go_id: str):
    """
    Resolve deprecated GO identifiers.
    - Prefer 'replaced_by' (exact replacement).
    - Otherwise try 'consider' and choose the first non-deprecated candidate.
    - If no valid replacement exists, return None so the term can be omitted.
    """
    cur = go_id
    seen = set()

    while True:
        node = id_index.get(cur)
        if not node:
            return cur  # Not found in the graph: keep as-is.

        if not is_deprecated(node):
            return cur  # Current term: keep.

        # 1) replaced_by (exact replacement)
        rep = get_first_basic(node, IAO_REPLACED_BY)
        if rep:
            nxt = iri_to_go(rep) if isinstance(rep, str) else None
            if nxt and nxt not in seen:
                seen.add(nxt)
                cur = nxt
                continue

        # 2) consider (suggestions; take the first current one)
        considers = get_all_basic(node, OBO_CONSIDER)
        for c in considers:
            cand = iri_to_go(c) if isinstance(c, str) else None
            if not cand:
                continue
            cand_node = id_index.get(cand)
            if cand_node and not is_deprecated(cand_node):
                return cand

        # 3) Deprecated with no valid replacement
        return None


# 4) Extract GO node metadata (mapped ID, label, aspect) for annotation

nodes = {}
for n in graph.get("nodes", []):
    nid = n.get("id", "")
    if isinstance(nid, str) and nid.startswith("http://purl.obolibrary.org/obo/GO_"):
        go_id = nid.replace("http://purl.obolibrary.org/obo/GO_", "GO:")

        mapped_id = resolve_go(go_id)
        if mapped_id is None:
            continue  # Omit deprecated terms with no replacement.

        mapped_node = id_index.get(mapped_id, n)
        nodes[go_id] = {
            "GO_ID_map": mapped_id,
            "name": mapped_node.get("lbl"),
            "aspect": get_aspect_from_node(mapped_node),
        }

meta = pd.DataFrame.from_dict(nodes, orient="index").reset_index().rename(
    columns={"index": "GO_ID", "name": "GO_name"}
)

meta_idx = meta.set_index("GO_ID")
df_go = df_go.join(meta_idx, on=col_gf)

df_go = df_go[[col_gf, "GO_ID_map", "GO_name", "aspect"]]
df_go = df_go.dropna(subset=["GO_ID_map", "GO_name", "aspect"])


# 5) Build the within-aspect 'is_a' DAG and retain only GO terms present in the data
IS_A = {
    "is_a",
    "rdfs:subClassOf",
    "http://www.w3.org/2000/01/rdf-schema#subClassOf",
}

parents = defaultdict(set)   # child -> {parent}
children = defaultdict(set)  # parent -> {child}

for e in graph.get("edges", []):
    child = e.get("sub")
    parent = e.get("obj")

    child = (
        child.replace("http://purl.obolibrary.org/obo/GO_", "GO:")
        if isinstance(child, str) and child.startswith("http://purl.obolibrary.org/obo/GO_")
        else None
    )
    parent = (
        parent.replace("http://purl.obolibrary.org/obo/GO_", "GO:")
        if isinstance(parent, str) and parent.startswith("http://purl.obolibrary.org/obo/GO_")
        else None
    )

    if not child or not parent:
        continue
    if child not in nodes or parent not in nodes:
        continue
    if nodes[child]["aspect"] != nodes[parent]["aspect"]:
        continue

    pred = e.get("pred") or ""
    if pred in IS_A or (isinstance(pred, str) and pred.endswith("is_a")):
        parents[child].add(parent)
        children[parent].add(child)


def ancestors(go_id):
    """Return the ancestor set of a GO term within the constructed DAG."""
    seen = set()
    q = deque([go_id])
    while q:
        cur = q.popleft()
        for p in parents.get(cur, []):
            if p not in seen:
                seen.add(p)
                q.append(p)
    return seen


present_go = {g for g in df_go[col_gf].unique() if g in nodes}

anc2children = defaultdict(set)
for g in present_go:
    for anc in ancestors(g):
        anc2children[anc].add(g)

children_count = {go: len(children.get(go, set())) for go in nodes}
# Define the maximum number of children allowed for a term to be treated as a “leaf” candidate.
MAX_CHILDREN = 0
candidate_go = [g for g in present_go if children_count.get(g, 0) <= MAX_CHILDREN]


def anc_popularity(anc_id):
    """How many observed GO terms are subsumed by the ancestor."""
    return len(anc2children.get(anc_id, set()))


def parents_same_aspect(go_id):
    asp = nodes[go_id]["aspect"]
    return [p for p in parents.get(go_id, []) if nodes.get(p, {}).get("aspect") == asp]


def ancestors_same_aspect(go_id):
    asp = nodes[go_id]["aspect"]
    return [a for a in ancestors(go_id) if nodes.get(a, {}).get("aspect") == asp]


def pick_immediate_ancestor(go_id):
    """Choose an immediate parent, prioritising broadly-used ancestors."""
    ps = parents_same_aspect(go_id)
    if not ps:
        return None
    return sorted(ps, key=lambda a: (-anc_popularity(a), nodes[a]["name"] or ""))[0]


def compute_depths(children_map, roots):
    """Compute depth-from-root (per aspect) using a breadth-first traversal."""
    depth_map = {}
    q = deque([(r, 0) for r in roots])
    seen = set(roots)
    while q:
        cur, d = q.popleft()
        depth_map[cur] = d
        for ch in children_map.get(cur, []):
            if ch not in seen:
                seen.add(ch)
                q.append((ch, d + 1))
    return depth_map


roots = ["GO:0008150", "GO:0003674", "GO:0005575"]  # BP, MF, CC
depth = compute_depths(children, roots)


def pick_most_global_ancestor(go_id, allow_root=True):
    """
    Choose a more global ancestor within the same aspect.
    Lower depth (closer to the root) is preferred; ties are broken by popularity then name.
    """
    ancs = ancestors_same_aspect(go_id)
    if not ancs:
        return None
    if not allow_root:
        ancs = [a for a in ancs if depth.get(a, 1) > 0]
        if not ancs:
            return None
    return sorted(
        ancs,
        key=lambda a: (depth.get(a, float("inf")), -anc_popularity(a), nodes[a]["name"] or ""),
    )[0]



# 6) Assemble the output table and export to CSV and Excel

rows = []
for g in sorted(candidate_go):
    asp = nodes[g]["aspect"]
    imm = pick_immediate_ancestor(g)
    glob = pick_most_global_ancestor(g, allow_root=True)

    rows.append(
        {
            "GO_ID": g,
            "GO_name": nodes[g]["name"],
            "aspect": asp,
            "children_count": children_count.get(g, 0),
            "depth": depth.get(g, 0),
            "Immediate_parent_GO": imm,
            "Immediate_parent_name": nodes[imm]["name"] if imm else None,
            "Global_Ancestor_GO": glob,
            "Global_Ancestor_name": nodes[glob]["name"] if glob else None,
        }
    )

df_go_leaves_with_anc = pd.DataFrame(rows)

out_csv = "Data/my_GO_leaves_with_immediate_and_global_ancestor.csv"
df_go_leaves_with_anc.to_csv(out_csv, index=False, encoding="utf-8-sig")

out_xlsx = "Data/my_GO_leaves_with_immediate_and_global_ancestor.xlsx"
with pd.ExcelWriter(out_xlsx, engine="xlsxwriter") as writer:
    df_go_leaves_with_anc.to_excel(writer, index=False, sheet_name="GO_leaves")
