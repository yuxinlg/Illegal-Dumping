"""
02b_network_connectivity.py
===========================
Download and preprocess Philadelphia street-network connectivity metrics.

This script is COMPUTATIONALLY EXPENSIVE (expect 2–6 hours on a laptop due to
the 1 000 m ego-graph step on ~121 000 nodes).  Run it ONCE; 02_covariates.py
loads the resulting file directly.

Data source
-----------
OpenStreetMap walking network, accessed via OSMnx (October 2025 snapshot used
in the paper).  Re-running will pull the current OSM state and may produce
slightly different results.

Algorithm
---------
  1. Download Philadelphia County walking network (OSMnx).
  2. Convert to igraph for fast GLOBAL metrics:
       - Degree centrality
       - PageRank   (damping = 0.85, length-weighted; Page & Brin 1998)
       - Closeness centrality  (normalized, length-weighted)
  3. Sampled betweenness centrality via NetworkX
       (k = 12 000 nodes ≈ 10 % sample; seed = 42 for reproducibility).
  4. LOCAL (1 000 m ego-graph) closeness and betweenness for each node,
       computed in parallel batches with joblib.
  5. Attach all node attributes to the graph.
  6. Convert graph to edge GeoDataFrame; each edge metric = mean of its
       two incident node values.
  7. Export to output/  (GeoJSON, CSV, Shapefile).

Output files
------------
  output/philadelphia_walk_network_connectivity.geojson  ← loaded by 02_covariates.py
  output/philadelphia_walk_network_connectivity.csv
  output/philadelphia_walk_network_connectivity.shp / .dbf / .shx / .prj

Dependencies
------------
  osmnx, networkx, igraph, geopandas, joblib, tqdm, numpy
"""

import os
import time
import pickle
import warnings

import numpy as np
import networkx as nx
import igraph as ig
import osmnx as ox
import geopandas as gpd
from tqdm import tqdm
from joblib import Parallel, delayed

warnings.filterwarnings("ignore")

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
REPL_DIR = os.path.dirname(os.path.abspath(__file__))   # replication/
OUT      = os.path.join(REPL_DIR, "output")             # replication/output/
os.makedirs(OUT, exist_ok=True)

GRAPH_PKL   = os.path.join(OUT, "_graph_walk_step2.pkl")   # temporary; deleted after use
OUT_GEOJSON = os.path.join(OUT, "philadelphia_walk_network_connectivity.geojson")
OUT_SHP     = os.path.join(OUT, "philadelphia_walk_network_connectivity.shp")
OUT_CSV     = os.path.join(OUT, "philadelphia_walk_network_connectivity.csv")


# ---------------------------------------------------------------------------
# Helper: NetworkX → igraph conversion
# ---------------------------------------------------------------------------
def networkx_to_igraph(G_nx):
    """Convert an undirected NetworkX graph to an igraph Graph.

    Returns
    -------
    G_ig       : igraph.Graph  (undirected, edge weight = 'length')
    node_to_idx: dict  {osmid: igraph_index}
    idx_to_node: dict  {igraph_index: osmid}
    """
    nodes = list(G_nx.nodes())
    node_to_idx = {node: idx for idx, node in enumerate(nodes)}
    idx_to_node = {idx: node for node, idx in node_to_idx.items()}

    edges   = [(node_to_idx[u], node_to_idx[v]) for u, v in G_nx.edges()]
    weights = [data.get("length", 1) for _, _, data in G_nx.edges(data=True)]

    G_ig = ig.Graph(n=len(nodes), edges=edges, directed=False)
    G_ig.es["weight"] = weights
    return G_ig, node_to_idx, idx_to_node


# ---------------------------------------------------------------------------
# Helper: ego-graph metrics for a batch of nodes  (used by joblib)
# ---------------------------------------------------------------------------
def _ego_batch(batch_data):
    """Compute 1 000 m closeness and betweenness for a list of nodes.

    For large local graphs (> 500 nodes) betweenness is approximated with
    k = 100 pivots; otherwise exact computation is used.
    """
    nodes_batch, G_nx = batch_data
    results = []
    for node in nodes_batch:
        try:
            ego = nx.ego_graph(G_nx, node, radius=1000, distance="length")
            if len(ego) > 1:
                cl = nx.closeness_centrality(ego, u=node, distance="length")
                if len(ego) > 500:
                    bt_dict = nx.betweenness_centrality(
                        ego, k=100, weight="length", normalized=True, seed=42
                    )
                else:
                    bt_dict = nx.betweenness_centrality(
                        ego, weight="length", normalized=True
                    )
                bt = bt_dict.get(node, 0)
            else:
                cl = bt = 0
        except Exception:
            cl = bt = 0
        results.append((node, cl, bt))
    return results


# ---------------------------------------------------------------------------
# Main pipeline
# ---------------------------------------------------------------------------
def main():
    total_start = time.time()
    print("=" * 60)
    print("Network Connectivity Preprocessing  (walking network)")
    print("=" * 60)

    # ------------------------------------------------------------------
    # STEP 1  Download Philadelphia walking network
    # ------------------------------------------------------------------
    print("\n[1/8] Downloading Philadelphia County walking network (OSMnx) ...")
    t0 = time.time()
    G = ox.graph_from_place(
        "Philadelphia County, Pennsylvania, USA", network_type="walk"
    )
    print(f"  {len(G.nodes):,} nodes, {len(G.edges):,} edges  "
          f"({time.time() - t0:.0f} s)")

    # ------------------------------------------------------------------
    # STEP 2  Global centrality metrics via igraph
    # ------------------------------------------------------------------
    print("\n[2/8] Global centrality (igraph) ...")
    t0 = time.time()
    G_ig, _, idx_to_node = networkx_to_igraph(G)

    print("  degree ...", end=" ", flush=True)
    degree = {idx_to_node[i]: v for i, v in enumerate(G_ig.degree())}
    print("done")

    print("  PageRank (damping=0.85, length-weighted) ...", end=" ", flush=True)
    pagerank = {
        idx_to_node[i]: v
        for i, v in enumerate(
            G_ig.pagerank(weights="weight", damping=0.85, implementation="prpack")
        )
    }
    print("done")

    print("  closeness (normalized, length-weighted) ...", end=" ", flush=True)
    closeness = {
        idx_to_node[i]: v
        for i, v in enumerate(G_ig.closeness(weights="weight", normalized=True))
    }
    print("done")

    # ------------------------------------------------------------------
    # STEP 3  Global betweenness via sampled NetworkX
    # ------------------------------------------------------------------
    print("\n[3/8] Global betweenness (NetworkX, k=12 000 sample, seed=42) ...")
    t0 = time.time()
    betweenness = nx.betweenness_centrality(
        G, k=12_000, weight="length", normalized=True, seed=42
    )
    print(f"  done ({time.time() - t0:.0f} s)")

    # ------------------------------------------------------------------
    # STEP 4  Save graph for ego-graph step
    # ------------------------------------------------------------------
    print(f"\n[4/8] Pickling graph for ego-graph step → {GRAPH_PKL}")
    with open(GRAPH_PKL, "wb") as f:
        pickle.dump(G, f)

    # ------------------------------------------------------------------
    # STEP 5  1 000 m ego-graph metrics (parallel, the slow step)
    # ------------------------------------------------------------------
    print("\n[5/8] 1 000 m ego-graph metrics (parallel) — this is the slow step ...")
    print("  Estimated runtime: 2–6 h depending on hardware.\n")
    t0 = time.time()

    nodes_list = list(G.nodes())
    n_cores    = os.cpu_count()
    batch_size = max(100, len(nodes_list) // (n_cores * 20))
    batches    = [
        (nodes_list[i : i + batch_size], G)
        for i in range(0, len(nodes_list), batch_size)
    ]

    raw_results = Parallel(n_jobs=n_cores)(
        delayed(_ego_batch)(b) for b in tqdm(batches, desc="  processing batches")
    )

    closeness_1000 = {}
    betweenness_1000 = {}
    for batch_res in raw_results:
        for node, cl, bt in batch_res:
            closeness_1000[node]  = cl
            betweenness_1000[node] = bt

    # fill any nodes missed by the parallel run
    for node in nodes_list:
        closeness_1000.setdefault(node, 0)
        betweenness_1000.setdefault(node, 0)

    print(f"  1 000 m metrics done  ({(time.time() - t0) / 60:.1f} min)")

    # ------------------------------------------------------------------
    # STEP 6  Attach all attributes to graph nodes
    # ------------------------------------------------------------------
    print("\n[6/8] Attaching node attributes to graph ...")
    nx.set_node_attributes(G, degree,           "degree")
    nx.set_node_attributes(G, pagerank,         "pagerank")
    nx.set_node_attributes(G, closeness,        "closeness")
    nx.set_node_attributes(G, betweenness,      "betweenness")
    nx.set_node_attributes(G, closeness_1000,   "closeness_1000")
    nx.set_node_attributes(G, betweenness_1000, "betweenness_1000")

    # ------------------------------------------------------------------
    # STEP 7  Convert to GeoDataFrame; compute edge-level metrics
    # ------------------------------------------------------------------
    print("\n[7/8] Converting to GeoDataFrame and computing edge metrics ...")
    gdf_nodes, gdf_edges = ox.graph_to_gdfs(G)
    gdf_edges = gdf_edges.reset_index().to_crs(epsg=32618)

    # each edge metric = mean of its two incident node values
    all_metrics = {
        "degree":           degree,
        "pagerank":         pagerank,
        "closeness":        closeness,
        "betweenness":      betweenness,
        "closeness_1000":   closeness_1000,
        "betweenness_1000": betweenness_1000,
    }
    for metric, mdict in all_metrics.items():
        u_arr = np.array([mdict.get(u, 0) for u in gdf_edges["u"]])
        v_arr = np.array([mdict.get(v, 0) for v in gdf_edges["v"]])
        gdf_edges[metric] = (u_arr + v_arr) / 2

    # ------------------------------------------------------------------
    # STEP 8  Export results
    # ------------------------------------------------------------------
    print(f"\n[8/8] Saving outputs ...")

    gdf_edges.to_file(OUT_GEOJSON, driver="GeoJSON")
    print(f"  GeoJSON   → {OUT_GEOJSON}")

    gdf_edges.to_file(OUT_SHP)
    print(f"  Shapefile → {OUT_SHP}")

    df_csv = gdf_edges.copy()
    df_csv["geometry_wkt"] = df_csv.geometry.to_wkt()
    df_csv.drop(columns="geometry").to_csv(OUT_CSV, index=False)
    print(f"  CSV       → {OUT_CSV}")

    # clean up temp pickle
    if os.path.exists(GRAPH_PKL):
        os.remove(GRAPH_PKL)

    elapsed = (time.time() - total_start) / 60
    print(f"\nDone. Total runtime: {elapsed:.1f} min")
    print(
        "\nNext step: run  02_covariates.py  — it will load\n"
        f"  {OUT_GEOJSON}"
    )


if __name__ == "__main__":
    main()
