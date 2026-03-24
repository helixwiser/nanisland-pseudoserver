# -*- coding: utf-8 -*-
import pandas as pd
import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
import community as community_louvain
import pickle
import time
import os


# === 打补丁开始 ===
def to_scipy_sparse_matrix(G, nodelist=None, dtype=None, weight='weight', format='lil'):
    M = nx.to_scipy_sparse_array(G, nodelist=nodelist, dtype=dtype, weight=weight, format=format)
    return M


nx.to_scipy_sparse_matrix = to_scipy_sparse_matrix


# 注意：不再导入 fa2，因为我们使用 kamada_kawai_layout


def log_step(message):
    print("Step: " + message)
    return time.time()


# === 创建输出目录 ===
output_path = "kk_layout/"
os.makedirs(output_path, exist_ok=True)  # 创建 kk_layout/ 文件夹（如果不存在）


start_time = time.time()
current_time = start_time
print("Pipeline started")




# --- Step 1: Load data ---
t0 = log_step("Loading network data...")
df = pd.read_csv('representation_protein_network.csv')
print("Loaded " + str(len(df)) + " edges")
t1 = time.time()
print("Time elapsed: " + str(t1 - t0))




## --- Step 2: Filter edges ---
#t0 = log_step("Filtering top 5000 edges by weight")
#weight_threshold = 0
#df = df[df['weight'] >= weight_threshold].nlargest(5000, 'weight')
#t1 = time.time()
#print("Time elapsed: " + str(t1 - t0))




# --- Step 3: Build full graph ---
t0 = log_step("Building full graph")
G_full = nx.Graph()
for _, row in df.iterrows():
    G_full.add_edge(row['family1'], row['family2'], weight=np.log(row['weight']))
print("Original network has " + str(G_full.number_of_nodes()) + " nodes and " + str(G_full.number_of_edges()) + " edges")
t1 = time.time()
print("Time elapsed: " + str(t1 - t0))




# --- Step 4: Compute MST skeleton ---
t0 = log_step("Computing maximum spanning tree for skeleton")
G_mst = nx.maximum_spanning_tree(G_full, weight='weight')
G_skeleton = G_mst.copy()
t1 = time.time()
print("Time elapsed: " + str(t1 - t0))
print("Skeleton network has " + str(G_skeleton.number_of_nodes()) + " nodes and " + str(G_skeleton.number_of_edges()) + " edges")




# --- Step 5: Compute centrality metrics ---
t0 = log_step("Computing degree, betweenness, pagerank")
degree_dict = dict(G_skeleton.degree())
betweenness = nx.betweenness_centrality(G_skeleton)
pagerank = nx.pagerank(G_skeleton, weight='weight')
t1 = time.time()
print("Time elapsed: " + str(t1 - t0))




# --- Step 6: Print top hubs ---
t0 = log_step("Printing top 10 hub nodes")
top_nodes = sorted(degree_dict.items(), key=lambda x: x[1], reverse=True)[:10]
print("Top 10 Hub Nodes in Skeleton:")
for i, (node, degree) in enumerate(top_nodes, 1):
    print("  " + str(i).rjust(2) + ". " + str(node) + ": " + str(degree) + " connections")
if len(degree_dict) >= 20:
    top20_degree_threshold = sorted(degree_dict.values(), reverse=True)[19]
else:
    top20_degree_threshold = 0
t1 = time.time()
print("Time elapsed: " + str(t1 - t0))




# --- Step 7: Set plot style ---
t0 = log_step("Setting matplotlib parameters")
plt.rcParams['font.size'] = 10
plt.rcParams['axes.facecolor'] = 'white'
plt.rcParams['savefig.facecolor'] = 'white'
plt.rcParams['savefig.bbox'] = 'tight'
plt.rcParams['figure.figsize'] = (12, 10)
t1 = time.time()
print("Time elapsed: " + str(t1 - t0))




# --- Step 8: Extract largest connected component ---
t0 = log_step("Extracting largest connected component for visualization")
if G_skeleton.number_of_nodes() > 20000:
    print("Network too large for full layout, using largest connected component...")
    G_plot = max(nx.connected_components(G_skeleton), default=[], key=len)
    G_vis = G_skeleton.subgraph(G_plot).copy()
else:
    G_vis = G_skeleton.copy()
t1 = time.time()
print("Time elapsed: " + str(t1 - t0))




# --- Step 9: Ensure connected for layout ---
t0 = log_step("Ensuring graph is connected for layout")
G_main = G_vis.copy()
if not nx.is_connected(G_main):
    G_main = G_main.subgraph(max(nx.connected_components(G_main), key=len)).copy()
t1 = time.time()
print("Time elapsed: " + str(t1 - t0))




# --- Step 10: Create dummy gene/protein maps ---
t0 = log_step("Creating dummy UID to gene/protein mapping")
uid_to_gene = dict(zip(list(G_main.nodes()), list(G_main.nodes())))
uid_to_protein = dict(zip(list(G_main.nodes()), list(G_main.nodes())))
t1 = time.time()
print("Time elapsed: " + str(t1 - t0))




# --- Step 11: Louvain community detection ---
t0 = log_step("Running Louvain community detection")
partition = community_louvain.best_partition(G_main, weight='weight', resolution = 0.5, random_state=42)
n_communities = len(set(partition.values()))
t1 = time.time()
print("Time elapsed: " + str(t1 - t0))




# --- Step 12: Compute node sizes ---
t0 = log_step("Computing node sizes from degree")
degree_dict = dict(G_main.degree())
max_degree = max(degree_dict.values()) if degree_dict else 1
node_size = [10 + (degree_dict[node] / max_degree) * 100 for node in G_main.nodes()]
t1 = time.time()
print("Time elapsed: " + str(t1 - t0))




# --- Step 13: Compute Kamada-Kawai layout (with weighted distances) ---
t0 = log_step("Computing Kamada-Kawai layout (weighted via shortest paths)")


if G_main.number_of_nodes() > 1000:
    print("Graph too large for Kamada-Kawai (>1000 nodes), falling back to spring layout.")
    pos = nx.spring_layout(G_main, weight='weight', seed=42, iterations=50)
else:
    try:
        # 计算加权最短路径长度作为距离
        dist = dict(nx.shortest_path_length(G_main, weight='weight'))
        dist_matrix = {}
        for u in G_main.nodes():
            dist_matrix[u] = {}
            for v in G_main.nodes():
                if u == v:
                    dist_matrix[u][v] = 0.0
                else:
                    try:
                        d = dist[u][v]
                        dist_matrix[u][v] = float(d)
                    except KeyError:
                        dist_matrix[u][v] = 10.0  # 不连通则设为大距离
        pos = nx.kamada_kawai_layout(G_main, dist=dist_matrix, weight=None)
    except Exception as e:
        print("Error in Kamada-Kawai layout: ", e)
        print("Falling back to spring layout.")
        pos = nx.spring_layout(G_main, weight='weight', seed=42, iterations=50)


t1 = time.time()
print("Time elapsed: " + str(t1 - t0))




# --- Step 14: Compute hub centralities ---
t0 = log_step("Computing hub metrics: degree, betweenness, closeness, pagerank")
degree_centrality = nx.degree_centrality(G_main)
betweenness = nx.betweenness_centrality(G_main, weight='weight')
closeness = nx.closeness_centrality(G_main)
pagerank = nx.pagerank(G_main, weight='weight')
t1 = time.time()
print("Time elapsed: " + str(t1 - t0))




# --- Step 15: Save hub node table ---
t0 = log_step("Saving hub node data to CSV")
hub_data = []
for node in G_main.nodes():
    gene_name = uid_to_gene.get(node, "")
    protein_name = uid_to_protein.get(node, "")
    hub_data.append({
        'Node': node,
        'Gene_Name': gene_name,
        'Protein_Name': protein_name,
        'Degree': str(degree_dict[node]),
        'Degree_Centrality': str(degree_centrality[node]),
        'Betweenness': str(betweenness[node]),
        'Closeness': str(closeness[node]),
        'PageRank': str(pagerank[node])
    })
hub_df = pd.DataFrame(hub_data)
hub_df['Degree_int'] = pd.to_numeric(hub_df['Degree'], errors='coerce')
hub_df = hub_df.sort_values('Degree_int', ascending=False).reset_index(drop=True)
hub_df.drop(columns=['Degree_int'], inplace=True)
hub_csv = os.path.join(output_path, "ppi_hub_nodes_annotated.csv")
hub_df.to_csv(hub_csv, index=False)
t1 = time.time()
print("Time elapsed: " + str(t1 - t0))




# --- Step 16: Build hub connection graph ---
t0 = log_step("Building hub connection graph for longest path")
top_hubs = hub_df['Node'].head(100).tolist()
hub_graph = nx.Graph()
for i, src in enumerate(top_hubs):
    for j, tgt in enumerate(top_hubs):
        if i < j and nx.has_path(G_main, src, tgt):
            try:
                path = nx.shortest_path(G_main, src, tgt)
                hub_graph.add_edge(src, tgt, weight=1.0 / len(path))
            except:
                pass
t1 = time.time()
print("Time elapsed: " + str(t1 - t0))




# --- Step 17: Find longest path among hubs ---
t0 = log_step("Finding longest path among top hubs")
visited = set()
path_order = []
current = top_hubs[0]
path_order.append(current)
visited.add(current)
while len(visited) < len(top_hubs):
    neighbors = [n for n in hub_graph.neighbors(current) if n not in visited]
    if not neighbors:
        break
    next_node = max(neighbors, key=lambda n: hub_graph.degree(n))
    path_order.append(next_node)
    visited.add(next_node)
    current = next_node
t1 = time.time()
print("Time elapsed: " + str(t1 - t0))




# --- Step 18: Reconstruct full path ---
t0 = log_step("Reconstructing full shortest-path sequence")
full_path_nodes = []
for i in range(len(path_order) - 1):
    src = path_order[i]
    tgt = path_order[i + 1]
    try:
        segment = nx.shortest_path(G_main, src, tgt)
        if i == 0:
            full_path_nodes.extend(segment)
        else:
            full_path_nodes.extend(segment[1:])
    except:
        continue
longest_path_str = " -> ".join(full_path_nodes)
full_hub_count = len(path_order)
t1 = time.time()
print("Time elapsed: " + str(t1 - t0))




# --- Step 19: Extract core subgraph ---
t0 = log_step("Extracting core subgraph along path")
core_nodes = set(n for n in longest_path_str.split(" -> "))
core_subgraph = G_main.subgraph(core_nodes).copy() if core_nodes else G_main
t1 = time.time()
print("Time elapsed: " + str(t1 - t0))




# --- Step 20a: Annotate core subgraph with UniProt details ---
t0 = log_step("Annotating core subgraph with UniProt details")


# 读取注释文件
annot_df = pd.read_csv('merged_details_unique.tsv', sep='\t', dtype=str).fillna("Unknown")
print(f"Loaded annotation for {len(annot_df)} proteins")


# 构建映射字典
uniprot_to_protein = pd.Series(annot_df['protein_name'].values, index=annot_df['uniprot_id']).to_dict()
uniprot_to_gene = pd.Series(annot_df['gene_name'].values, index=annot_df['uniprot_id']).to_dict()
uniprot_to_organism = pd.Series(annot_df['Organism'].values, index=annot_df['uniprot_id']).to_dict()


# 准备注释数据
core_annotation = []
for node in core_subgraph.nodes():
    protein_name = uniprot_to_protein.get(node, "Unknown")
    gene_name = uniprot_to_gene.get(node, "Unknown")
    organism = uniprot_to_organism.get(node, "Unknown")
    core_annotation.append({
        'UniProt_ID': node,
        'Gene_Name': gene_name,
        'Protein_Name': protein_name,
        'Organism': organism,
        'Degree_in_Main': degree_dict.get(node, 0),
        'Betweenness': betweenness.get(node, 0),
        'PageRank': pagerank.get(node, 0)
    })


# 转为 DataFrame 并排序
core_annot_df = pd.DataFrame(core_annotation)
core_annot_df = core_annot_df.sort_values('Degree_in_Main', ascending=False).reset_index(drop=True)


# 保存为 CSV
core_annot_csv = os.path.join(output_path, "core_path_subgraph_annotated.csv")
core_annot_df.to_csv(core_annot_csv, index=False)
print(f"Annotated {len(core_annot_df)} nodes in core path")
print(f"Top annotated nodes:")
print(core_annot_df[['UniProt_ID', 'Gene_Name', 'Protein_Name']].head(10))
t1 = time.time()
print("Time elapsed: " + str(t1 - t0))




# --- Step 20: Create 2x2 visualization ---
t0 = log_step("Creating 2x2 subplot visualization")
fig, axes = plt.subplots(2, 2, figsize=(16, 14))
axes = axes.ravel()


node_color_degree = [degree_dict[node] for node in G_main.nodes()]
cmap = plt.cm.plasma


node_color_comm = [partition[node] for node in G_main.nodes()]
cmap_comm = plt.cm.tab20 if n_communities <= 20 else plt.cm.tab20b
norm_comm = plt.Normalize(vmin=min(node_color_comm), vmax=max(node_color_comm))
colors_comm = [cmap_comm(norm_comm(c)) for c in node_color_comm]


# Plot 1: Layout + degree
nx.draw_networkx_nodes(G_main, pos, ax=axes[0], node_size=node_size, node_color=node_color_degree, cmap=cmap, alpha=0.8)
nx.draw_networkx_edges(G_main, pos, ax=axes[0], edge_color='gray', alpha=0.3, width=0.5)
axes[0].set_title("1. Kamada-Kawai Layout\nColor: degree", fontsize=14)
axes[0].axis('off')


# Plot 2: Layout + communities
nx.draw_networkx_nodes(G_main, pos, ax=axes[1], node_size=node_size, node_color=colors_comm, alpha=0.8)
nx.draw_networkx_edges(G_main, pos, ax=axes[1], edge_color='gray', alpha=0.3, width=0.5)
axes[1].set_title("2. Kamada-Kawai + Louvain\n" + "Communities: " + str(n_communities), fontsize=14)
axes[1].axis('off')


# Plot 3: Longest path
node_colors_path = ['red' if node in path_order else 'lightblue' for node in G_main.nodes()]
nx.draw_networkx_nodes(G_main, pos, ax=axes[2], node_size=node_size, node_color=node_colors_path, alpha=0.8)
nx.draw_networkx_edges(G_main, pos, ax=axes[2], edge_color='lightgray', alpha=0.3, width=0.5)
nx.draw_networkx_edges(core_subgraph, pos, ax=axes[2], edge_color='red', alpha=0.8, width=1.5)
axes[2].set_title("3. Longest Path Among Top Hubs\n" + "Hub path length: " + str(full_hub_count), fontsize=14)
axes[2].axis('off')


# Plot 4: Core subgraph with annotated labels
nx.draw_networkx_nodes(core_subgraph, pos, ax=axes[3], node_size=50, node_color='orange', alpha=0.9)
nx.draw_networkx_edges(core_subgraph, pos, ax=axes[3], edge_color='red', width=2.0, alpha=0.9)


labels = {}
for node in core_subgraph.nodes():
    gene_name = uniprot_to_gene.get(node, "Unknown")
    if isinstance(gene_name, str) and gene_name != "Unknown" and gene_name.strip() != "" and len(gene_name) <= 15:
        labels[node] = gene_name


nx.draw_networkx_labels(
    core_subgraph,
    pos,
    labels,
    ax=axes[3],
    font_size=9,
    font_weight='bold',
    font_color='darkblue'
)
axes[3].set_title("4. Core Path Subgraph\n" + "Nodes: " + str(core_subgraph.number_of_nodes()) +
                  " Edges: " + str(core_subgraph.number_of_edges()), fontsize=14)
axes[3].axis('off')


# Colorbar
sm_degree = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(vmin=min(node_color_degree), vmax=max(node_color_degree)))
sm_degree.set_array(node_color_degree)
cbar1 = plt.colorbar(sm_degree, ax=axes[1], shrink=0.6, pad=0.02, location='right')
cbar1.set_label('Node Degree', rotation=270, labelpad=15)


fig.suptitle("Protein Interaction Network Analysis (Kamada-Kawai Layout)", fontsize=18, y=0.95)
output_file = os.path.join(output_path, "ppi_network_four_panel.png")
plt.savefig(output_file, dpi=300, bbox_inches='tight')
plt.close()
t1 = time.time()
print("Time elapsed: " + str(t1 - t0))




# --- Step X: Interactive PyVis visualization of core subgraph ---
try:
    from pyvis.network import Network
    import json


    t0 = log_step("Creating interactive PyVis visualization of core subgraph")


    pyvis_net = Network(height="800px", width="100%", bgcolor="#ffffff", font_color="#000000")
    pyvis_net.set_options("""
    var options = {
      "physics": {
        "enabled": true,
        "repulsion": {
          "nodeDistance": 120,
          "centralGravity": 0.1,
          "springLength": 150,
          "springStrength": 0.05,
          "damping": 0.2
        },
        "maxVelocity": 50,
        "minVelocity": 0.75,
        "solver": "repulsion"
      },
      "interaction": {
        "hover": true,
        "tooltipDelay": 100
      }
    }
    """)


    core_degrees = dict(core_subgraph.degree())
    max_degree = max(core_degrees.values()) if core_degrees else 1
    core_partition = {node: partition.get(node, -1) for node in core_subgraph.nodes()}
    communities = set(core_partition.values())
    cmap = plt.cm.tab20 if len(communities) <= 20 else plt.cm.Set3
    color_map = {}
    for i, comm in enumerate(sorted(communities)):
        color = tuple(int(x * 255) for x in cmap(i % cmap.N)[:3])
        color_map[comm] = f"rgb{color}"


    for node in core_subgraph.nodes():
        gene_name = uniprot_to_gene.get(node, "Unknown")
        protein_name = uniprot_to_protein.get(node, "Unknown")
        organism = uniprot_to_organism.get(node, "Unknown")
        degree = core_degrees[node]
        comm_id = core_partition[node]
        size = 10 + (degree / max_degree) * 30
        color = color_map.get(comm_id, "lightblue")
        title = f"""
        UniProt ID: {node}
        Gene Name: {gene_name}
        Protein Name: {protein_name}
        Organism: {organism}
        Degree: {degree}
        Community: {comm_id}
        """
        label = gene_name if gene_name != "Unknown" and len(gene_name) <= 15 else ""
        pyvis_net.add_node(node, label=label, title=title, size=int(size), color=color, shape="dot")


    for u, v, data in core_subgraph.edges(data=True):
        weight = data.get("weight", 1.0)
        width = 0.5 + (weight / max(df['weight']) * 5) if 'weight' in df.columns else 1.0
        pyvis_net.add_edge(u, v, width=width, color="lightgray")


    pyvis_output = os.path.join(output_path, "core_subgraph_interactive.html")
    pyvis_net.write_html(pyvis_output)
    t1 = time.time()
    print(f"Interactive visualization saved to: {pyvis_output}")
    print("Time elapsed: " + str(t1 - t0))


except Exception as e:
    print("Error creating PyVis visualization:")
    print(e)




# --- Step 21: Save hub path ---
t0 = log_step("Saving longest hub path to text file")
path_txt = os.path.join(output_path, "ppi_longest_hub_path.txt")
with open(path_txt, "w") as f:
    f.write("Longest connecting path among top hubs:\n")
    f.write(longest_path_str + "\n")
    f.write("Number of hubs connected: " + str(full_hub_count) + "\n")
t1 = time.time()
print("Time elapsed: " + str(t1 - t0))




# --- Step 22: Save skeleton graph ---
t0 = log_step("Saving skeleton graph to pickle")
graph_pickle = os.path.join(output_path, "protein_network_skeleton.pkl")
with open(graph_pickle, 'wb') as f:
    pickle.dump(G_skeleton, f)
t1 = time.time()
print("Time elapsed: " + str(t1 - t0))



#
## --- Step 23: Export skeleton edges ---
#t0 = log_step("Exporting skeleton edges to CSV")
#edges_data = []
#for u, v, data in G_skeleton.edges(data=True):
#    edges_data.append({"source": str(u), "target": str(v), "weight": data["weight"]})
#df_skeleton_edges = pd.DataFrame(edges_data)
#skeleton_csv = os.path.join(output_path, "skeleton_edges.csv")
#df_skeleton_edges.to_csv(skeleton_csv, index=False)
#t1 = time.time()
#print("Time elapsed: " + str(t1 - t0))

# --- Step 23: Export skeleton edges with node communities ---
t0 = log_step("Exporting skeleton edges with community labels")
# Ensure partition includes all skeleton nodes
partition_full = {node: partition.get(node, -1) for node in G_skeleton.nodes()}


edges_data = []
for u, v, data in G_skeleton.edges(data=True):
    edges_data.append({
        "source": str(u),
        "target": str(v),
        "weight": str(data["weight"]),
        "source_community": str(partition_full[u]),
        "target_community": str(partition_full[v])
    })


df_skeleton_edges = pd.DataFrame(edges_data)
skeleton_csv = os.path.join(output_path, "skeleton_edges.csv")
df_skeleton_edges.to_csv(skeleton_csv, index=False)
t1 = time.time()
print("Time elapsed: " + str(t1 - t0))

# --- Final Summary: Print + Save to Log File ---
total_time = time.time() - start_time


log_lines = [
    "=" * 60,
    "All tasks completed",
    "Total runtime: " + str(round(total_time, 2)) + " seconds",
    "=" * 60,
    "Output files saved to folder: " + output_path,
    "Output files:",
    "  - Four-panel plot: " + os.path.join(output_path, "ppi_network_four_panel.png"),
    "  - Hub node table: " + os.path.join(output_path, "ppi_hub_nodes_annotated.csv"),
    "  - Longest hub path: " + os.path.join(output_path, "ppi_longest_hub_path.txt"),
    "  - Skeleton graph (pickle): " + os.path.join(output_path, "protein_network_skeleton.pkl"),
    "  - Skeleton edge list: " + os.path.join(output_path, "skeleton_edges.csv"),
    "Load skeleton with:",
    "    import pickle",
    f"    with open('{os.path.join(output_path, 'protein_network_skeleton.pkl')}', 'rb') as f:",
    "        G = pickle.load(f)",
    "=" * 60
]


# 打印到控制台
for line in log_lines:
    print(line)


# 写入日志文件
log_file_path = os.path.join(output_path, "pipeline_log.txt")
with open(log_file_path, "w", encoding="utf-8") as log_file:
    log_file.write("\n".join(log_lines))


print(f"Pipeline log saved to: {log_file_path}")