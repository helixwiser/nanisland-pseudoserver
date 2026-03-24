# -*- coding: utf-8 -*-
import pandas as pd
import networkx as nx
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import community as community_louvain
import requests
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
import os


# === Step 1: PPI Network Analysis ===
print("Loading network data...")
df = pd.read_csv('representation_protein_network.csv')


weight_threshold = 0
df = df[df['weight'] >= weight_threshold].nlargest(5000, 'weight')


G = nx.Graph()
for _, row in df.iterrows():
    G.add_edge(row['family1'], row['family2'], weight=row.get('weight', 1))
print("Network built.")


# Extract UniProt IDs from nodes
uniprot_ids = list(G.nodes())
print(f"Found {len(uniprot_ids)} unique UniProt IDs for annotation.")


# Save IDs to temporary file for UniProt query
input_id_file = "temp_uniprot_ids.txt"
with open(input_id_file, "w") as f:
    for uid in uniprot_ids:
        f.write(str(uid).strip() + "\n")


# === Step 2: Fetch Protein and Gene Names from UniProt ===
def fetch_single(uniprot_id):
    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}"
    try:
        response = requests.get(url, timeout=5)
        if response.status_code != 200:
            return [uniprot_id, "", "", ""]
        data = response.json()


        protein_name = ""
        if "proteinDescription" in data:
            rec_name = data["proteinDescription"].get("recommendedName")
            if rec_name and "fullName" in rec_name:
                protein_name = rec_name["fullName"].get("value", "")


        gene_name = ""
        if data.get("genes") and len(data["genes"]) > 0:
            gene = data["genes"][0]
            gene_name = gene.get("geneName", {}).get("value", "")


        org = data.get("organism", {})
        scientific = org.get("scientificName", "")
        common = org.get("commonName", "")
        taxon_id = str(org.get("taxonId", ""))
        if scientific and common:
            organism = scientific + " (" + common + ")"
        else:
            organism = scientific or "Unknown"
        if taxon_id:
            organism = organism + " [TaxID: " + taxon_id + "]"


        return [uniprot_id, protein_name, gene_name, organism]
    except Exception as e:
        return [uniprot_id, "", "", ""]


def write_tsv(data, output_path):
    with open(output_path, "w") as f:
        f.write("uniprot_id\tprotein_name\tgene_name\tOrganism\n")
        for row in data:
            f.write("\t".join(str(x) for x in row) + "\n")


print("Querying UniProt for protein and gene names...")
start_fetch = time.time()
all_data = []
max_workers = 10


with ThreadPoolExecutor(max_workers=max_workers) as executor:
    future_to_id = {executor.submit(fetch_single, uid): uid for uid in uniprot_ids}
    for future in as_completed(future_to_id):
        try:
            result = future.result()
            all_data.append(result)
        except Exception as e:
            uid = future_to_id[future]
            all_data.append([uid, "", "", ""])


output_tsv = "uniprot_detailed_annotations.tsv"
write_tsv(all_data, output_tsv)
print("UniProt annotation completed: " + output_tsv)


# Create mapping dictionary
annotation_df = pd.read_csv(output_tsv, sep="\t", dtype=str)
annotation_df.fillna("", inplace=True)
uid_to_gene = dict(zip(annotation_df['uniprot_id'], annotation_df['gene_name']))
uid_to_protein = dict(zip(annotation_df['uniprot_id'], annotation_df['protein_name']))


# Clean up temporary file
if os.path.exists(input_id_file):
    os.remove(input_id_file)


# === Step 3: Continue Network Analysis ===
largest_cc = max(nx.connected_components(G), key=len)
G_main = G.subgraph(largest_cc).copy()


# Topological metrics
n_nodes = str(G.number_of_nodes())
n_edges = str(G.number_of_edges())
n_components = str(nx.number_connected_components(G))
main_size = str(len(largest_cc))


avg_clust_val = nx.average_clustering(G)
try:
    avg_path_val = nx.average_shortest_path_length(G_main) if len(G_main) > 1 else float('nan')
except Exception:
    avg_path_val = float('nan')


partition = community_louvain.best_partition(G_main)
modularity_val = community_louvain.modularity(partition, G_main)


avg_clust = str(avg_clust_val)
avg_path = str(avg_path_val) if not pd.isna(avg_path_val) else "N/A"
modularity = str(modularity_val)


# Centrality measures
degree_centrality = nx.degree_centrality(G_main)
betweenness = nx.betweenness_centrality(G_main)
closeness = nx.closeness_centrality(G_main)
pagerank = nx.pagerank(G_main, weight='weight')


hub_data = []
for node in G_main.nodes():
    gene_name = uid_to_gene.get(node, "")
    protein_name = uid_to_protein.get(node, "")
    hub_data.append({
        'Node': node,
        'Gene_Name': gene_name,
        'Protein_Name': protein_name,
        'Degree': str(G_main.degree(node)),
        'Degree_Centrality': str(degree_centrality[node]),
        'Betweenness': str(betweenness[node]),
        'Closeness': str(closeness[node]),
        'PageRank': str(pagerank[node])
    })


hub_df = pd.DataFrame(hub_data)
hub_df['Degree_int'] = pd.to_numeric(hub_df['Degree'], errors='coerce')
hub_df = hub_df.sort_values('Degree_int', ascending=False).reset_index(drop=True)
hub_df.drop(columns=['Degree_int'], inplace=True)
hub_df.to_csv('ppi_hub_nodes_annotated.csv', index=False)


# Core path: Longest connecting path among top 10 hubs
top_hubs = hub_df['Node'].head(10).tolist()


hub_graph = nx.Graph()
for i, src in enumerate(top_hubs):
    for j, tgt in enumerate(top_hubs):
        if i < j and nx.has_path(G_main, src, tgt):
            try:
                path = nx.shortest_path(G_main, src, tgt)
                hub_graph.add_edge(src, tgt, weight=1.0 / len(path))
            except:
                pass


if len(hub_graph.edges) == 0:
    longest_path_str = " -> ".join(top_hubs[:2]) + " (no full path found)"
    full_hub_count = 2
else:
    from collections import deque
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


# === Step 4: Build Core Subgraph ===
core_nodes = set(n for path in [longest_path_str.split(" -> ")] for n in path)
core_subgraph = G_main.subgraph(core_nodes).copy() if core_nodes else G_main


# Create node labels: Gene (Protein) or fallback to UniProt ID
def make_node_label(uid):
    gene = uid_to_gene.get(uid, "").strip()
    prot = uid_to_protein.get(uid, "").strip()
    if gene and prot:
        return gene + "\n(" + prot + ")"
    elif gene:
        return gene
    elif prot:
        return prot
    else:
        return uid  # fallback


# Generate PDF report
pdf_path = "PPI_Network_Key_Structure_Report_annotated.pdf"
with PdfPages(pdf_path) as pdf:
    # Page 1: Topology stats
    fig1, ax1 = plt.subplots(figsize=(8, 6))
    ax1.axis('off')
    stats_text = (
        "PPI Network Key Structure Analysis Report\n" +
        "="*60 + "\n\n"
        "Basic Topology\n"
        "  Nodes:                 " + n_nodes + "\n"
        "  Edges:                 " + n_edges + "\n"
        "  Connected Components:  " + n_components + "\n"
        "  Largest Component:     " + main_size + "\n"
        "  Avg Clustering:        " + avg_clust + "\n"
        "  Avg Shortest Path:     " + avg_path + "\n"
        "  Modularity:            " + modularity + "\n\n"
        "Top 5 Hub Nodes by Degree:\n"
    )
    for i in range(min(5, len(hub_df))):
        node = hub_df.iloc[i]['Node']
        deg = hub_df.iloc[i]['Degree']
        gene = hub_df.iloc[i]['Gene_Name'] or node
        stats_text += "  " + str(i+1) + ". " + gene + " (" + node + ") [degree=" + deg + "]\n"
    ax1.text(0.02, 0.98, stats_text, fontsize=11, va='top', family='monospace')
    pdf.savefig(fig1, bbox_inches='tight')
    plt.close(fig1)


    # Page 2: Longest path among hubs
    fig2, ax2 = plt.subplots(figsize=(9, 6))
    ax2.axis('off')
    path_text = "Longest Connecting Path Among Top 10 Hubs\n"
    path_text += "="*60 + "\n\n"
    path_text += "Hub coverage: " + str(full_hub_count) + " out of 10 hubs\n\n"
    path_text += "Full path:\n"
    if len(longest_path_str) > 500:
        import textwrap
        lines = textwrap.wrap(longest_path_str, width=80)
        for line in lines:
            path_text += line + "\n"
    else:
        path_text += longest_path_str + "\n"
    ax2.text(0.02, 0.98, path_text, fontsize=10, va='top', family='monospace', linespacing=1.3)
    pdf.savefig(fig2, bbox_inches='tight')
    plt.close(fig2)


    # Page 3: Core Network Backbone with gene/protein labels
    fig3, ax3 = plt.subplots(figsize=(12, 9))
    pos = nx.spring_layout(core_subgraph, k=0.5, seed=42, iterations=30)


    # Node colors and sizes
    colors = [partition[n] for n in core_subgraph.nodes()]
    sizes = [G_main.degree(n) * 20 for n in core_subgraph.nodes()]


    nx.draw_networkx_nodes(core_subgraph, pos, node_color=colors, cmap='tab20', node_size=sizes, alpha=0.85)
    nx.draw_networkx_edges(core_subgraph, pos, alpha=0.3, edge_color='gray')


    # Label with gene (protein) or fallback
    labels = {}
    for node in core_subgraph.nodes():
        labels[node] = make_node_label(node)


    nx.draw_networkx_labels(core_subgraph, pos, labels, font_size=9, font_color='darkred', bbox=dict(boxstyle="round,pad=0.2", facecolor="white", alpha=0.7))


    ax3.set_title("Core Network Backbone (Gene and Protein Names)", fontsize=14)
    ax3.text(0.01, 0.01,
             "Node size: degree, Color: community, Labels: Gene (Protein)",
             transform=ax3.transAxes, fontsize=9, color='black')
    ax3.axis('off')
    pdf.savefig(fig3, bbox_inches='tight')
    plt.close(fig3)


# Final output
elapsed = time.time() - start_fetch
print("Final report generated: " + pdf_path)
print("Annotated hub data saved: ppi_hub_nodes_annotated.csv")
print("Total annotation time: " + str(elapsed) + " seconds")