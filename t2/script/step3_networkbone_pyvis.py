# -*- coding: utf-8 -*-
import pandas as pd
import networkx as nx
from pyvis.network import Network
import numpy as np
from scipy.linalg import expm




print("Loading network data...")
df = pd.read_csv('representation_protein_network.csv')
print("Loaded " + str(len(df)) + " edges")
##########################
#weight_threshold = 0  # Can be adjusted
#df = df[df['weight'] >= weight_threshold].nlargest(5000, 'weight')
####################

# Build weighted undirected graph
G_full = nx.Graph()
for _, row in df.iterrows():
    G_full.add_edge(row['family1'], row['family2'], weight=row['weight'])
print("Original network has " + str(G_full.number_of_nodes()) + " nodes and " + str(G_full.number_of_edges()) + " edges")




## Function to extract communicability backbone
#def extract_communicability_backbone(G, alpha=0.05):
#    A = nx.to_numpy_array(G, weight='weight')
#    n = A.shape[0]
#    nodes = list(G.nodes())
#
#
#
#
#    # Compute communicability matrix: C = exp(A)
#    C = expm(A)
#
#
#
#
#    # Compute expected values based on degree product: E_ij = (k_i * k_j) / (2m)
#    k = np.array([G.degree(n, weight='weight') for n in nodes])
#    m = G.size(weight='weight')
#    if m == 0:
#        m = 1
#    E = np.outer(k, k) / (2 * m)
#
#
#
#
#    # Create backbone graph
#    G_backbone = nx.Graph()
#    G_backbone.add_nodes_from(G.nodes())
#
#
#
#
#    for i in range(n):
#        for j in range(i + 1, n):
#            if G.has_edge(nodes[i], nodes[j]):
#                if C[i, j] > alpha * E[i, j]:
#                    weight = G[nodes[i]][nodes[j]]['weight']
#                    G_backbone.add_edge(nodes[i], nodes[j], weight=weight)
#
#
#
#
#    return G_backbone
#
#
#
#
## Extract skeleton
#print("Extracting communicability backbone...")
#G_skeleton = extract_communicability_backbone(G_full, alpha=0.05)
#
#
#
#
## Fallback to MST if too sparse
#if G_skeleton.number_of_edges() < G_skeleton.number_of_nodes():
#    print("Skeleton too sparse, using Maximum Spanning Tree instead")
#    G_mst = nx.maximum_spanning_tree(G_full, weight='weight')
#    G_skeleton = G_mst.copy()

G_mst = nx.maximum_spanning_tree(G_full, weight='weight')
G_skeleton = G_mst.copy()


print("Skeleton network has " + str(G_skeleton.number_of_nodes()) + " nodes and " + str(G_skeleton.number_of_edges()) + " edges")




# Compute node metrics
degree_dict = dict(G_skeleton.degree())
betweenness = nx.betweenness_centrality(G_skeleton)
pagerank = nx.pagerank(G_skeleton, weight='weight')




top_nodes = sorted(degree_dict.items(), key=lambda x: x[1], reverse=True)[:10]
print("Top 10 Hub Nodes in Skeleton:")
for i, (node, degree) in enumerate(top_nodes, 1):
    print("  " + str(i).rjust(2) + ". " + str(node) + ": " + str(degree) + " connections")




# Threshold for labeling top nodes
if len(degree_dict) >= 20:
    top20_degree_threshold = sorted(degree_dict.values(), reverse=True)[19]
else:
    top20_degree_threshold = 0




# Create PyVis network for skeleton
print("Creating interactive skeleton visualization...")
net = Network(
    height='900px',
    width='100%',
    bgcolor='#ffffff',
    font_color='#000000',
    notebook=False,
    cdn_resources='in_line'
)




# Physics configuration (plain string)
net.set_options("""
{
  "physics": {
    "forceAtlas2Based": {
      "gravitationalConstant": -50,
      "centralGravity": 0.01,
      "springLength": 100,
      "springConstant": 0.08,
      "damping": 0.8,
      "avoidOverlap": 0.5
    },
    "maxVelocity": 2,
    "minVelocity": 0.01,
    "solver": "forceAtlas2Based",
    "timestep": 0.35,
    "stabilization": {
      "enabled": true,
      "iterations": 100,
      "updateInterval": 50,
      "onlyDynamicEdges": false,
      "fit": true
    }
  },
  "interaction": {
    "hover": true,
    "tooltipDelay": 100,
    "navigationButtons": true,
    "keyboard": true
  },
  "nodes": {
    "font": {
      "size": 12,
      "face": "arial"
    }
  },
  "edges": {
    "smooth": {
      "type": "continuous"
    }
  }
}
""")




# Normalize for visual scaling
max_degree = max(degree_dict.values()) if degree_dict else 1
edge_weights = [G_skeleton[u][v]['weight'] for u, v in G_skeleton.edges()]
max_weight = max(edge_weights) if edge_weights else 1
min_weight = min(edge_weights) if edge_weights else 0




# Add nodes
print("Adding skeleton nodes...")
for node in G_skeleton.nodes():
    degree = degree_dict[node]
    bet = betweenness[node]
    pr = pagerank[node]
    node_size = 12 + (degree / max_degree) * 35
    color_intensity = int((degree / max_degree) * 255)
    r = 200 - color_intensity
    g = color_intensity // 2
    b = color_intensity
    node_color = f"#{r:02x}{g:02x}{b:02x}"




    title = "<b>" + str(node) + "</b><br>Degree: " + str(degree) + \
            "<br>Betweenness: " + str(bet) + \
            "<br>PageRank: " + str(pr)
    label = str(node) if degree >= top20_degree_threshold else ""




    net.add_node(
        node,
        label=label,
        title=title,
        size=node_size,
        color=node_color,
        borderWidth=2,
        borderWidthSelected=4
    )




# Add edges
print("Adding skeleton edges...")
for u, v, data in G_skeleton.edges(data=True):
    weight = data['weight']
    edge_width = 1.0 + (weight - min_weight) / (max_weight - min_weight + 1e-8) * 4
    opacity = 0.5 + (weight - min_weight) / (max_weight - min_weight + 1e-8) * 0.4
    edge_color = "rgba(100, 100, 100, " + str(opacity) + ")"




    net.add_edge(
        u, v,
        width=edge_width,
        title="Weight: " + str(round(weight, 3)),
        color=edge_color
    )




# Save to file
output_file = "protein_network_skeleton.html"
print("Generating " + output_file + "...")
net.save_graph(output_file)




# Read HTML content
with open(output_file, 'r', encoding='utf-8') as f:
    html_content = f.read()




# Build plain ASCII header with stats
custom_header = """
<style>
    body { margin: 0; padding: 0; font-family: Arial, sans-serif; }
    #header { background: linear-gradient(135deg, #667eea 0%, #764ba2 100%); color: white; padding: 20px; text-align: center; box-shadow: 0 2px 10px rgba(0,0,0,0.1); }
    #header h1 { margin: 0; font-size: 28px; }
    #header p { margin: 5px 0 0 0; font-size: 14px; opacity: 0.9; }
    #stats-panel { position: fixed; top: 120px; left: 20px; width: 250px; background: white; padding: 15px; border-radius: 8px; box-shadow: 0 2px 10px rgba(0,0,0,0.2); z-index: 1000; font-size: 13px; }
    .stat-row { display: flex; justify-content: space-between; margin: 5px 0; padding: 5px 0; border-bottom: 1px solid #eee; }
    .stat-label { font-weight: bold; color: #555; }
    .stat-value { color: #667eea; }
</style>
<div id="header">
    <h1>Protein Family Network - Core Skeleton</h1>
    <p>Interactive visualization of the core backbone extracted via communicability analysis</p>
</div>
<div id="stats-panel">
    <h3>Skeleton Statistics</h3>
    <div class="stat-row"><span class="stat-label">Nodes:</span><span class="stat-value">""" + str(G_skeleton.number_of_nodes()) + """</span></div>
    <div class="stat-row"><span class="stat-label">Edges:</span><span class="stat-value">""" + str(G_skeleton.number_of_edges()) + """</span></div>
    <div class="stat-row"><span class="stat-label">Sparsity Reduction:</span><span class="stat-value">""" + str(round(1 - (G_skeleton.number_of_edges() / (G_full.number_of_edges() + 1e-8)), 2)) + """</span></div>
    <div class="stat-row"><span class="stat-label">Density:</span><span class="stat-value">""" + str(round(nx.density(G_skeleton), 4)) + """</span></div>
    <div class="stat-row"><span class="stat-label">Components:</span><span class="stat-value">""" + str(nx.number_connected_components(G_skeleton)) + """</span></div>
    <div style="margin-top: 15px; padding-top: 10px; border-top: 2px solid #eee;">
        <div style="font-weight: bold; color: #764ba2; margin-bottom: 8px;">Top Hubs:</div>
"""


for i, (node, degree) in enumerate(top_nodes[:5], 1):
    custom_header += "        <div style=\"margin: 3px 0; font-size: 11px;\">#" + str(i) + " " + str(node) + " (" + str(degree) + ")</div>\n"


custom_header += """
    </div>
</div>
"""




# Inject header
html_content = html_content.replace("<body>", "<body>\n" + custom_header)




# Write back
with open(output_file, 'w', encoding='utf-8') as f:
    f.write(html_content)


# Save the skeleton graph object (for later use in Python)
graph_output_file = "protein_network_skeleton.pkl"
print("Saving skeleton graph object to " + graph_output_file + "...")
import pickle
with open(graph_output_file, 'wb') as f:
    pickle.dump(G_skeleton, f)
print("Skeleton graph saved: " + graph_output_file)




# Export skeleton edges to CSV
edges_data = []
for u, v, data in G_skeleton.edges(data=True):
    row = {
        "source": str(u),
        "target": str(v),
        "weight": data["weight"]
    }
    edges_data.append(row)




df_skeleton_edges = pd.DataFrame(edges_data)
skeleton_csv_file = "skeleton_edges.csv"
print("Exporting skeleton edge list to " + skeleton_csv_file + "...")
df_skeleton_edges.to_csv(skeleton_csv_file, index=False)
print("Skeleton edge table saved: " + skeleton_csv_file)






# Final output
print("=" * 60)
print("Interactive skeleton network visualization created")
print("=" * 60)
print("File: " + output_file)
print("Nodes: " + str(G_skeleton.number_of_nodes()))
print("Edges: " + str(G_skeleton.number_of_edges()))
print("Sparsity Reduction: " + str(round(1 - (G_skeleton.number_of_edges() / (G_full.number_of_edges() + 1e-8)), 2)))
print("Density: " + str(round(nx.density(G_skeleton), 4)))
print("Open '" + output_file + "' in a web browser to view the core structure")
print("This skeleton highlights central pathways in the network")
print("=" * 60)
# Final summary
print("=" * 60)
print("All tasks completed")
print("=" * 60)
print("Output files:")
print("  - Interactive visualization: " + output_file)
print("  - Skeleton graph (pickle): " + graph_output_file)
print("  - Skeleton edge table (CSV): " + skeleton_csv_file)
print("You can load the pickle file later with:")
print("    import pickle")
print("    with open('protein_network_skeleton.pkl', 'rb') as f:")
print("        G = pickle.load(f)")
print("=" * 60)