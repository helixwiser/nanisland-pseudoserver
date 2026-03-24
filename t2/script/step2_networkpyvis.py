# -*- coding: utf-8 -*-
import pandas as pd
import networkx as nx
from pyvis.network import Network


print("Loading network data...")
df = pd.read_csv('representation_protein_network.csv')
print(f"Loaded {len(df)} edges")


# Filter for top 500 edges by weight
weight_threshold = 0  # Can be adjusted
df_filtered = df[df['weight'] >= weight_threshold].nlargest(500, 'weight')
print(f"Visualizing {len(df_filtered)} edges with weight >= {weight_threshold}")


# Create NetworkX graph
G = nx.Graph()
for _, row in df_filtered.iterrows():
    G.add_edge(row['family1'], row['family2'], weight=row['weight'])


print(f"Network has {G.number_of_nodes()} nodes and {G.number_of_edges()} edges")


# Calculate node metrics
degree_dict = dict(G.degree())
betweenness = nx.betweenness_centrality(G)
# Note: Using PageRank on undirected graph (symmetric relationships)
pagerank = nx.pagerank(G, weight='weight')


print("Calculating network statistics...")
top_nodes = sorted(degree_dict.items(), key=lambda x: x[1], reverse=True)[:10]


print("\nTop 10 Hub Nodes:")
for i, (node, degree) in enumerate(top_nodes, 1):
    print(f"  {i:2d}. {node}: {degree} connections")


# Precompute threshold for top 20 nodes (if exist)
if len(degree_dict) >= 20:
    top20_degree_threshold = sorted(degree_dict.values(), reverse=True)[19]
else:
    top20_degree_threshold = 0


# Create PyVis network
print("\nCreating interactive visualization...")
net = Network(
    height='900px',
    width='100%',
    bgcolor='#ffffff',
    font_color='#000000',
    notebook=False,
    cdn_resources='in_line'  # Ensure offline compatibility
)


# Configure physics
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


# Normalize metrics
max_degree = max(degree_dict.values()) if degree_dict else 1
max_weight = df_filtered['weight'].max()
min_weight = df_filtered['weight'].min()


# Add nodes
print("Adding nodes...")
for node in G.nodes():
    degree = degree_dict[node]
    bet = betweenness[node]
    pr = pagerank[node]


    node_size = 10 + (degree / max_degree) * 40
    color_intensity = int((degree / max_degree) * 255)
    node_color = f'#{255-color_intensity:02x}{color_intensity//2:02x}{color_intensity:02x}'


    connections = list(G.neighbors(node))
    top_connections = sorted(
        [(n, G[node][n]['weight']) for n in connections],
        key=lambda x: x[1],
        reverse=True
    )[:5]


    title = f"{node}  "
    title += f"Degree: {degree}"
#    title += f"Betweenness: {bet:.4f}<br>"
#    title += f"PageRank: {pr:.4f}<br>"
#    title += "<br><b>Top Connections:</b><br>"
#    for conn, weight in top_connections:
#        title += f"  → {conn} (weight: {weight})<br>"


    label = node if degree >= top20_degree_threshold else ""


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
print("Adding edges...")
for _, row in df_filtered.iterrows():
    source = row['family1']
    target = row['family2']
    weight = row['weight']


    edge_width = 0.5 + (weight - min_weight) / (max_weight - min_weight) * 5
    opacity = 0.3 + (weight - min_weight) / (max_weight - min_weight) * 0.5
    edge_color = f'rgba(150, 150, 150, {opacity})'


    net.add_edge(
        source,
        target,
        width=edge_width,           # ← Now correctly controls line width
        title=f"Weight: {weight}",
        color=edge_color
    )


# Save and enhance HTML
output_file = 'protein_network_pyvis.html'
print(f"\nGenerating {output_file}...")
net.save_graph(output_file)


with open(output_file, 'r', encoding='utf-8') as f:
    html_content = f.read()


# Build custom header with proper structure
custom_header = f"""
<style>
    body {{ margin: 0; padding: 0; font-family: Arial, sans-serif; }}
    #header {{ background: linear-gradient(135deg, #667eea 0%, #764ba2 100%); color: white; padding: 20px; text-align: center; box-shadow: 0 2px 10px rgba(0,0,0,0.1); }}
    #header h1 {{ margin: 0; font-size: 28px; }}
    #header p {{ margin: 5px 0 0 0; font-size: 14px; opacity: 0.9; }}
    #info-panel, #stats-panel, .legend {{
        position: fixed; background: white; padding: 15px; border-radius: 8px;
        box-shadow: 0 2px 10px rgba(0,0,0,0.2); z-index: 1000;
    }}
    #info-panel {{ top: 120px; right: 20px; max-width: 300px; }}
    #stats-panel {{ top: 120px; left: 20px; max-width: 250px; font-size: 13px; }}
    .stat-row {{ display: flex; justify-content: space-between; margin: 5px 0; padding: 5px 0; border-bottom: 1px solid #eee; }}
    .stat-label {{ font-weight: bold; color: #555; }}
    .stat-value {{ color: #667eea; }}
    .legend {{ bottom: 20px; left: 20px; font-size: 12px; }}
    .legend-color {{ width: 20px; height: 20px; border-radius: 50%; margin-right: 10px; border: 2px solid white; }}
</style>
<div id="header">
    <h1>Protein Family Network Visualization</h1>
    <p>Interactive network showing protein family interactions | Node size and color represent connection degree</p>
</div>
<div id="stats-panel">
    <h3>Network Statistics</h3>
    <div class="stat-row"><span class="stat-label">Nodes:</span><span class="stat-value">{G.number_of_nodes()}</span></div>
    <div class="stat-row"><span class="stat-label">Edges:</span><span class="stat-value">{G.number_of_edges()}</span></div>
    <div class="stat-row"><span class="stat-label">Density:</span><span class="stat-value">{nx.density(G):.4f}</span></div>
    <div class="stat-row"><span class="stat-label">Avg Degree:</span><span class="stat-value">{sum(degree_dict.values()) / len(degree_dict):.1f}</span></div>
    <div class="stat-row"><span class="stat-label">Components:</span><span class="stat-value">{nx.number_connected_components(G)}</span></div>
    <div style="margin-top: 15px; padding-top: 10px; border-top: 2px solid #eee;">
        <div style="font-weight: bold; color: #764ba2; margin-bottom: 8px;">Top Hub Nodes:</div>
"""
for i, (node, degree) in enumerate(top_nodes[:5], 1):
    custom_header += f'        <div style="margin: 3px 0; font-size: 11px;">#{i} {node} ({degree})</div>\n'
custom_header += """
    </div>
</div>
"""


# Insert custom header after <body>
html_content = html_content.replace('<body>', '<body>\n' + custom_header)


with open(output_file, 'w', encoding='utf-8') as f:
    f.write(html_content)


print(f"\n{'='*60}")
print("✅ Interactive visualization created successfully!")
print(f"{'='*60}")
print(f"File: {output_file}")
print(f"Nodes: {G.number_of_nodes()}")
print(f"Edges: {G.number_of_edges()}")
print(f"Density: {nx.density(G):.4f}")
print(f"\nOpen '{output_file}' in your web browser to explore!")
print(f"{'='*60}")
