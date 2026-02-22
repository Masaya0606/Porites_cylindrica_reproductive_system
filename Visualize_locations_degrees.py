import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

# **Read Files**
degree1 = pd.read_csv("king.kin_degree1_PoritesPcylInoha2503_rm57_Pcyl_24_40_Pcyl_Pc24_37_PcylYellow43_geno005_ed", sep="\t", engine="python")
degree2 = pd.read_csv("king.kin_degree2_PoritesPcylInoha2503_rm57_Pcyl_24_40_Pcyl_Pc24_37_PcylYellow43_geno005_ed", sep="\t", engine="python")
degree3 = pd.read_csv("king.kin_degree3_PoritesPcylInoha2503_rm57_Pcyl_24_40_Pcyl_Pc24_37_PcylYellow43_geno005_ed", sep="\t", engine="python")
all_kinship = pd.read_csv("king.kin_all_PoritesPcylInoha2503_rm57_Pcyl_24_40_Pcyl_Pc24_37_PcylYellow43_geno005_ed", sep="\t", engine="python")
location_info = pd.read_csv("ColonySize_locationInfo_ID_color_matchedWithGPS_ed_incl21", sep="\t", engine="python")
additional_location_info = pd.read_csv("ColonySize_locationInfo_ID_color_matchedWithGPS_ed_incl21", sep="\t", engine="python")

# **Define criteria for Clone (only)**
criteria = {
    "Clone": {"HetHet": 0.1, "IBS0": 0.001, "Kinship": 0.45}
}

# **Separate Clone from Degree1**
clones = degree1[
    (degree1["HetHet"] >= criteria["Clone"]["HetHet"]) &
    (degree1["IBS0"] <= criteria["Clone"]["IBS0"]) &
    (degree1["Kinship"] >= criteria["Clone"]["Kinship"])
].copy()

# **Assign Clone Relationship Label**
clones["Relationship"] = "Clone"
print("Total Degree1 Entries:", degree1.shape)
print("Total Clone Entries:", clones.shape)  # This should be a subset of Degree1

# Preview the extracted Clone relationships
print("Clones DataFrame:")
print(clones.head())

# Remove clone pairs from Degree1
filtered_degree1 = degree1.drop(clones.index).copy()
filtered_degree1["Relationship"] = "Degree1"

# **Degree1 pairs**
degree1_pairs = set(zip(filtered_degree1["ID1"], filtered_degree1["ID2"]))

# **Remove Degree1 pairs from Degree2**
filtered_degree2 = degree2[
    ~degree2.apply(lambda row: (row["ID1"], row["ID2"]) in degree1_pairs or (row["ID2"], row["ID1"]) in degree1_pairs, axis=1)
].copy()
filtered_degree2["Relationship"] = "Degree2"

# **Degree2 pairs**
degree2_pairs = set(zip(filtered_degree2["ID1"], filtered_degree2["ID2"]))

# **Remove Degree1 and Degree2 pairs from Degree3**
filtered_degree3 = degree3[
    ~degree3.apply(lambda row: (row["ID1"], row["ID2"]) in degree1_pairs or (row["ID2"], row["ID1"]) in degree1_pairs, axis=1)
].copy()
filtered_degree3 = filtered_degree3[
    ~filtered_degree3.apply(lambda row: (row["ID1"], row["ID2"]) in degree2_pairs or (row["ID2"], row["ID1"]) in degree2_pairs, axis=1)
].copy()
filtered_degree3["Relationship"] = "Degree3"

# **Concatenate all relationships**
all_relationships = pd.concat([clones, filtered_degree1, filtered_degree2, filtered_degree3])

#Debugging Step: Verify relationships are included correctly
print("All relationships counts:\n", all_relationships["Relationship"].value_counts())

# **Check Data**
print("Clone:", clones.shape)
print("Filtered Degree1:", filtered_degree1.shape)
print("Filtered Degree2:", filtered_degree2.shape)
print("Filtered Degree3:", filtered_degree3.shape)
print("All relationships counts:\n", all_relationships["Relationship"].value_counts())

# **Color mapping**
color_map = {"Clone": "red", "Degree1": "blue", "Degree2": "green", "Degree3": "orange"}
sex_color_map = {"M": "lightblue", "F": "pink", "MH": "#6FA8DC", "FH": "magenta"}
edge_color_map = {"Y": "yellow", "B": "brown"}

# **Create network graph**
G = nx.Graph()

# **Ensure Clone edges are added first**
for _, row in clones.iterrows():
    G.add_edge(row["ID1"], row["ID2"], relationship="Clone", color=color_map["Clone"])

# **Then add other relationships**
for _, row in all_relationships.iterrows():
  if not G.has_edge(row["ID1"], row["ID2"]):  # Prevent overwriting Clone edges
       G.add_edge(row["ID1"], row["ID2"], relationship=row["Relationship"], color=color_map[row["Relationship"]])

# Debugging step: Print all Clone edges
clone_edges = [(u, v) for u, v, d in G.edges(data=True) if d.get("relationship") == "Clone"]
print("Clone edges in graph:", clone_edges)


# **Merge all individual location data**
combined_location_info = pd.concat([location_info, additional_location_info]).drop_duplicates(subset="ID").set_index("ID")

# **Add all individuals to G**
all_individuals = set(all_kinship["ID1"]).union(set(all_kinship["ID2"]), set(location_info["ID"]), set(additional_location_info["ID"]))
G.add_nodes_from(all_individuals)

# **Obtain position information (Set default position if missing)**
positions = {
    node: (combined_location_info.loc[node, "Longitude"], combined_location_info.loc[node, "Latitude"])
    if node in combined_location_info.index and "Longitude" in combined_location_info.columns and "Latitude" in combined_location_info.columns
    else (None, None)
    for node in all_individuals
}

# **Ensure no None values in positions**
positions = {k: v if v != (None, None) else (0, 0) for k, v in positions.items()}

# **Set node colors based on sex**
node_colors = [
    sex_color_map.get(combined_location_info.loc[node, "sex"], "gray") if node in combined_location_info.index else "gray"
    for node in G.nodes
]

# **Adjust node sizes**
node_sizes = [
    (combined_location_info.loc[node, "MeanDiameter"] / combined_location_info["MeanDiameter"].max()) * 1000
    if node in combined_location_info.index and pd.notna(combined_location_info.loc[node, "MeanDiameter"])
    else 100
    for node in G.nodes
]

# **Set node border colors based on `color`**
node_border_colors = [
    edge_color_map.get(combined_location_info.loc[node, "color"], "gray") if node in combined_location_info.index and pd.notna(combined_location_info.loc[node, "color"])
    else "gray"
    for node in G.nodes
]

# **Set edge colors**
edge_colors = [G[u][v]["color"] if "color" in G[u][v] else "black" for u, v in G.edges]

# **Plot**
plt.figure(figsize=(16, 12))
ax = plt.gca()

# **Draw nodes with borders**
for node, (x, y) in positions.items():
    idx = list(G.nodes).index(node)  # Get node index
    ax.scatter(
        x, y, 
        s=node_sizes[idx], 
        facecolors=node_colors[idx],  # **Node interior color**
        edgecolors=node_border_colors[idx],  # **Border color**
        linewidths=2,
        zorder=3
    )

# **Draw normal edges**
nx.draw_networkx_edges(G, pos=positions, edge_color=edge_colors, width=2, alpha=0.6)

# **Add labels**
nx.draw_networkx_labels(G, pos=positions, font_size=10)

# **Create legends**
degree_legend = [Line2D([0], [0], color=color, lw=4, label=degree) for degree, color in color_map.items()]
sex_legend = [Line2D([0],[0],marker='o',color='w',markersize=10, markerfacecolor=color,label=sex) for sex, color in sex_color_map.items()]
border_legend =[Line2D([0],[0],marker='o',color='w',markersize=10, markeredgewidth=2, markeredgecolor=color, label=color_name) for color_name, color in edge_color_map.items()]

legend_elements = degree_legend+sex_legend+border_legend

# **Place legend in upper right**
plt.legend(handles=legend_elements, loc='upper left', bbox_to_anchor=(1.05, 1), fontsize=12)

# **Add Scale Bar (10m)**
from matplotlib_scalebar.scalebar import ScaleBar
import numpy as np

# 研究対象の緯度 (例: データの緯度範囲で中央付近の値を指定)
latitude = 26.651738  # データの中心付近の緯度

# この緯度における1度の緯度と経度の距離 (メートル) を計算
meters_per_degree_lat = 111320  # 1度の緯度 = 約111.32 km
meters_per_degree_lon = 111320 * np.cos(np.radians(latitude))  # 緯度に依存する1度の経度の距離

# 今回は経度方向に対するスケールを使用する
scalebar = ScaleBar(meters_per_degree_lon, units='m', location='lower left', length_fraction=0.2)
ax.add_artist(scalebar)


plt.title("Geographic Network with Kinship, sex, and color of Porites cylindrica")
plt.xlabel("Longitude")
plt.ylabel("Latitude")
plt.tight_layout()
plt.savefig("sex_clone_gps.pdf", format='pdf', bbox_inches='tight')
plt.show()
