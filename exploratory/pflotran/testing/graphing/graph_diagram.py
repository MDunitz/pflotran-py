import re
import networkx as nx
import matplotlib.pyplot as plt

def clean_species_name(name):
    if 'Tracer' in name:
        return None
    if name == "SO4--":
        return "SO₄²⁻"
    if name == "Fe++":
        return "Fe(II)"
    if name == "Fe+++":
        return "Fe(III)"
    return name

# Load reaction text
with open("reactions.txt", "r") as f:
    text = f.read()

# Extract reactions
reaction_blocks = re.findall(r"(MICROBIAL_REACTION|GENERAL_REACTION).*?REACTION\s+([^\n]+)", text, re.DOTALL)

G = nx.DiGraph()

species_pattern = re.compile(r"(\d*\.\d+e[+-]?\d+|\d+\.?\d*)\s+([A-Za-z0-9\-\+\(\)]+)")

for _, reaction_line in reaction_blocks:
    if '<->' in reaction_line:
        lhs, rhs = reaction_line.split('<->')
        reversible = True
    elif '->' in reaction_line:
        lhs, rhs = reaction_line.split('->')
        reversible = False
    else:
        continue

    reactants = species_pattern.findall(lhs)
    products = species_pattern.findall(rhs)

    reactant_names = [clean_species_name(name) for _, name in reactants]
    product_names = [clean_species_name(name) for _, name in products]

    reactant_names = [r for r in reactant_names if r]
    product_names = [p for p in product_names if p]

    for r in reactant_names:
        for p in product_names:
            G.add_edge(r, p)
            if reversible:
                G.add_edge(p, r)

# Purple for CH4(aq)
node_colors = ["purple" if node == "CH4(aq)" else "lightblue" for node in G.nodes()]

# Spring layout
pos = nx.spring_layout(G, k=2.5, iterations=500, seed=42)

# Manual node repulsion
def repel_close_nodes(pos, min_dist=0.15, iterations=10):
    nodes = list(pos.keys())
    for _ in range(iterations):
        for i, n1 in enumerate(nodes):
            for n2 in nodes[i+1:]:
                dx = pos[n2][0] - pos[n1][0]
                dy = pos[n2][1] - pos[n1][1]
                dist = (dx**2 + dy**2)**0.5
                if dist < min_dist and dist > 0:
                    shift = (min_dist - dist) / 2
                    dx_norm = dx / dist
                    dy_norm = dy / dist
                    pos[n1] = (pos[n1][0] - dx_norm*shift, pos[n1][1] - dy_norm*shift)
                    pos[n2] = (pos[n2][0] + dx_norm*shift, pos[n2][1] + dy_norm*shift)
    return pos

pos = repel_close_nodes(pos)

plt.figure(figsize=(20, 16))
nx.draw(
    G, pos,
    with_labels=True,
    node_color=node_colors,
    node_size=3000,
    font_size=12,
    font_weight='bold',
    edge_color='gray',
    arrows=True,
    arrowsize=30 
)
plt.title("Species Flow Diagram (cleaned)", fontsize=16)
plt.tight_layout()
plt.show()
