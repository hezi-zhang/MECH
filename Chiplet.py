import networkx as nx
import matplotlib.pyplot as plt

def gen_square_lattice(x_dim, y_dim):
    return nx.grid_2d_graph(x_dim, y_dim)

def gen_hexagon_lattice(x_dim, y_dim):
    G = nx.grid_2d_graph(x_dim, y_dim)
    for node in G.nodes:
        if node[1] < y_dim-1 and (node[0] + node[1])%2 == 1:
            G.remove_edge(node, (node[0], node[1]+1))
    return G

def gen_heavy_square_lattice(x_dim, y_dim):
    G = nx.grid_2d_graph(x_dim, y_dim)
    nodes_to_remove = []
    for node in G.nodes:
        if node[0]%2==1 and node[1]%2==1:
           nodes_to_remove.append(node)
    G.remove_nodes_from(nodes_to_remove)
    return G 

def gen_heavy_hexagon_lattice(x_dim, y_dim):
    G=gen_heavy_square_lattice(x_dim, y_dim)
    nodes_to_remove = []
    for x,y in G.nodes:
        if y<y_dim-1 and x%2==0 and y%2==0 and (x+y)%4==2:
            nodes_to_remove.append((x,y+1))
    G.remove_nodes_from(nodes_to_remove)
    return G 

def gen_node_ids(lattice, chiplet_x_size, chiplet_y_size):
    for x, y in lattice.nodes:
        lattice.nodes[(x,y)]['id'] = (x//chiplet_x_size, y//chiplet_y_size, x%chiplet_x_size, y%chiplet_y_size)
        
def gen_link_types(lattice):
    for u,v in lattice.edges:
        if lattice.nodes[u]['id'][:2] == lattice.nodes[v]['id'][:2]:
            lattice.edges[(u,v)]['type'] = 'on_chip'
        else:
            lattice.edges[(u,v)]['type'] = 'cross_chip'

def gen_sparse_chip_links(lattice, rows, cols):
    edges_to_remove = []
    for u,v in lattice.edges:
        u_id, v_id = lattice.nodes[u]['id'], lattice.nodes[v]['id']
        if lattice.edges[(u,v)]['type'] == 'cross_chip':
            if u_id[2] == v_id[2] and u_id[2] not in cols:
                edges_to_remove.append((u,v))
            elif u_id[3] == v_id[3] and u_id[3] not in rows:
                edges_to_remove.append((u,v))
    lattice.remove_edges_from(edges_to_remove)

def gen_chiplet_array(geometry,array_x_dim, array_y_dim, chiplet_x_size, chiplet_y_size, cross_link_rows=None, cross_link_cols=None):
    if geometry == 'square':
        G=gen_square_lattice(array_x_dim * chiplet_x_size, array_y_dim * chiplet_y_size)
    if geometry == 'hexagon':
        G=gen_hexagon_lattice(array_x_dim * chiplet_x_size, array_y_dim * chiplet_y_size)
    if geometry == 'heavy_square':
        G=gen_heavy_square_lattice(array_x_dim * chiplet_x_size, array_y_dim * chiplet_y_size)
    if geometry == 'heavy_hexagon':
        G=gen_heavy_hexagon_lattice(array_x_dim * chiplet_x_size, array_y_dim * chiplet_y_size)
    gen_node_ids(G, chiplet_x_size, chiplet_y_size)
    gen_link_types(G)
    if cross_link_rows is None:
        cross_link_rows = list(range(chiplet_y_size))
    if cross_link_cols is None:
        cross_link_cols = list(range(chiplet_y_size))
    gen_sparse_chip_links(G,cross_link_rows, cross_link_cols)
    return G


def draw_lattice(graph, size=8, with_labels=True, border=True, data_color='#99CCFF', highway_color='#004C99'):
        edge_color = ['red' if graph.edges[edge].get('type', None) == 'cross_chip' else 'black' for edge in graph.edges]
        edge_width = [8 if graph.edges[edge].get('type', None) == 'cross_chip' else 1 for edge in graph.edges]

        if border:
            edgecolors = ['black' for node in graph.nodes]
        else:
            edgecolors = data_color
        f = plt.figure(figsize=(size,size))
        labels = {node: graph.nodes[node].get('id', node)[-2:] for node in graph.nodes}
        nx.draw(graph, pos={(i,j): (i,j) for i, j in graph.nodes()}, with_labels=with_labels, labels=labels, font_size=10, node_color=data_color, edgecolors=edgecolors, linewidths=0.5, edge_color=edge_color, width=edge_width)
        