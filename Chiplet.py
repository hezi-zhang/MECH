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

def gen_chiplet_array(geometry, array_x_dim, array_y_dim, chiplet_x_size, chiplet_y_size, cross_link_rows=None, cross_link_cols=None):

    if geometry == 'square':
        G=gen_square_lattice(array_x_dim * chiplet_x_size, array_y_dim * chiplet_y_size)
    if geometry == 'hexagon':
        G=gen_hexagon_lattice(array_x_dim * chiplet_x_size, array_y_dim * chiplet_y_size)
    if geometry == 'heavy_square':
        G=gen_heavy_square_lattice(array_x_dim * chiplet_x_size, array_y_dim * chiplet_y_size)
    if geometry == 'heavy_hexagon':
        G=gen_heavy_hexagon_lattice(array_x_dim * chiplet_x_size, array_y_dim * chiplet_y_size)
    G.array_x_dim = array_x_dim
    G.array_y_dim = array_y_dim
    G.chiplet_x_size = chiplet_x_size
    G.chiplet_y_size = chiplet_y_size
    for node in G.nodes:
        G.nodes[node]['type'] = 'data'
        
    gen_node_ids(G, chiplet_x_size, chiplet_y_size)
    gen_link_types(G)
    if cross_link_rows is None:
        cross_link_rows = list(range(chiplet_y_size))
    if cross_link_cols is None:
        cross_link_cols = list(range(chiplet_y_size))
    gen_sparse_chip_links(G,cross_link_rows, cross_link_cols)
    return G

def gen_interleaving_path_between(G, source, target, adhoc_dense=True, offset=None):
    if offset == 'left' or offset == 'down':
        coor_offset = -1
    elif offset == 'right' or offset == 'up':
        coor_offset = 1
    else:
        coor_offset = 0
    all_shortest_paths = list(nx.all_shortest_paths(G, source, target))
    shortest_path = min(all_shortest_paths, key=lambda path: sum(abs(x-(source[0] + 0.5*coor_offset)) for x,_ in path))
    left_pointer = 0
    right_pointer = len(shortest_path) - 1

    while right_pointer - left_pointer > 2:
        left_node, right_node = shortest_path[left_pointer], shortest_path[right_pointer]
        G.nodes[left_node]['type'] = 'highway'
        G.nodes[right_node]['type'] = 'highway'

        left_pointer += 2
        right_pointer -= 2
    if left_pointer > right_pointer:
        left_pointer, right_pointer = right_pointer, left_pointer
    for node in shortest_path[left_pointer:right_pointer + 1]:
        G.nodes[node]['type'] = 'undetermined'
    return shortest_path[left_pointer:right_pointer + 1]

def potential_highway_nodes_within_radius(G, node, radius):
    neighborhood = nx.ego_graph(G, node, radius=radius, center=False)
    filtered_nodes = [node for node in neighborhood.nodes if G.nodes[node].get('type') in {'highway', 'undetermined'}]
    return filtered_nodes

def deal_with_undetermined_nodes(G, adhoc_dense=True):
    for node in G.nodes:
        if G.nodes[node]['type'] == 'undetermined':
            neighborhood = potential_highway_nodes_within_radius(G, node, 2)
            for nei in neighborhood:
                if not adhoc_dense:
                    x, y = nei
                    max_x, max_y = G.array_x_dim * G.chiplet_x_size, G.array_y_dim * G.chiplet_y_size
                else:
                    x, y = G.nodes[nei]['id'][-2:]
                    max_x, max_y = G.chiplet_x_size, G.chiplet_y_size
                is_on_edge = x==0 or x==max_x-1 or y==0 or y==max_y-1

                if G.nodes[nei]['type'] == 'highway':
                    least_neighbors = 2 - is_on_edge + 1
                if G.nodes[nei]['type'] == 'undetermined':
                    least_neighbors = 4 - is_on_edge
                # print(nei, is_on_edge, len(potential_highway_nodes_within_radius(G, nei, 2)), least_neighbors)
                if len(potential_highway_nodes_within_radius(G, nei, 2)) < least_neighbors:
                        G.nodes[node]['type'] = 'highway'
                        break
                else:
                    neighborhood = potential_highway_nodes_within_radius(G, nei, 2)
                    neighborhood.remove(node)
                    this_row = [n for n in neighborhood if n[1] == nei[1]]
                    if len(this_row) == 1:
                        G.nodes[node]['type'] = 'highway'
                        break
            if G.nodes[node]['type'] == 'undetermined': 
                G.nodes[node]['type'] = 'data' 

def gen_road_along(G, row=None, col=None, offset=None, adhoc_dense=True):
    max_x, max_y = G.array_x_dim * G.chiplet_x_size, G.array_y_dim * G.chiplet_y_size
    if not adhoc_dense:
        if row:
            gen_interleaving_path_between(G, (0, row), (max_x-1, row))
        if col:
            gen_interleaving_path_between(G, (col, 0), (col, max_y-1), offset=offset)
    else:
        if row:
            left_most_id, right_most_id = G.nodes[(0, row)]['id'], G.nodes[(max_x-1, row)]['id']
            for array_x in range(left_most_id[0], right_most_id[0]+1):
                array_y = left_most_id[1]
                gen_interleaving_path_between(G, (array_x * G.chiplet_x_size, array_y * G.chiplet_y_size + row), ((array_x+1) * G.chiplet_x_size -1, array_y * G.chiplet_y_size + row))
        if col:
            bottom_most_id, up_most_id = G.nodes[(col, 0)]['id'], G.nodes[(col, max_y-1)]['id']
            for array_y in range(bottom_most_id[1], up_most_id[1]+1):
                array_x = left_most_id[0]
                gen_interleaving_path_between(G, (array_x * G.chiplet_x_size + col, array_y * G.chiplet_y_size), (array_x * G.chiplet_x_size + col, (array_y+1) * G.chiplet_y_size -1), offset=offset)
     


def draw_lattice(graph, size=8, with_labels=True, border=True, data_color='#99CCFF', highway_color='#004C99'):
        node_colors = []
        for node in graph.nodes:
            node_type = graph.nodes[node].get('type', '')
            if node_type == 'highway':
                node_colors.append(highway_color)
            elif node_type == 'undetermined':
                node_colors.append('yellowgreen')
            else:
                node_colors.append(data_color)

        edge_color = ['red' if graph.edges[edge].get('type', None) == 'cross_chip' else 'black' for edge in graph.edges]
        edge_width = [8 if graph.edges[edge].get('type', None) == 'cross_chip' else 1 for edge in graph.edges]

        if border:
            edgecolors = ['black' for node in graph.nodes]
        else:
            edgecolors = data_color
        f = plt.figure(figsize=(size,size))
        ids = {node: graph.nodes[node].get('id', node)[-2:] for node in graph.nodes}
        nx.draw(graph, pos={(i,j): (i,j) for i, j in graph.nodes()}, with_labels=with_labels, font_size=8, node_color=node_colors, edgecolors=edgecolors, linewidths=0.5, edge_color=edge_color, width=edge_width)
        