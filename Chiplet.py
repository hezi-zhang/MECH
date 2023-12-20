import networkx as nx
import matplotlib.pyplot as plt
from CouplingGraph import *

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

def gen_chiplet_array(geometry, array_x_dim, array_y_dim, chiplet_x_size, chiplet_y_size, cross_link_rows=None, cross_link_cols=None, auto=True):

    if geometry == 'square':
        G=gen_square_lattice(array_x_dim * chiplet_x_size, array_y_dim * chiplet_y_size)
        G.cross_link_rows = set(range(chiplet_y_size))
        G.cross_link_cols = set(range(chiplet_x_size))
    if geometry == 'hexagon':
        if auto:
            if chiplet_x_size % 2 == 1:
                print('WARNING: hexagon chiplets should have even columns')
                chiplet_x_size -= 1
        G=gen_hexagon_lattice(array_x_dim * chiplet_x_size, array_y_dim * chiplet_y_size)
        G.cross_link_rows = set(range(chiplet_y_size))
        G.cross_link_cols = set(range((chiplet_y_size+1)%2,chiplet_x_size,2))
    if geometry == 'heavy_square':
        if auto:
            if chiplet_x_size % 2 == 1:
                print('WARNING: heavy square chiplets should have even columns')
                chiplet_x_size -= 1
            if chiplet_y_size % 2 == 1:
                print('WARNING: heavy square chiplets should have even rows')
                chiplet_y_size -= 1
        G=gen_heavy_square_lattice(array_x_dim * chiplet_x_size, array_y_dim * chiplet_y_size)
        G.cross_link_rows = set(range(0,chiplet_y_size,2))
        G.cross_link_cols = set(range(0,chiplet_x_size,2))
    if geometry == 'heavy_hexagon':
        if auto:
            if chiplet_x_size % 4 != 0:
                print('WARNING: heavy hexagon chiplets should have 4n columns')
                chiplet_x_size = chiplet_x_size // 4 * 4
            if chiplet_y_size % 4 != 0:
                print('WARNING: heavy hexagon chiplets should have 4n rows')
                chiplet_y_size = chiplet_y_size // 4 * 4
        G=gen_heavy_hexagon_lattice(array_x_dim * chiplet_x_size, array_y_dim * chiplet_y_size)
        G.cross_link_rows = set(range(0,chiplet_y_size,2))
        G.cross_link_cols = set(range(2,chiplet_x_size,4))
    G.structure = geometry
    G.array_x_dim = array_x_dim
    G.array_y_dim = array_y_dim
    G.chiplet_x_size = chiplet_x_size
    G.chiplet_y_size = chiplet_y_size
    G.highway_qubits = None
    G.highway_distance_dict = None
    G.local_coupling_graph = None
    G.highway_coupling_graph = None
    for node in G.nodes:
        G.nodes[node]['type'] = 'data'
        
    gen_node_ids(G, chiplet_x_size, chiplet_y_size)
    gen_link_types(G)
    if cross_link_rows:
        G.cross_link_rows &= set(cross_link_rows)
    if cross_link_cols:
        G.cross_link_cols &= set(cross_link_cols)
    gen_sparse_chip_links(G,G.cross_link_rows, G.cross_link_cols)
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
    
    even_flag = 1
    for snode in shortest_path:
        if 'used_times' in G.nodes[snode].keys():
            G.nodes[snode]['used_times'] = G.nodes[snode]['used_times'] + 1
        else:
            G.nodes[snode]['used_times'] = 1
        cross_chip_flag = 0
        for neigh_node in G.neighbors(snode):
            if G[snode][neigh_node]['type'] == 'cross_chip':
                cross_chip_flag = 1
                break  
        
        if cross_chip_flag:
            G.nodes[snode]['type'] = 'highway'
        
        if even_flag:
            G.nodes[snode]['type'] = 'highway'
            even_flag = 0
        else:
            if G.nodes[snode]['type'] == 'data':
                even_flag = 1
    begin_node = shortest_path[0]
    end_node = shortest_path[-1]
    G.nodes[begin_node]['type'] = 'highway'
    G.nodes[end_node]['type'] = 'highway'
    return


def potential_highway_nodes_within_radius(G, node, radius):
    neighborhood = nx.ego_graph(G, node, radius=radius, center=False)
    filtered_nodes = [node for node in neighborhood.nodes if G.nodes[node].get('type') in {'highway', 'undetermined'}]
    return filtered_nodes

def deal_with_undetermined_nodes(G, adhoc_dense=True):
    for node in G.nodes:
        if G.nodes[node]['type'] == 'highway':
            index = 0
            for neigh_node in G.neighbors(node):
                if G.nodes[neigh_node]['type'] == 'highway':
                    index += 1
            if index >= G.nodes[node]['used_times'] * 2:
                cross_chip_flag = 0
                for neigh_node in G.neighbors(node):
                    if G[node][neigh_node]['type'] == 'cross_chip':
                        cross_chip_flag = 1
                        break
                if cross_chip_flag == 0:  
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
                # array_y = left_most_id[1]
                # print(array_x, array_y,'?',(array_x * G.chiplet_x_size, array_y * G.chiplet_y_size + row), ((array_x+1) * G.chiplet_x_size -1, array_y * G.chiplet_y_size + row))
                
                gen_interleaving_path_between(G, (array_x * G.chiplet_x_size, row), ((array_x+1) * G.chiplet_x_size -1, row))
        if col:
            bottom_most_id, up_most_id = G.nodes[(col, 0)]['id'], G.nodes[(col, max_y-1)]['id']
            for array_y in range(bottom_most_id[1], up_most_id[1]+1):
                # array_x = left_most_id[0]
                gen_interleaving_path_between(G, (col, array_y * G.chiplet_y_size), (col, (array_y+1) * G.chiplet_y_size -1), offset=offset)

def gen_highway_layout(G):
    highway_row = min(G.cross_link_rows, key=lambda x: abs(x - G.chiplet_y_size/2))
    highway_col = min(G.cross_link_cols, key=lambda x: abs(x - G.chiplet_x_size/2))
    if highway_col <  G.chiplet_x_size/2:
        offset = 'right'
    else:
        offset = 'left'
    for y in range(G.array_y_dim):
        row = y * G.chiplet_y_size + highway_row
        gen_road_along(G, row=row)
        
    for x in range(G.array_x_dim):
        col = x * G.chiplet_x_size + highway_col
        gen_road_along(G, col=col, offset=offset)
        
    deal_with_undetermined_nodes(G)


def gen_idx_qubit_dict(chiplet):
    idx_qubit_dict = dict()
    data_idx, highway_idx = 0, len(get_highway_qubits(chiplet))
    for node in sorted(chiplet.nodes):
        if get_node_type(chiplet, node)  == 'data':
            idx_qubit_dict[data_idx] = node
            data_idx += 1
        else:
            idx_qubit_dict[highway_idx] = node
            highway_idx += 1
    return idx_qubit_dict

def gen_qubit_idx_dict(chiplet):
    qubit_idx_dict = dict()
    data_idx, highway_idx = 0, len(get_highway_qubits(chiplet))
    for node in sorted(chiplet.nodes):
        if get_node_type(chiplet, node)  == 'data':
            qubit_idx_dict[node] = data_idx
            data_idx += 1
        else:
            qubit_idx_dict[node] = highway_idx
            highway_idx += 1
    return qubit_idx_dict


def draw_lattice(graph, size=8, node_size=300, with_labels=True, border=True, data_color='#99CCFF', highway_color='#004C99', on_chip_width=1, cross_link_width=8, fig_name=None, save=False):
        node_lable_dict = gen_qubit_idx_dict(graph)
        
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
        edge_width = [cross_link_width if graph.edges[edge].get('type', None) == 'cross_chip' else on_chip_width for edge in graph.edges]

        if border:
            edgecolors = ['black' for node in graph.nodes]
        else:
            edgecolors = data_color
        f = plt.figure(figsize=(size,size))
        ids = {node: graph.nodes[node].get('id', node)[-2:] for node in graph.nodes}
        nx.draw(graph, pos={(i,j): (i,j) for i, j in graph.nodes()}, with_labels=with_labels, labels=node_lable_dict, font_size=8, node_size=node_size, node_color=node_colors, edgecolors=edgecolors, linewidths=0.5, edge_color=edge_color, width=edge_width)
        if save:
            if fig_name is None:
                fig_name = '{}x{}_{}x{}_'.format(graph.array_x_dim, graph.array_y_dim, graph.chiplet_x_size, graph.chiplet_y_size) + graph.structure + '_' + data_color
            print('./figures/'+ fig_name + '.pdf')
            plt.savefig('./figures/'+ fig_name + '.pdf')


def get_local_coupling_graph(G):
    if G.local_coupling_graph is None:
        coupling_graph = G.copy()
        nodes_to_remove = []
        for node in G.nodes:
            if G.nodes[node]['type'] == 'highway':
                nodes_to_remove.append(node)
        coupling_graph.remove_nodes_from(nodes_to_remove)
        for node in coupling_graph.nodes:
            coupling_graph.nodes[node]['pos'] = node
        G.local_coupling_graph = CouplingGraph(coupling_graph)
    return G.local_coupling_graph


def get_highway_coupling_graph(G):
    if G.local_coupling_graph is None:
        coupling_graph = nx.Graph()
        highway_qubits = set()
        for node in G.nodes:
            if G.nodes[node]['type'] == 'highway':
                coupling_graph.add_node(node, pos=node)
                highway_qubits.add(node)
        
        for node in highway_qubits:
            neighbors_of_higway_neighbors = set()
            highway_neighbors_in_1_step = potential_highway_nodes_within_radius(G, node, radius=1)
            for nei in highway_neighbors_in_1_step:
                coupling_graph.add_edge(node, nei)
                neighbors_of_higway_neighbors = neighbors_of_higway_neighbors.union(set(G.neighbors(nei)))
            highway_neighbors_in_2_step = potential_highway_nodes_within_radius(G, node, radius=2)
            filtered_highway_neighbors = set(highway_neighbors_in_2_step) - set(highway_neighbors_in_1_step) - neighbors_of_higway_neighbors
            for nei in filtered_highway_neighbors:
                coupling_graph.add_edge(node, nei)
        G.local_coupling_graph = CouplingGraph(coupling_graph)
    return G.local_coupling_graph


def are_regularly_connected(G, node_1, node_2):
    if get_node_type(G, node_1) == 'highway' or get_node_type(G, node_2) == 'highway':
        return False
    else:
        return G.has_edge(node_1, node_2)

def get_node_type(G, node):
    return G.nodes[node].get('type')

def get_highway_qubits(G):
    if G.highway_qubits is None:
        G.highway_qubits = set([node for node in G.nodes if get_node_type(G, node) == 'highway'])
    return G.highway_qubits

def is_qubit_next_to_highway(G, node):
    if get_node_type(G, node) == 'highway':
        return False
    for nei in list(G.neighbors(node)):
        if G.nodes[nei].get('type') == 'highway':
            return True
    return False

def qubits_next_to_highway(G):
    result = set()
    highway_qubtis = get_highway_qubits(G)
    for highway_qubit in highway_qubtis:
        for nei in list(G.neighbors(highway_qubit)):
            if G.nodes[nei].get('type') == 'data':
                result.add(nei)
    return result


def get_distance_to_highway(G, node):
    if G.highway_distance_dict is None: 
        G.highway_distance_dict = {}
        for source_node in G.nodes():
            shortest_paths = nx.single_source_dijkstra_path_length(G, source_node)
            G.highway_distance_dict[source_node] = min(shortest_paths[highway_qubit] for highway_qubit in get_highway_qubits(G))

    return G.highway_distance_dict[node]


def find_possible_entrances_with_paths(G, node, range_beyond_closest=2):
        max_steps=int(get_distance_to_highway(G, node)) + range_beyond_closest
        entrance_with_paths = []
        nearby_entrance_cadidates = set(potential_highway_nodes_within_radius(G, node, radius=max_steps))
        for entrance in nearby_entrance_cadidates:
            all_shortest_paths = nx.all_shortest_paths(G, node, entrance)
            for path in all_shortest_paths:
                if len(set(path).intersection(get_highway_qubits(G))) == 1:
                    entrance_with_paths.append((entrance, tuple(path[:-1])))
        return entrance_with_paths
