from collections import deque, defaultdict, Counter
from functools import *
from numpy import *
import networkx as nx
from Chiplet import *

class HighwayOccupancy:
    def __init__(self, G):
        self.chiplet_array = G
        self.highway_qubits = get_highway_qubits(G)
        self.occupied_highway_qubits = set()
        self.occupied_highway_segments = dict() # key = sorted_data_qubits, val = path_on_highway
        self.source_occupied_highway = defaultdict(set) # key = source_data_qubit, val = set_of_paths
        self.source_targets_dict = defaultdict(set) # key = source, val = set_of_targets
        self.source_target_counter = Counter() # key = (source, target)
        self.highway_entrance_workload = Counter() # key = entrance_on_highway
        self.target_entrance_dict = dict() # key = target, val = entrance_on_highway
        self.entrance_data_dict = defaultdict(set)
    def __copy__(self):
        new_instance = copy.deepcopy(self)
        return new_instance
    
    @property
    def free_highway_qubits(self):
        return self.highway_qubits - self.occupied_highway_qubits
    
    def occupy_highway_path(self, source, path):
        path_nodes = set(path)
        data_qubits = path_nodes - self.highway_qubits
        sorted_data_qubits = tuple(sorted(list(data_qubits)))

        self.source_occupied_highway[source] |= path_nodes - data_qubits
        self.occupied_highway_segments[sorted_data_qubits] = path_nodes - data_qubits
        self.occupied_highway_qubits |= path_nodes - data_qubits

        target = sorted_data_qubits[1] if sorted_data_qubits[0] == source else sorted_data_qubits[0]
    
        self.source_targets_dict[source].add(target)
        self.source_target_counter[(source, target)] += 1
        
        highway_entrance = path[-2]
        self.target_entrance_dict[target] = highway_entrance
        self.entrance_data_dict[highway_entrance].add(target)
        self.highway_entrance_workload[highway_entrance] += 1


    def find_free_highway_path(self, source, target, allow_sharing_source=True):
        if target in self.source_targets_dict.keys():
            return []
        existing_targets = set()
        for t in self.source_targets_dict.values():
            existing_targets = existing_targets.union(set(t))
        if source in existing_targets:
            return []
        
        subnodes = self.free_highway_qubits.union(set([source, target]))
        occupied_highway_by_same_source = set()
        if allow_sharing_source:
            occupied_highway_by_same_source = self.source_occupied_highway[source]
            subnodes = subnodes.union(occupied_highway_by_same_source)
        highway_coupling_graph = get_highway_coupling_graph(self.chiplet_array)
        subgraph = nx.subgraph(highway_coupling_graph.graph, subnodes).copy()
        subgraph.add_nodes_from([source, target])
        for nei in self.chiplet_array.neighbors(source):
            if nei in self.highway_qubits:
                subgraph.add_edge(source, nei)
        for nei in self.chiplet_array.neighbors(target):
            if nei in self.highway_qubits:
                subgraph.add_edge(target, nei)
        if subgraph.has_edge(source, target):
            subgraph.remove_edge(source, target)
        for n1,n2 in subgraph.edges:
            if n1 in occupied_highway_by_same_source and n2 in occupied_highway_by_same_source:
                subgraph.edges[(n1,n2)]['weight'] = 0
            else:
                subgraph.edges[(n1,n2)]['weight'] = 1
        if not nx.has_path(subgraph, source, target):
            return []
        else:
            path = nx.dijkstra_path(subgraph, source, target)
            if not is_qubit_next_to_highway(self.chiplet_array, source) or not is_qubit_next_to_highway(self.chiplet_array, target):
                return []
            else:
                self.occupy_highway_path(source, path)
                return path
        
    def find_free_highway_path_for_control_multi_targets(self, source, targets, allow_sharing_source=True):
        path = []
        assert not source in self.highway_qubits, "source qubit appears on highway"
        for target in targets:
            assert not target in self.highway_qubits, "target qubits appear on highway"
        for target in targets:
            path += self.find_free_highway_path(source, target, allow_sharing_source)
        return path
    
    def draw(self, size=8, with_labels=True, border=True, data_color='#99CCFF', highway_color='#004C99'):
        node_colors = []
        graph = self.chiplet_array
        end_nodes = []
        for key in self.occupied_highway_segments.keys():
            end_nodes += list(key)
        end_nodes = set(end_nodes)
        for node in graph.nodes:
            if node in self.occupied_highway_qubits:
                node_colors.append('orange')
            elif node in end_nodes:
                node_colors.append('yellow')
            elif get_node_type(self.chiplet_array, node) == 'highway':
                node_colors.append(highway_color)
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


