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
        assert not source in self.highway_qubits, "source qubit {} appears on highway".format(source)
        for target in targets:
            assert not target in self.highway_qubits, "target qubit {} appear on highway".format(target)
        for target in targets:
            path += self.find_free_highway_path(source, target, allow_sharing_source)
        return path
    
    def draw(self, size=8, with_labels=True, border=True, data_color='#99CCFF', highway_color='#004C99'):
        node_lable_dict = gen_qubit_idx_dict(self.chiplet_array)

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
        nx.draw(graph, pos={(i,j): (i,j) for i, j in graph.nodes()}, with_labels=with_labels, labels=node_lable_dict, font_size=8, node_color=node_colors, edgecolors=edgecolors, linewidths=0.5, edge_color=edge_color, width=edge_width)


class HighwayManager:
    def __init__(self, G, prep_period, meas_period):
        self.chiplet_array = G
        self.prep_period = prep_period
        self.meas_period = meas_period
        self.highway_schedule = defaultdict(dict) # shuttle idx --> {prep_start_time, exec_start_time, meas_start_time, end_time}
        self.shuttle_stack = []                   # a list of HighwayOccupancy
        self.qubit_idx_dict = gen_qubit_idx_dict(self.chiplet_array)
        self.idx_qubit_dict = gen_idx_qubit_dict(self.chiplet_array)
    def __copy__(self):
        new_instance = copy.deepcopy(self)
        return new_instance
    
    def get_shuttle_prep_start_time(self, shuttle_idx):
        return self.highway_schedule[shuttle_idx]['prep_start_time']
    def get_shuttle_exec_start_time(self, shuttle_idx):
        return self.highway_schedule[shuttle_idx]['exec_start_time']
    def get_shuttle_meas_start_time(self, shuttle_idx):
        return self.highway_schedule[shuttle_idx]['meas_start_time']
    def get_shuttle_end_time(self, shuttle_idx):
        return self.highway_schedule[shuttle_idx]['end_time']

    def set_shuttle_prep_start_time(self, shuttle_idx, time):
        self.highway_schedule[shuttle_idx]['prep_start_time'] = time
    def set_shuttle_exec_start_time(self, shuttle_idx, time):
        self.highway_schedule[shuttle_idx]['exec_start_time'] = time
    def set_shuttle_meas_start_time(self, shuttle_idx, time):
        self.highway_schedule[shuttle_idx]['meas_start_time'] = time
    def set_shuttle_end_time(self, shuttle_idx, time):
        self.highway_schedule[shuttle_idx]['end_time'] = time

    def end_shuttle(self, shuttle_idx, circuit):
        mop_depth = circuit.mop_depth
        self.set_shuttle_meas_start_time(shuttle_idx, mop_depth)
        self.set_shuttle_end_time(shuttle_idx, mop_depth + self.meas_period - 1)

    def allocate_shuttle(self, circuit):
        if not self.shuttle_stack:
            start_idx = 0
            new_shuttle_idx = 0
        else:
            last_shuttle_idx = len(self.shuttle_stack)-1
            self.end_shuttle(last_shuttle_idx, circuit)
            start_idx = self.get_shuttle_end_time(last_shuttle_idx) + 1
            new_shuttle_idx = last_shuttle_idx + 1

        self.shuttle_stack.append(HighwayOccupancy(self.chiplet_array))
        self.set_shuttle_prep_start_time(new_shuttle_idx, start_idx)
        self.set_shuttle_exec_start_time(new_shuttle_idx, start_idx + self.prep_period)
        self.set_shuttle_meas_start_time(new_shuttle_idx, None)
        self.set_shuttle_end_time(new_shuttle_idx, None)
        

    def which_shuttle(self, target_time):
        last_shuttle_idx = len(self.shuttle_stack) - 1
        if last_shuttle_idx < 0:
            return -1
        if self.get_shuttle_meas_start_time(last_shuttle_idx) is None and target_time >= self.get_shuttle_prep_start_time(last_shuttle_idx):
            return last_shuttle_idx

        low, high = 0, len(self.shuttle_stack) - 1
        while low <= high:
            mid = (low + high) // 2

            if self.get_shuttle_prep_start_time(mid) <= target_time <= self.get_shuttle_end_time(mid):
                if target_time < self.get_shuttle_meas_start_time(mid):
                    return mid
                elif mid + 1 < len(self.shuttle_stack):
                    return mid + 1
                else:
                    return -1
            elif target_time < self.get_shuttle_prep_start_time(mid):
                high = mid - 1
            else:
                low = mid + 1
        return -1
    
    def which_shuttle_exec_period(self, target_time):
        last_shuttle_idx = len(self.shuttle_stack) - 1
        if last_shuttle_idx < 0:
            return -1
        if self.get_shuttle_meas_start_time(last_shuttle_idx) is None and target_time >= self.get_shuttle_exec_start_time(last_shuttle_idx):
            return last_shuttle_idx

        low, high = 0, len(self.shuttle_stack) - 1
        while low <= high:
            mid = (low + high) // 2

            if self.get_shuttle_exec_start_time(mid) <= target_time < self.get_shuttle_meas_start_time(mid):
                return mid  # Found a matching interval
            elif target_time < self.get_shuttle_exec_start_time(mid):
                high = mid - 1
            else:
                low = mid + 1
        return -1
    
    def is_in_exec_period(self, target_time):
        return self.which_shuttle_exec_period(target_time) >= 0
    
    
    def is_entrance_idle_at_idx(self, circuit, shuttle_idx, target_entrance, idx):
        highway_shuttle = self.shuttle_stack[shuttle_idx]
        for data_qubit in highway_shuttle.entrance_data_dict[target_entrance]:
            
            line = self.qubit_idx_dict[data_qubit]
            if circuit.take_role(line, idx) in {'mc', 'mt'}:
                return False
        return True

    def get_earliest_entrance_idle_time_on_shuttle(self, circuit, shuttle_idx, target_entrance):
        exec_start_time = self.get_shuttle_exec_start_time(shuttle_idx)
        meas_time = self.get_shuttle_meas_start_time(shuttle_idx)
        if meas_time is not None:
            latest_idx = meas_time - 1
        else:
            latest_idx = circuit.depth - 1

        earliest_idx = None
        for idx in range(latest_idx, exec_start_time-1, -1):
            if self.is_entrance_idle_at_idx(circuit, shuttle_idx, target_entrance, idx):
                earliest_idx = idx

        if earliest_idx is not None:
            return earliest_idx
        else:
            if meas_time is None:
                return latest_idx + 1
            else:
                return -1
            
    def get_earliest_index_for_2qubit_component_on_shuttle(self, shuttle_idx, circuit, op, auto_commuting=True, auto_cancellation=False):
        control, target = op.control, op.target
        
        control_qubit, target_qubit = self.idx_qubit_dict[op.control], self.idx_qubit_dict[op.target]
        highway_shuttle = self.shuttle_stack[shuttle_idx]
        if highway_shuttle.source_target_counter[(control_qubit, target_qubit)] > 0:
            target_entrance = highway_shuttle.target_entrance_dict[op.target]
        else:
            path = highway_shuttle.find_free_highway_path(control_qubit, target_qubit)
            if not path:
                return -1
            target_entrance = path[-2]

        earliest_entrance_idle_idx = self.get_earliest_entrance_idle_time_on_shuttle(circuit, shuttle_idx, target_entrance)
        if earliest_entrance_idle_idx < 0:
            return -1

        meas_start_time = self.get_shuttle_meas_start_time(shuttle_idx)
        if meas_start_time is not None:
            latest_idx = meas_start_time - 1
        else:
            latest_idx = max(circuit.get_line_depth(control), circuit.get_line_depth(target))
            exec_start_time = self.get_shuttle_exec_start_time(shuttle_idx)
            if latest_idx <= exec_start_time:
                return exec_start_time
        
        earliest_idx = -1
        for idx in range(latest_idx, -1, -1):
            c_node, c_role = circuit.take_node(control, idx), circuit.take_role(control, idx)
            t_node, t_role = circuit.take_node(target, idx), circuit.take_role(target, idx)
            if circuit.is_node_the_same_gate(c_node, control, target) or circuit.is_component_contained_in_mop(c_node, control, target):
                if auto_cancellation:
                    return idx
                elif c_role == 'mc':
                    return earliest_idx
            
            is_target_idle = self.is_entrance_idle_at_idx(circuit, shuttle_idx, target_entrance, idx)
            if c_node is None and t_node is None and is_target_idle:
                earliest_idx = idx
            elif c_role == 'mc' and t_node is None and is_target_idle:
                if not auto_commuting:
                    return idx
                else:
                    earliest_idx = idx
            elif not auto_commuting or c_role not in {None, 'c', 'mc'}:
                return -1           # not commutable at control line
            elif t_role not in {None, 't', 'mt'}:
                return earliest_idx # not commutable at control line
                     
        return earliest_idx
            
    def get_highway_aware_earliest_index_for_2qubit_component(self, circuit, op, auto_commuting=True, auto_cancellation=False):
        highway_agnostic_earliest_idx = circuit.get_earliest_index_for_2qubit_component(op, auto_commuting=auto_commuting, auto_cancellation=auto_cancellation)
        earliest_shuttle_idx = self.which_shuttle(highway_agnostic_earliest_idx)
        print('earliest_shuttle = ', earliest_shuttle_idx)
        for shuttle_idx in range(earliest_shuttle_idx, len(self.shuttle_stack) - 1):
            this_shuttle_possible = True
            for idx in range(highway_agnostic_earliest_idx, self.get_shuttle_exec_start_time(shuttle_idx), -1):
                prior_role = circuit.take_role(op.control, idx - 1)
                if prior_role not in {None, 'c', 'mc'}:
                    this_shuttle_possible = False

            if this_shuttle_possible:
                exec_idx = self.get_earliest_index_for_2qubit_component_on_shuttle(shuttle_idx, circuit, op, auto_commuting=auto_commuting, auto_cancellation=auto_cancellation)
                
                if exec_idx >= 0:
                    return exec_idx
        return -1
    