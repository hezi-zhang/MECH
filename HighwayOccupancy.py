from collections import deque, defaultdict, Counter
from functools import *
from numpy import *
import networkx as nx
from Chiplet import *
from Circuit import *

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
        self.highway_status_record = []           # record the highway status (which shuttle, which period) at each idx
        self.shuttle_stack = []                   # a list of HighwayOccupancy
        self.qubit_idx_dict = gen_qubit_idx_dict(G)
        self.idx_qubit_dict = gen_idx_qubit_dict(G)

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

    def set_highway_status(self, idx, shuttle_idx, period):
        assert idx <= len(self.highway_status_record)
        if idx == len(self.highway_status_record):
            self.highway_status_record.append((shuttle_idx, period))
        else:
            self.highway_status_record[idx] = (shuttle_idx, period)
    
    def get_highway_status(self, time):
        if time < 0 or not self.shuttle_stack:
            return (-1, None)
        if time < len(self.highway_status_record):
            return self.highway_status_record[time]
        
        last_shuttle_idx = len(self.shuttle_stack) - 1
        if self.get_shuttle_end_time(last_shuttle_idx) is None:
            return (last_shuttle_idx, 'exec')
        else:
            return (-1, None)


    def end_shuttle(self, shuttle_idx, circuit):
        # shuttle measurement start condition: after the last mop
        # This is because the Pauli correction on the control qubit in the measurement period can be propagate afterward
        exec_start_time = self.get_shuttle_exec_start_time(shuttle_idx)
        meas_start_time = circuit.mop_depth
        end_time = meas_start_time + self.meas_period - 1
        self.set_shuttle_meas_start_time(shuttle_idx, meas_start_time)
        self.set_shuttle_end_time(shuttle_idx, end_time)
        for idx in range(exec_start_time, meas_start_time):
            self.set_highway_status(idx, shuttle_idx, 'exec')
        for idx in range(meas_start_time, end_time+1):
            self.set_highway_status(idx, shuttle_idx, 'meas')

    def allocate_shuttle(self, circuit):
        # new shuttle start condition: 
        # 1. after the last shuttle
        # 2. after all gates on the control qubit and data qubits in the critical paths
        if not self.shuttle_stack:
            prep_start_idx = 0
            new_shuttle_idx = 0
        else:
            last_shuttle_idx = len(self.shuttle_stack)-1
            self.end_shuttle(last_shuttle_idx, circuit)
            last_end_idx = self.get_shuttle_end_time(last_shuttle_idx)
            prep_start_idx = max(self.get_shuttle_end_time(last_shuttle_idx) + 1, circuit.depth)  #TODO: distinguish data qubits on and off critical paths
            for idx in range(last_end_idx, prep_start_idx):
                self.set_highway_status(idx, -1, 'local')
            new_shuttle_idx = last_shuttle_idx + 1

        exec_start_time = prep_start_idx + self.prep_period
        self.shuttle_stack.append(HighwayOccupancy(self.chiplet_array))
        self.set_shuttle_prep_start_time(new_shuttle_idx, prep_start_idx)
        self.set_shuttle_exec_start_time(new_shuttle_idx, exec_start_time)
        self.set_shuttle_meas_start_time(new_shuttle_idx, None)
        self.set_shuttle_end_time(new_shuttle_idx, None)
        for idx in range(prep_start_idx, exec_start_time):
            self.set_highway_status(idx, new_shuttle_idx, 'prep')
        self.set_highway_status(exec_start_time, new_shuttle_idx, 'exec')

    def which_shuttle(self, time):
        shuttle_idx, period = self.get_highway_status(time)
        if period in {'meas', 'local'}:
            if shuttle_idx + 1 < len(self.shuttle_stack):
                return shuttle_idx + 1
            else:
                return -1
        else:
            return shuttle_idx
    
    def is_entrance_idle_at_idx(self, circuit, shuttle_idx, target_entrance, idx):
        highway_shuttle = self.shuttle_stack[shuttle_idx]
        for data_qubit in highway_shuttle.entrance_data_dict[target_entrance]:
            
            line = self.qubit_idx_dict[data_qubit]
            if circuit.take_role(line, idx) in {'mc', 'mt'}:
                return False
        return True
            
    def get_earliest_index_for_2qubit_component_on_shuttle(self, shuttle_idx, circuit, op, auto_commuting=True, auto_cancellation=False):
        # check highway path availability       
        control_qubit, target_qubit = self.idx_qubit_dict[op.control], self.idx_qubit_dict[op.target]
        highway_shuttle = self.shuttle_stack[shuttle_idx]
        if highway_shuttle.source_target_counter[(control_qubit, target_qubit)] > 0:
            target_entrance = highway_shuttle.target_entrance_dict[target_qubit]
        else:
            path = highway_shuttle.find_free_highway_path(control_qubit, target_qubit)
            if not path:
                return -1
            target_entrance = path[-2]

        # check whether it is the last shuttle
        exec_start_time = self.get_shuttle_exec_start_time(shuttle_idx)
        meas_start_time = self.get_shuttle_meas_start_time(shuttle_idx)
        if meas_start_time is not None:
            latest_idx = meas_start_time - 1
        else:
            line_depth = max(circuit.get_line_depth(op.control), circuit.get_line_depth(op.target))
            latest_idx = max(line_depth, exec_start_time)
            
        # check control data viability
        for idx in range(latest_idx, exec_start_time, -1):
            prior_role = circuit.take_role(op.control, idx - 1)
            if prior_role not in {None, 'c', 'mc'}: # control data has been changed
                print(op, 'impossible: ', shuttle_idx)
                return -1
        
        # check entrance availability
        is_entrance_idle = partial(self.is_entrance_idle_at_idx, circuit, shuttle_idx, target_entrance)
        earliest_idx = circuit.get_earliest_index_for_2qubit_component(op, auto_commuting, auto_cancellation, min_idx=exec_start_time, max_idx=latest_idx, addition_condition=is_entrance_idle)
                     
        return earliest_idx
            
    def get_highway_aware_earliest_index_for_2qubit_component(self, circuit, op, auto_commuting=True, auto_cancellation=False):
        if not self.shuttle_stack:
            return -1
        highway_agnostic_earliest_idx = circuit.get_earliest_index_for_2qubit_component(op, auto_commuting=auto_commuting, auto_cancellation=auto_cancellation)
        earliest_shuttle_idx = self.which_shuttle(highway_agnostic_earliest_idx)
        for shuttle_idx in range(earliest_shuttle_idx, len(self.shuttle_stack)):
            exec_idx = self.get_earliest_index_for_2qubit_component_on_shuttle(shuttle_idx, circuit, op, auto_commuting=auto_commuting, auto_cancellation=auto_cancellation)
            if exec_idx >= 0:
                return exec_idx
        return -1
    
    def execute_on_highway(self, circuit, op, shuttle_idx=None, exec_idx=None, auto_commuting=True, auto_cancellation=False):
        if shuttle_idx is None:
            if exec_idx is None:
                exec_idx = self.get_highway_aware_earliest_index_for_2qubit_component(circuit, op, auto_commuting=auto_commuting, auto_cancellation=auto_cancellation)
            shuttle_idx = self.which_shuttle(exec_idx)
        if shuttle_idx < 0:
            self.allocate_shuttle(circuit)
            shuttle_idx = len(self.shuttle_stack) - 1
            exec_idx = self.get_earliest_index_for_2qubit_component_on_shuttle(shuttle_idx, circuit, op, auto_commuting=auto_commuting, auto_cancellation=auto_cancellation)

        highway_shuttle = self.shuttle_stack[shuttle_idx]
        
        control_qubit, target_qubit = self.idx_qubit_dict[op.control], self.idx_qubit_dict[op.target]
        if highway_shuttle.source_target_counter[(control_qubit, target_qubit)] > 0:
            highway_shuttle.source_target_counter[(control_qubit, target_qubit)] += 1
        else:
            path = highway_shuttle.find_free_highway_path(control_qubit, target_qubit)
            highway_shuttle.occupy_highway_path(control_qubit, path)
        depth = exec_idx + 1
        circuit.add_mop_2qubit_component(op, depth, auto_commuting=auto_commuting, auto_cancellation=auto_cancellation)

    def get_highway_aware_earliest_index_for_1qubit_op(self, circuit, op):
        # highway aware: local gates should avoid occupying highway preparation periods 
        line_depth = circuit.get_line_depth(op.q)
        
        earliest_idx = -1
        for idx in range(line_depth, -1, -1):
            if len(self.shuttle_stack) == 0:
                is_in_prep_period = False
            else:
                is_in_prep_period = self.get_highway_status(idx)[1] == 'prep'
            if circuit.is_position_empty(op.q, idx) and not is_in_prep_period: #TODO: distinguish data qubits on and off critical paths
                earliest_idx = idx
        return earliest_idx

    def get_highway_aware_earliest_index_for_2qubit_op(self, circuit, op, auto_commuting, auto_cancellation):
        # highway aware: local gates should avoid occupying highway preparation periods 
        control, target = op.control, op.target
        line_depth = max(circuit.get_line_depth(control), circuit.get_line_depth(target))
        def is_not_in_prep_period(idx): #TODO: distinguish data qubits on and off critical paths
                if len(self.shuttle_stack) == 0:
                    return True
                else:
                    return self.get_highway_status(idx)[1] != 'prep'
        earliest_idx = circuit.get_earliest_index_for_2qubit_op(op, auto_commuting, auto_cancellation, min_idx=0, max_idx=line_depth, addition_condition=is_not_in_prep_period)
        return earliest_idx
            
    def execute_on_local(self, circuit, op, auto_commuting=True, auto_cancellation=False):
        if op.qubit_num == 1:
            earliest_idx = self.get_highway_aware_earliest_index_for_1qubit_op(circuit, op)
            if earliest_idx < 0:
                earliest_idx = self.get_shuttle_exec_start_time(len(self.shuttle_stack)-1)
            circuit.add_1qubit_op(op, depth=earliest_idx + 1)
        else:
            earliest_idx = self.get_highway_aware_earliest_index_for_2qubit_op(circuit, op, auto_commuting=auto_commuting, auto_cancellation=auto_cancellation)
            if earliest_idx < 0:
                earliest_idx = self.get_shuttle_exec_start_time(len(self.shuttle_stack)-1)
            circuit.add_2qubit_op(op, depth=earliest_idx + 1)
        