import numpy as np
from qiskit.converters import circuit_to_dag
import copy
from time import time
from scipy.special import comb
import math
from itertools import combinations

def chi2_distance(target,obs):
    assert len(target)==len(obs)
    obs = np.absolute(obs)
    obs = obs / sum(obs)
    # assert abs(sum(obs)-1)<1e-5
    distance = 0
    for x,y in zip(target,obs):
        if x-y==0:
            distance += 0
        else:
            distance += np.power(x-y,2)/(x+y)
    # distance /= len(target)
    return distance

def fidelity(target,obs):
    assert len(target)==len(obs)
    epsilon = 1e-20
    obs = np.absolute(obs)
    obs = obs / sum(obs)
    # assert abs(sum(obs)-1)<1e-5
    fidelity = 0
    for t,o in zip(target,obs):
        if t > 1e-16:
            fidelity += o
    return fidelity

def reliability_score(circuit,backend_prop):
    gate_reliab = {}
    for ginfo in backend_prop.gates:
        for item in ginfo.parameters:
            if item.name == 'gate_error':
                g_reliab = 1.0 - item.value
            elif item.name == 'gate_length':
                g_len = item.value
        if len(ginfo.qubits) not in gate_reliab:
            gate_reliab[len(ginfo.qubits)] = {}
        if ginfo.gate not in gate_reliab[len(ginfo.qubits)]:
            gate_reliab[len(ginfo.qubits)][ginfo.gate] = {}
            gate_reliab[len(ginfo.qubits)][ginfo.gate][tuple(ginfo.qubits)] = g_reliab
        else:
            gate_reliab[len(ginfo.qubits)][ginfo.gate][tuple(ginfo.qubits)] = g_reliab
    for num_qubit in gate_reliab:
        if num_qubit<=1:
            continue
        for gate_type in gate_reliab[num_qubit]:
            for edge in gate_reliab[num_qubit][gate_type]:
                assert tuple(list(edge)[::-1]) in gate_reliab[num_qubit][gate_type]

    # [print(x,gate_reliab[x]) for x in gate_reliab]
    reliability = 1
    dag = circuit_to_dag(circuit)
    for vertex in dag.topological_op_nodes():
        qargs = [x.index for x in vertex.qargs]
        reliability *= gate_reliab[len(qargs)][vertex.name][tuple(qargs)]
    return reliability

def count_dirty_qubits(circuit):
    dag = circuit_to_dag(circuit)
    dirty_qubits = set()
    for vertex in dag.topological_op_nodes():
        qargs = [x.index for x in vertex.qargs]
        for qarg in qargs:
            dirty_qubits.add(qarg)
    return len(list(dirty_qubits))

def count_qubit_mappings(coupling_map,circ_size):
    def count(graph,curr_subset,remaining_size,dirty_subsets,mappings):
        # print('Looking for curr_subset = {}, remaining_size = {:d}'.format(curr_subset,remaining_size))
        if tuple(curr_subset) in dirty_subsets:
            # print('Already explored')
            return mappings, dirty_subsets
        elif remaining_size==0:
            # print('Found',curr_subset)
            mappings.append(tuple(curr_subset))
            dirty_subsets.add(tuple(curr_subset))
            return mappings, dirty_subsets
        else:
            if len(curr_subset)==0:
                candidates = set(graph.keys())
            else:
                candidates = set()
                for vertex in curr_subset:
                    children = graph[vertex]
                    for child in children:
                        if child not in curr_subset:
                            candidates.add(child)
            for candidate in candidates:
                new_subset = copy.deepcopy(curr_subset)
                new_subset.add(candidate)
                mappings, dirty_subsets = count(graph=graph,curr_subset=new_subset,remaining_size=remaining_size-1,dirty_subsets=dirty_subsets,mappings=mappings)
            dirty_subsets.add(tuple(curr_subset))
            return mappings, dirty_subsets

    device_size = coupling_map.size()
    # print(device_size)
    edges = coupling_map.get_edges()
    graph = {}
    for edge in edges:
        src, dest = edge
        if src in graph:
            graph[src].add(dest)
        else:
            graph[src] = set()
            graph[src].add(dest)
    # [print(x,graph[x]) for x in graph]
    count_begin = time()
    mappings, dirty_subsets = count(graph=graph,curr_subset=set(),remaining_size=circ_size,dirty_subsets=set(),mappings=[])
    num_qubit_mappings = len(mappings)
    # print(mappings)
    elapsed = time()-count_begin
    print('%d-on-%d problem : found %d mappings in %f seconds.'%(circ_size,device_size,num_qubit_mappings,elapsed))
    return num_qubit_mappings