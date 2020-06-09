from utils.MIP_searcher import find_cuts
from utils.cutter import cut_circuit
from qcg.generators import gen_qft
from utils.helper_fun import get_evaluator_info
from utils.passes import noiseadaptive
from utils.metrics import reliability_score

import qiskit.circuit.library as library
from qiskit.converters import circuit_to_dag, dag_to_circuit
from qiskit.circuit.quantumregister import QuantumRegister
from qiskit import QuantumCircuit
import pickle
import matplotlib.pyplot as plt
import numpy as np
import math
import os
import networkx as nx
from qiskit.compiler import transpile

num_qubits = 10

approximation_degree=int(math.log(num_qubits,2)+2)
circ = library.HiddenLinearFunction(adjacency_matrix= [[1, 1, 0], [1, 0, 1], [0, 1, 1]])
circ = library.PauliFeatureMap(feature_dimension=5,reps=1)
circ = library.QFT(num_qubits=num_qubits,approximation_degree=num_qubits-approximation_degree,do_swaps=False)
circ = library.QuantumVolume(num_qubits=num_qubits)

device_evaluator_info = get_evaluator_info(circ=circ,device_name='ibmq_johannesburg',fields=['device','initial_layout','properties','basis_gates','coupling_map'])
mapped_circuit = transpile(circ,backend=device_evaluator_info['device'],initial_layout=device_evaluator_info['initial_layout'])
reliability = reliability_score(circuit=mapped_circuit,backend_prop=device_evaluator_info['properties'])
print(reliability)

na_circ = noiseadaptive(circ,device_evaluator_info)
reliability = reliability_score(circuit=na_circ,backend_prop=device_evaluator_info['properties'])
print(reliability)

decomposed_circ = transpile(circ,basis_gates=device_evaluator_info['basis_gates'])
solution_dict = find_cuts(circ=decomposed_circ, max_cluster_qubit=8)
if solution_dict != {}:
    solution_dict['model'].print_stat()