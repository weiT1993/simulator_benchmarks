from utils.MIP_searcher import find_cuts
from utils.cutter import cut_circuit
from qcg.generators import gen_qft

# import qiskit.circuit.library as library
from qiskit.converters import circuit_to_dag, dag_to_circuit
from qiskit.circuit.quantumregister import QuantumRegister
from qiskit import QuantumCircuit
import pickle
import matplotlib.pyplot as plt
import numpy as np
import math
import os
import networkx as nx

num_qubits = 5
approximation_degree=int(math.log(num_qubits,2))
print(approximation_degree)
# circ = library.HiddenLinearFunction(adjacency_matrix= [[1, 1, 0], [1, 0, 1], [0, 1, 1]])
# circ = library.PauliFeatureMap(feature_dimension=5,reps=1)
# circ = library.QFT(num_qubits=num_qubits,approximation_degree=num_qubits-approximation_degree,do_swaps=False)
# print(circ)

circ = gen_qft(width=num_qubits,approximation_degree=approximation_degree,barriers=False)
print(circ)

# circ = library.Permutation(num_qubits=6,pattern=[2,3,4,5,0,1]).decompose()
# pickle.dump(circ, open('circ_eg.pckl','wb'))
# dag = circuit_to_dag(circ)
# dag.draw()
# for node in dag._multi_graph.nodes():
#     data = dag._multi_graph.get_node_data(node._node_id)
#     print(type(data),data)
# dag.draw()
# print('%d two-qubit gates'%circ.num_nonlocal_gates())

solution_dict = find_cuts(circ=circ,max_cluster_qubit=int(num_qubits/1.5))
solution_dict['model'].print_stat()
cluster_circs, complete_path_map, K, d = cut_circuit(circ=circ,positions=solution_dict['positions'])
print('MIP: {:d}-qubit full circuit {:d} cuts, clusters: {}'.format(num_qubits,len(solution_dict['positions']),solution_dict['num_d_qubits']))
print('Cutter: {}'.format(d))

# approximation_degree = np.array(range(num_qubits))
# prob = [8/math.pi**2*math.sin(math.pi/4*(num_qubits-m)/num_qubits)**2 for m in approximation_degree]

# plt.figure()
# plt.plot(approximation_degree,prob)
# plt.axhline(y=4/math.pi**2,color='k',linestyle='--')
# plt.xlabel('approximation degree')
# plt.ylabel('prob')
# plt.yscale('log')
# plt.title('%d-qubit AQFT'%num_qubits)
# plt.show()