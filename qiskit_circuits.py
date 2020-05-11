from qiskit.circuit.library.arithmetic import WeightedAdder
import qiskit.circuit.library as library
from qiskit.converters import circuit_to_dag, dag_to_circuit
from qiskit.circuit.library.standard_gates import XGate
import pickle

circ = library.HiddenLinearFunction(adjacency_matrix= [[1, 1, 0], [1, 0, 1], [0, 1, 1]])
circ = library.PauliFeatureMap(feature_dimension=5,reps=1)
circ = library.Permutation(num_qubits=5,pattern=[2,4,3,0,1])
# pickle.dump(circ, open('circ_eg.pckl','wb'))
dag = circuit_to_dag(circ)
dag.draw()
print(circ)