from qiskit.transpiler.passes import BasicSwap, NoiseAdaptiveLayout, LookaheadSwap, StochasticSwap, CXCancellation
from qiskit.converters import circuit_to_dag, dag_to_circuit
from qiskit.compiler import transpile

def cancellations(dag):
    cx_cancellation = CXCancellation()
    dag = cx_cancellation.run(dag)
    circ = dag_to_circuit(dag)
    return circ

def basic(circ,device_info):
    basic_circ = transpile(circ,basis_gates=device_info['basis_gates'],
    coupling_map=device_info['coupling_map'],
    backend_properties=device_info['properties'],
    routing_method='basic')
    return basic_circ

def lookahead(circ,device_info):
    lookahead_circ = transpile(circ,basis_gates=device_info['basis_gates'],
    coupling_map=device_info['coupling_map'],
    backend_properties=device_info['properties'],
    routing_method='lookahead')
    return lookahead_circ

def stochastic(circ,device_info):
    stochastic_circ = transpile(circ,basis_gates=device_info['basis_gates'],
    coupling_map=device_info['coupling_map'],
    backend_properties=device_info['properties'],
    routing_method='stochastic')
    return stochastic_circ

def noiseadaptive(circ,device_info):
    noise_adaptive_circ = transpile(circ,basis_gates=device_info['basis_gates'],
    coupling_map=device_info['coupling_map'],
    backend_properties=device_info['properties'],
    layout_method='noise_adaptive')
    return noise_adaptive_circ