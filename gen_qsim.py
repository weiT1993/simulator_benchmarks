import os
from utils.helper_fun import generate_circ
from qiskit.converters import circuit_to_dag

def write_circuit(txt_file,circ_size,nodes):
    txt_file.write('%d\n'%circ_size)
    for node_ctr in nodes:
        node = nodes[node_ctr]
        txt_file.write('%d '%node['time_step'])
        txt_file.write('%s '%node['op_name'])
        for i, qubit in enumerate(node['qargs']):
            if i != len(node['qargs'])-1:
                txt_file.write('%d '%qubit)
            else:
                txt_file.write('%d'%qubit)
        if len(node['params'])>0:
            txt_file.write(' %f\n'%node['params'][0])
        else:
            txt_file.write('\n')

def order_circuit(circuit):
    dag = circuit_to_dag(circuit)
    all_nodes = {}
    node_ctr = 0
    for vertex in dag.topological_op_nodes():
        op_name = vertex.op.name
        params = vertex.op.params
        qargs = [x.index for x in vertex.qargs]
        all_nodes[node_ctr] = {'op_name':vertex.op.name if vertex.op.name!='tdg' else 't',
        'params':vertex.op.params,
        'qargs':qargs,
        'assigned':0}
        node_ctr += 1
    
    qubit_time_step = {x:0 for x in range(circuit.num_qubits)}
    for node_ctr in all_nodes:
        node = all_nodes[node_ctr]
        time_step = 0
        for qubit in node['qargs']:
            time_step = max(time_step,qubit_time_step[qubit])
        for qubit in node['qargs']:
            qubit_time_step[qubit] = time_step+1
        all_nodes[node_ctr]['time_step'] = time_step
    return all_nodes

if __name__ == '__main__':
    dirname = './circuits'
    if not os.path.exists(dirname):
        os.makedirs(dirname)

    for circuit_type in ['supremacy_linear', 'supremacy_grid', 'hwea', 'bv', 'adder']:
        for full_circ_size in range(5,31):
            circ = generate_circ(full_circ_size=full_circ_size,circuit_type=circuit_type)
            if circ.size()>0:
                all_nodes = order_circuit(circuit=circ)
                txt_file = open('./circuits/%s_q%d'%(circuit_type,full_circ_size),'w')
                write_circuit(txt_file=txt_file,circ_size=full_circ_size,nodes=all_nodes)
                txt_file.close()

            # txt_file = open('./circuits/bitstrings_q%d'%(full_circ_size),'w')
            # for state in range(2**25):
            #     bin_state = bin(state)[2:].zfill(full_circ_size)
            #     txt_file.write('%s\n'%bin_state)
            # txt_file.close()