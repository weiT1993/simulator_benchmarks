import quimb as qu
import quimb.tensor as qtn
import random
from tqdm import tqdm
from time import time
import matplotlib.pyplot as plt
import pickle
import os

def plot_runtime(runtimes,num_threads):
    min_full_circ_size = 100
    max_full_circ_size = 0
    plt.figure()
    for num_threads in runtimes:
        for circuit_type in runtimes[num_threads]:
            full_circ_sizes = list(runtimes[num_threads][circuit_type].keys())
            times = list(runtimes[num_threads][circuit_type].values())
            full_circ_sizes, times = zip(*sorted(zip(full_circ_sizes,times)))
            plt.plot(full_circ_sizes,times,'*-',label='%s, %d-thread'%(circuit_type,num_threads))
            min_full_circ_size = min(min(full_circ_sizes),min_full_circ_size)
            max_full_circ_size = max(max(full_circ_sizes),max_full_circ_size)
    plt.xlabel('Number of qubits',fontsize=15)
    plt.xticks(range(min_full_circ_size,max_full_circ_size+1,2),fontsize=15)
    plt.ylabel('Runtime (s)',fontsize=15)
    plt.yticks(fontsize=15)
    plt.legend(fontsize=15)
    plt.title('QUIMB runtime benchmark',fontsize=15)
    plt.tight_layout()
    plt.savefig('./quimb/quimb_runtime.pdf',dpi=400)
    plt.close()

if __name__ == '__main__':
    runtimes = {}
    for num_threads in [1]:
        os.environ['QUIMB_NUM_THREAD_WORKERS'] = '%d'%num_threads
        os.environ['QUIMB_NUM_PROCS'] = '%d'%num_threads
        os.environ['OMP_NUM_THREADS'] = '%d'%num_threads
        runtimes[num_threads] = {}
        # for circuit_type in ['supremacy_linear', 'supremacy_grid', 'bv']:
        for circuit_type in ['bv']:
            runtimes[num_threads][circuit_type] = {}
            for full_circ_size in range(20,28):
                if not os.path.isfile('./circuits/%s_q%d'%(circuit_type,full_circ_size)):
                    continue
                circ = qtn.Circuit.from_qasm_file('./circuits/%s_q%d'%(circuit_type,full_circ_size), tags='PSI_i')
                print(circ)
                # circ.psi.graph(color=['H', 'RX','RY','RZ','CZ','T'])

                begin = time()
                psi = circ.psi  # get the tensor network describing the full wavefunction
                psi.full_simplify_()  # try and simplify the network
                print(type(psi))
                T = psi.contract(all, optimize='random-greedy', output_inds=psi.site_inds)
                T.to_dense(psi.site_inds)
                elapsed = time() - begin

                output = T.data.flatten()
                print(output.shape)
                print('%s: %d-threads %d-qubits took %f seconds'%(circuit_type,num_threads,full_circ_size,elapsed))
                runtimes[num_threads][circuit_type][full_circ_size] = elapsed
                plot_runtime(runtimes=runtimes,num_threads=num_threads)
                pickle.dump(runtimes, open('./quimb/quimb_benchmark.pckl','wb'))