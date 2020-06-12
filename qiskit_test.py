from time import time
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import numpy as np
import pickle

from utils.helper_fun import generate_circ, evaluate_circ

def plot_runtime(runtimes):
    min_full_circ_size = 100
    max_full_circ_size = 0
    plt.figure()
    for circuit_type in runtimes:
        full_circ_sizes = list(runtimes[circuit_type].keys())
        times = list(runtimes[circuit_type].values())
        full_circ_sizes, times = zip(*sorted(zip(full_circ_sizes,times)))
        plt.plot(full_circ_sizes,times,'*-',label='%s'%circuit_type)
        min_full_circ_size = min(min(full_circ_sizes),min_full_circ_size)
        max_full_circ_size = max(max(full_circ_sizes),max_full_circ_size)
    plt.xlabel('Number of qubits',fontsize=15)
    plt.xticks(range(min_full_circ_size,max_full_circ_size+1,2),fontsize=15)
    plt.ylabel('Runtime (s)',fontsize=15)
    plt.yticks(fontsize=15)
    plt.legend(fontsize=15)
    plt.title('QISKIT runtime benchmark',fontsize=15)
    plt.tight_layout()
    plt.savefig('qiskit_runtime.pdf',dpi=400)
    plt.close()

if __name__ == '__main__':
    runtimes = {}
    for circuit_type in ['supremacy_linear', 'supremacy_grid', 'hwea', 'bv', 'adder', 'aqft']:
        runtimes[circuit_type] = {}
        for full_circ_size in range(10,29):
            if circuit_type != 'supremacy_grid' and full_circ_size%2!=0:
                continue
            circ = generate_circ(full_circ_size=full_circ_size,circuit_type=circuit_type)
            if circ.size()>0:
                begin = time()
                sv = evaluate_circ(circ=circ,backend='statevector_simulator',device_name='null')
                elapsed = time() - begin
                print('%s: %d-qubits %d-gates took %f seconds'%(circuit_type,full_circ_size,circ.size(),elapsed))
                runtimes[circuit_type][full_circ_size] = elapsed
                plot_runtime(runtimes=runtimes)
                pickle.dump(runtimes, open('qiskit_benchmark.pckl','wb'))
