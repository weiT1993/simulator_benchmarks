import pickle
from utils.helper_fun import read_file
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import numpy as np

def func(x, a, b, c):
    return a * np.exp(b * x) + c

if __name__ == '__main__':
    times = read_file(filename='./qiskit_benchmark.pckl')
    plt.figure()
    x_ticks = set()
    for circuit_type in times:
        sizes = list(times[circuit_type].keys())
        ydata = list(times[circuit_type].values())
        x_fit = np.linspace(min(sizes), 36, 50)
        [x_ticks.add(x) for x in range(int(min(x_fit)), int(max(x_fit))+1) if x%5==0]

        popt, pcov = curve_fit(func, sizes, ydata)
        print(popt)

        plt.plot(sizes,ydata, '*', label='%s'%circuit_type)
        # plt.plot(x_fit, func(x_fit, *popt), 'r-',label=r'$O(e^{%.2f \times q})$'%(popt[1]))
        plt.plot(x_fit, func(x_fit, *popt), '-')
        
    x_ticks = list(x_ticks)
    plt.axhline(y=60*60,color='k',linestyle='--')
    plt.text(x=max(x_ticks)-2, y=60*60*1.2, s='1h', color='k',fontsize=15)
    plt.axhline(y=60,color='k',linestyle='--')
    plt.text(x=max(x_ticks)-2, y=60*1.2, s='1min', color='k',fontsize=15)
    plt.xticks(x_ticks,x_ticks,fontsize=15)
    plt.yticks(fontsize=15)
    plt.yscale('log')
    plt.xlabel('Number of qubits',fontsize=15)
    plt.ylabel('Simulation Time (s)',fontsize=15)
    plt.legend(fontsize=15)
    plt.title('Classical Simulation of Quantum Circuits\nScales Exponentially',fontsize=15)
    plt.tight_layout()
    plt.savefig('./qiskit_runtime.pdf',transparent=False,dpi=400)
    plt.close()