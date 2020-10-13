
'''Generates convergence data of conformer RMSD or CCS at a single Monte Carlo 
step size--single point convergence.
'''

import numpy as np
import pandas as pd
from time import time
from numpy import random
from os.path import join
import argparse

def spConvergence(rmsds, S=500, mx_iter=100000):
    """Generates convergence data of conformer RMSD or CCS at a single Monte Carlo 
    step size--single point convergence.

    Args:
       rmsds (np.array):  list of conformer RMSDs for a molecule
       S (int): Sample size, default 500
       mx_iter (int): maximum iteration
    """
    random.seed(2) # Set seed for reproducibility 
    
    ave = []
    std = []
    iters = []
    simple_ave = []
    
    for iter in range(1, mx_iter+1):
        if iter % 10000 == 0:
            print('iter', iter)

        sample = random.choice(rmsds, S, replace=False)
        simple_ave.append(np.mean(sample))

        iters.append(iter)
        ave.append(np.mean(simple_ave))
        std.append(np.std(simple_ave))
        
    df = pd.DataFrame([iters, ave, std], index=['Iteration', 'Average', 'Standard Dev']).T
    return df
        
if __name__ == '__main__':
    start = time()

    parser = argparse.ArgumentParser()
    parser.add_argument('path', type=str, help='path to .txt file of conformer RMSD previously selected and averaged')
    parser.add_argument('molid', type=str, help='Molecular identifier, e.g. \'molid03\'')
    parser.add_argument('-s', '--ss', type=int, default='500', help='sample size step for building nlst, e.g. 50')
    parser.add_argument('-i', '--itr', type=int, default=100000, help='Number of Monte Carlo iterations')
    parser.add_argument('-w', '--writepath', type=str, default='.', help='write path')

    args = parser.parse_args()


    # `path` is the path to conformer RMSD previously randomly 
     # selected across or within cycles in groups of 50 and averaged
    molid = args.molid
    dr = args.path
    writepath = args.writepath
    rmsds = pd.read_csv(join(dr, f'{molid}_rmsd_acrossc_50.txt'), header=None)[0].values
    S = args.ss
    mx_iter = args.itr

    df = spConvergence(rmsds, S, mx_iter)
        
    df.to_csv(join(writepath, f'{molid}_single_point_MCconverg_S{S}_iter{mx_iter}.csv'), index=False)

    print((time()-start)/60, ' min')
