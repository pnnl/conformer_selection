"""Runs Monte Carlo convergence simulation for conformer geometry RMSD.
"""

from os.path import *
import numpy as np
from numpy import random
import argparse

def mcConvergence(sample_rmsdave, runs):
    """
      Args:
        sample_rmsdave (list): List of rmsd averages developed from sample pairwise rmsd values between conformer 
                              or molecule geometries.
        runs (int): number of simulation runs (iterations) to do for every sample size
    """
    random.seed(2)
    stdev_ave = []
    max_ave = []
    ave_stdev = []
    ave_max = []
    ave_ave = []
    #collect_aves = []
    
    N = len(sample_rmsdave) #number of samples or number of annealing cycles (population size)
    for i in range(1, N+1):
        #stdevs = []
        aves = []
        #maxes = []
        
        for r in range(0, runs):
            item = random.choice(sample_rmsdave, i, replace=False)
            #stdevs.append(np.std(item))
            aves.append(np.mean(item))
            #maxes.append(max(item))
        #collect_aves.append(aves)
        stdev_ave.append(np.std(aves))
        #ave_stdev.append(np.mean(stdevs))
        ave_ave.append(np.mean(aves))
        #ave_max.append(np.mean(maxes))
        max_ave.append(max(aves))
    #return stdev_ave, max_ave, ave_stdev, ave_max, ave_ave, collect_aves
    return stdev_ave, max_ave, ave_ave

if __name__ == '__main__':

    import time
    starttime = time.time()

    parser = argparse.ArgumentParser()
    parser.add_argument('txtpath', type=str, help='path to .txt file of conformer RMSD previously selected and averaged')
    parser.add_argument('molid', type=str, help='Molecular identifier, e.g. \'molid03\'')
    parser.add_argument('-w', '--writepath', type=str, default='.', help='write path')
    parser.add_argument('-i', '--itr', type=int, default=10000, help='Number of Monte Carlo iterations')

    args = parser.parse_args()

    molid = args.molid
    path = args.txtpath

    rmsds = []
    with open(path, 'r') as f:
        for rmsd in f.read().splitlines():
            rmsds.append(float(rmsd))

    stdev_ave, max_ave, ave_ave = mcConvergence(rmsds, args.itr)

    # This code is so old... Why did I right it this way? Why?
    with open(join(dir, 'molid{}_ave_ave{}.txt'.format(molid, extw)), 'w') as f:
        for item in ave_ave:
            f.write('%s\n' % item)
    with open(join(dir, 'molid{}_stdev_ave{}.txt'.format(molid, extw)), 'w') as f:
        for item in stdev_ave:
            f.write('%s\n' % item)
    with open(join(dir, 'molid{}_max_ave{}.txt'.format(molid, extw)), 'w') as f:
        for item in max_ave:
            f.write('%s\n' % item)

    print((time.time()-starttime)/3600, ' hours in, time time')
