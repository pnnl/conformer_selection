"""Calculates and averages pairwise RMSD in groups of 50 for either across or within ("per") cycles,
as used in the Monte Carlo convergence simulation for RMSD-based geometry variability analysis.

Assumes conformers are represented as .xyz files in a directory named xyz/
"""

from isicle.scripts import rmsd
import glob
from os.path import *
import numpy as np
from numpy import random


# molid02+H1_18927493.xyz
# moid02+H1000


def percycle_rmsd(molid, annealcycles):
    """Compares conformer geometries within a cycle for the specified number of cycles.
    
    Args:
        molid (string): molecular identification as given on file name.
        annealcycles (int): number of annealing cycles.
    Returns:
        list of floats: RMSD between conformer geometries per cycle.
    """
    rmsd_percycle = []
    for cycle in range(1, annealcycles+1):
        confs = glob.glob(abspath(join('xyz','{}{}_*.xyz'.format(molid, cycle))))
        confs.sort()
        confs = confs[:50] #take only the first 50 geoms        
        rmsd_percycle.append(pw_rmsd(confs))
    return rmsd_percycle


def pw_rmsd(mols):
    """Calculates pairwise RMSD between molecular geometries.

    Args:
        mols (list): paths to molecular geometry xyz files 
    Returns:
        average pairwise RMSD of the geometries.
    """
    m = len(mols)
    k = 0
    pw = []
    for mol1 in mols:
        k += 1
        if k > m:
            break
        for i in range(k, m):
            mol2 = mols[i]
            pw.append(rmsd.rmsd(mol1, mol2))
    ave_rmsd = np.mean(pw)
    return ave_rmsd

def acrosscycle_rmsd(molid, annealcycles, samplesize):
    """Compares conformers geometries across annealing cycles using 
    Monte Carlo methods (random sampling, and in this case without replacement). 

    Args:
        molid (string): molecular identification as given on file name.
        annealcycles (int): total number of annealing cycles.
        samplesize (int): number of cycles to randomly sample.
    Returns:
        list of floats: RMSD between conformer geometries across cycles. 
    """
    random.seed(0) # fix seed for reproducibility
    rmsd_acrosscycle = []
    
    for i in range(0, annealcycles):
        # sample without replacement
        cycles = random.choice(range(1, annealcycles+1), size=samplesize, replace=False)
        sample = []
        rlst = []
        for cycle in cycles:
            confs = glob.glob(abspath(join('xyz', '{}{}_*.xyz'.format(molid, cycle))))
            #confs.sort()
            #confs = confs[:50] #take only the first 50 geoms            
            r = random.randint(0, len(confs))
            while r in rlst: # sample random without replacement
                r = random.randint(0, len(confs))
            rlst.append(r)
            sample.append(confs[r])
        rlst = []

        rmsd_acrosscycle.append(pw_rmsd(sample))
    return rmsd_acrosscycle


if __name__ == '__main__':
    from time import process_time
    start = process_time()

    molid = 
    annealcycles = 1000 #1000
    samplesize = 50 #50
    cycle_type = 'perc' #perc or acrossc

    if cycle_type == 'acrossc':
        rmsd = acrosscycle_rmsd(molid, annealcycles, samplesize) # across cycles
    if cycle_type == 'perc':
        rmsd = percycle_rmsd(molid, annealcycles) # per cycles

    std = np.std(rmsd)
    ave = np.mean(rmsd)
    mx = max(rmsd)
    sam = [std, ave, mx]
    ext = '_sorted'
    # Write to txt
    with open('{}_rmsd_{}{}.txt'.format(molid, cycle_type, ext), 'w') as f:
        for item in rmsd:
            f.write('%s\n' % item)

    with open('{}_rmsd_{}_sam{}.txt'.format(molid, cycle_type, ext), 'w') as f:
        for item in sam:
            f.write('%s\n' % item)

    with open('time_SA_rmsd_{}{}.txt'.format(cycle_type, ext), 'w') as f:
        f.write(str(process_time()-start))

    print(process_time()-start)
