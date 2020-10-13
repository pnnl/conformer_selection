'''Finds the convergence point of the Monte Carlo (MC) convergence simulation, e.g. when the MC standard deviation 
comes within 0.01% of the final converged value.
'''
import numpy as np
import pandas as pd
from os.path import join
import glob
import argparse

def criterion(convals, stdevs, crit, molids, method='non_normalized'):
    """Finds when convergence is met based on a given criterion. Traces backwards so no stdev values above the criterion are past the cutoff.
    Args: 
      masses (list, float): Masses of molecules
      convals (list, float): The value all sample increments are converging to. Also the value of the last increment
                             (when the full population is selected as the sample). 
                             There is a convergence value for each molecule.
      stdevs (list, str): List of stdev arrays.
      crit (float): The criterion for determining convergence. 
                    When standard deviation is within this fraction of the convergence value, we say "convergence" is met.
                    Given in decimal (not percent) form.
      molids (list, str): list of molecule identifiers
      method (str): 

    Returns:
      incr (list, int): Increment numbers corresponding to the cutoff, when the criterion is met. 
    """
    print('Criterion: {}, or {}%'.format(crit, crit*100))
    increments = []

    for i in range(0,len(stdevs)):
        incr = 0

        if method == 'forward':
            ran = [x for x in range(0,len(stdevs[i]))]
            for j in ran:
                if stdevs[i][j]/convals[i] < crit:
                    #print('length stdev: ', len(stdevs[i]))
                    #print('j: ', j)                                                                                                 
                    incr = j+1
                    #print(incr)
                    increments.append(incr)
                    break

            if incr == 0 :
                increments.append(incr)
                print('Warning: No increment value found for {}sth molecule, {}.'.format(i, molid[i]))

        elif method == 'backwards':
            ran = [-x for x in range(1, len(stdevs[i]))]
            for j in ran:
                if stdevs[i][j]/convals[i] > crit:
                    incr = len(stdevs[i])+j+1
                    increments.append(incr)
                    break

            if incr == 0 :
                increments.append(incr)
                print('Warning: No increment value ABOVE the cutoff found for {}sth molecule, {}.'.format(i, molid[i]))
            if incr == 1001:
                increments.append('N/A')
                print('Warning: No increment value BELOW the cutoff found for {}sth molecule, {}.'.format(i, molid[i]))

        elif method == 'non_normalized':
            for j, std in enumerate(stdevs[i]):
                if std < crit:
                    increments.append(j+1)
                    break                    
        else:
            print(method, ' is not an accepted method')

    return increments


def array_build(path, ext, cycle_type):
    """Builds standard deviation matrix and convergence value array as output by ave-rmsd-MC.py
    Args:
      path (str): path to directory containing standard deviation and Monte Carlo averages
      ext (str): file name extension for the stdev and MC averages
      cycle_type (str):
    Returns:
      standard deviation matrix, convergence values, max standard deviations, and molecule ids
    """
    molids = []
    stdevs = []
    convals = []
    max_stdevs = []

    for x in ['stdev', 'ave']:
        paths = glob.glob(join(path, '*_{}_ave_{}{}.txt'.format(x, cycle_type, ext))) #molid02_stdev_ave_acrossc_50.txt
        paths.sort() # make sure ids in numerical order
        for p in paths:
            with open(p, 'r') as f:
                file = f.read().splitlines()
                if x == 'stdev':
                    std = [float(line) for line in file]
                    stdevs.append(std)
                    max_stdevs.append(np.max(std))

                    slash = p.split('/')
                    score = slash[-1].split('_')
                    molids.append(score[0])
                elif x == 'ave':
                    convals.append(float(file[-1]))
                else:
                    print('Error: %s is not an option.' % x)
                    break

    return convals, stdevs, max_stdevs, molids



if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('mcdir', type=str, help='path to directory containing standard deviation and Monte Carlo averages ')
    #parser.add_argument('masspath', type=str, help='path to .txt of masses')
    #parser.add_argument('writepath', type=str, help='write path for the output')
    parser.add_argument('-c', '--crit', type=float, default='0.001', help='')
    parser.add_argument('--extr', type=str, default='', help='file name extension on ave-rmsd-MC.py files')
    parser.add_argument('--extw', type=str, default='', help='write file name extension')
    parser.add_argument('-t', '--cycle_type', type=str, default='acrossc', help='either \'acrossc\' or \'perc\' (across cycles or per cycle)')
    parser.add_argument('-m', '--method', default='non_normalized', help='\'non_normalized\', \'forward\', or \'backward\'')

    args = parser.parse_args()


    # Build standard deviation and convergence value arrays
    mcdir = args.mcdir
    convals, stdevs, max_stdevs, molids = array_build(mcdir, args.extr, args.cycle_type)
    
    masspath = join(mcdir, 'mass-v-rmsdconverg/masses.txt') #Personal
    with open(masspath, 'r') as f:
        masses = f.read().splitlines()

    # Calculate sample size for convergence
    # Note that 'sample size' reflects 50n pairwise RMSDs, where n=1-1000. 
    # So if sampsize=102, that is 102 samples of 50 pairwise RMSD averages.
    sampsize = criterion(convals, stdevs, args.crit, molids, method=args.method)

    # Save dataframe
    df = pd.DataFrame([molids, masses, sampsize, convals, max_stdevs], 
                      index=['molid', 'mass', 'RMSD converg sample size', 'RMSD converg value', 'RMSD max stdev']).T
    print("\n", df) 
    df.to_csv(join(mcdir, 'mass-v-rmsdconverg/mass-v-converg-cutoff-{}-{}-{}{}.csv'.format(args.crit, args.cycle_type, args.method, args.extw)))
