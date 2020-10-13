"""Combines individual .tsv into a single .csv file. .tsv and .csv files contain information for cycle, geometry, ccs, and dft energy
as output from ISiCLE calculated ccs files. I wrote this code when I was an infant and did not understand the power of 
pandas dataframe software for simple .tsv/.csv file type reading and conversion. So this could easily be refactored to be much better.

"""

import pandas as pd
import matplotlib.pyplot as plt
import glob
from os.path import *
import argparse


def ccs_energy_df(molid, dirpath, writepath=''):
    """Parses MOBCAL ccs output file and writes pandas dataframe for cycle, geometry, ccs, and dft energy.
       Args:
         dirpath (str): path to directory containing .tsv files of ccs, dft energy, etc values as outputed by MOBCAL through ISiCLE.
         writepath (str): path to write .csv and/or .txt files.
         csv (bool): default True. Writes energy and ccs pandas dataframe to .csv (recommended)
         txt (bool): default True. Writes energy and ccs lists to .txt (optional)
    """
    ccs = []
    energy = []
    cycle = []
    geom = []

    tsvs = glob.glob(join(dirpath, '*.tsv'))
    print('Number of  tsvs found: ', len(tsvs))
    tsvs.sort()

    for file in tsvs:
        with open(file, 'r') as f:
            split = f.read().splitlines()
            if len(split) == 1:
                print('Error in file: ', file)
            val = split[1].split('\t')
            ccs.append(val[1])
            energy.append(val[3])
        slash = file.split('/')
        ident = slash[-1].split('_')
        cycle.append(ident[2])
        ident3 = ident[3].split('.')
        ident30 = ident3[0].split('m')
        geom.append(ident30[1])

    # Write dataframe
    df = pd.DataFrame([cycle, geom, ccs, energy], index=['cycle','geometry','ccs','dft_energy']).T
    df.to_csv(join(writepath, '{}_ccs{}.csv'.format(molid, ext)), index=False)

if __name__ == '__main__':
    import time
    starttime = time.time()

    parser = argparse.ArgumentParser()
    parser.add_argument('molid', type=str, help='Molecular identifier, e.g. \'molid03\'')
    parser.add_argument('dirpath', type=str, help='path to directory containing conformer ccs .tsv files as outputed by MOBCAL in ISiCLE')
    parser.add_argument('-w', '--writepath', type=str, default='.', help='write directory path')

    args = parser.parse_args()

    molid = args.molid
    dirpath = args.dirpath
    writepath = args.writepath

    ccs_energy_df(molid, dirpath, writepath)
    
    print((time.time()-starttime)/60, ' minutes')
