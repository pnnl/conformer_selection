"""Similarity Down-selection (SDS) heuristically finds the most dissimilar set of size `n`
out of a population represented as an NxN matrix where the ith row (as well as the
ith column) belong to item i, and the element (i,j) is some pairwise relation between 
items i and j. 

The pairwise relation is a floating point value where larger values indicate
greater dissimilarity, and where the pairwise relation between item i and itself is represented
as np.nan (this allows the program to work with numpy log conversion and sums).

SDS was originally written to find the subset of n most dissimliar conformers 
(similarity being determined by average pairwise RMSD between atoms), and relevant
functions for that particular application can also be found in this script.

In the original implemenation finding the n most dissimlar conformers 
for 50000x50000 matricies, SDS substationally outperformed a benchmark random sampling method in both 
time and accuracy. 

The pairwise RMSD matricies used for the original implemenation can be created using pwrmsd_writer.py

Note that because SDS finds the set of size `n` by building off of set `n-1`, finding set `n` also finds
also previous set sizes from 2-n. 



Author: Felicity Nielson
2019-2020
"""

from time import time
import numpy as np
import pandas as pd
from os.path import *
import argparse
import os
from statsmodels.stats.weightstats import DescrStatsW
from math import log


def SDS(df, n=3):
    """Finds the `n` most dissimilar items. In the input matrix, the ith row (and 
    ith column) is an array belonging to item i. The matrix element (i,j) would then be the pairwise 
    dissimilarity metric between items i and j (for example, geometric RMSD between molecular conformers i and j).
    Args:
      df (pandas.DataFrame): Square matrix where each row (and by symmetry, column) is an array 
                             corresponding to a specific item or object, and each element (i,j) the 
                             floating point dissimilarity between items i and j. The element (i,i) must
                             be represented as np.nan (for log-summing).

                             If pairwise data between two                                         
                             items is missing, this can also be represented as np.nan. What will  
                             happen is the second item will automatically be set to np.nan in the 
                             log-summation array once the first of one of the two items is chosen.
                             Thus, the second item will never be chosen. In this same manner,     
                             entire missing items can be represented as arrays of np.nan,         
                             as a trick to preserve externally related indexing.                  

      n (int): Dissimilar set size to find.
               1 < n < N, where N is the full population size.
    Returns:
      pandas.DataFrame containing indices of the items found in the dissimilar set of size `n`. 
    """
    # Check df matrix is square
    N = len(df.index)
    assert N == len(df.columns)
    
    # Reduce n to maximum number of items if n is over, otherwise search will fail with an erorr.
    M = len(df.dropna(how='all'))
    if n > M:
        n = M

    print(f'Starting SDS search for most dissimilar set of size n = {n}...')

    # First grab matrix indices of the two most dissimilar geometries
    row_mx = []
    
    for i in range(N):
        row_mx.append(np.nanmax(df.loc[i]))
    ind1 = np.nanargmax(row_mx)
    ind2 = ind1 + 1 + np.nanargmax(row_mx[ind1+1:])

    # Initialize the dissimilar matrix with the two most dissimilar
    disarray = [np.array(df.loc[ind1]), np.array(df.loc[ind2])]
    indices = [ind1, ind2]
    
    # Find n-2 other most dissimilar
     # Multiply the rows of the n-1 dissimilar set. Or,
     # use log summing if N is large (e.g. 50000) to avoid 
     # exceeding floating point machine precision.
     # This script uses log summing.
     # The index of the largest value is the index of the nth
     # item which makes the nth dissimilar set.

    # Initialize array for log summing
    logsum = [0 for x in range(N)]
    logsum += np.log(disarray[0])
    
    for i in range(n-2):
        logsum += np.log(disarray[-1])
        indn = np.nanargmax(logsum)
        indices.append(indn)
        disarray.append(np.array(df.loc[indn]))

    return_df = pd.DataFrame([indices], index=['matrix index']).T

    print('Finished')
    return return_df   


def conf_to_ccs(conformers, infodf):
    """Converts list of conformer indexes (from a 50x50 matrix) to Boltzmann weighted CCS average and Lowest Energy CCS.
    Args:
      conformers (np.array): Array of conformer indexes selected from a 50x50 rmsd matrix as returned in the 
                             SDS() dataframe.

    Returns:
      df with new columns: DFT Energy, CCS, BW CCS, and Lowest Energy CCS
    """
    print('Calculating CCSs...')
    
    csvdf = infodf.loc[conformers].copy()

    columns = ['n Dissimilar', 'Conformer', 'CCS', 'DFT Energy', 'BW CCS', 'Lowest Energy CCS']
    finaldf = pd.DataFrame(columns=columns)

    narray = np.array([x for x in range(1, len(conformers)+1)])

    for n in narray:
        # Calculate Boltzmann Weighted CCS
        bw =  bw_ccs(csvdf.iloc[0:n])

        # Find CCS with the lowest energy
        #idx = csvdf.loc[1:n]['dft_energy'].idxmin()
        idx = csvdf.iloc[0:n]['dft_energy'].idxmin()
        lec = csvdf['ccs'][idx]
        finaldf = finaldf.append(pd.DataFrame([[n,conformers[n-1], csvdf['ccs'].iloc[n-1], csvdf['dft_energy'].iloc[n-1], bw, lec]],
                                              columns=columns), ignore_index=True)

    return finaldf


def bw_ccs(df):
    """Calculates the Boltzmann weighted average. Note DFT energies are
    expected to be given in Hartree/mol. This function converts to kcal/mol before weighting.
    Modified from ISiCLE
    Args:
      df (pandas.DataFrame): contains 'dft_energy' and 'ccs' columns
    Returns:
      Boltzmann weighted mean
    """

    if len(df['dft_energy']) == 1:
        return df['ccs'].values[0]

    g = df['dft_energy'].values * 627.503
    mn = g.min()
    relG = g - mn
    b = np.exp(-relG / 0.5924847535)
    w = (b / b.sum()) * len(b)

    ws = DescrStatsW(df['ccs'], weights=w, ddof=0)

    return ws.mean



if __name__ == '__main__':
    start = time()
    
    parser = argparse.ArgumentParser()
    parser.add_argument('-m', '--mtrx', type=str,
                        help='Path to .pkl file containing an NxN matrix')    
    parser.add_argument('-n', '--ndis', type=int, default=3, 
                        help='Find the n most dissimilar items')
    parser.add_argument('-c', '--conf_info', type=str, default=None,
                        help='tsv file containing ccs and energy info for conformers. \
                        Note that the index of the .tsv file must exactly match the index of the matrix')    
    
    args = parser.parse_args()
    n = args.ndis
    mtrx = args.mtrx

    # If SDS is not already a directory, make it
    directory = 'SDS'
    if not exists(directory):
        os.makedirs(directory)


    df = pd.read_pickle(mtrx)
    SDSdf = SDS(df, n=n)
    narray = np.array([x for x in range(1, n+1)])

    # If comparing conformers, calculate Boltzmann weighted CCS.
    if args.conf_info != None:
        csvdf = pd.read_csv(args.conf_info)
        writedf = conf_to_ccs(SDSdf['matrix index'].values, csvdf)
    else:
        writedf = SDSdf
        writedf['n Dissimilar'] = narray


    writedf.to_csv(f'SDS/SDS_{n}_dissimilar.csv', index=False)
    print((time()-start)/60, 'min')
        
        
