import numpy as np
import pandas as pd
from statsmodels.stats.weightstats import DescrStatsW
from os.path import *
import argparse

# bolt_ccs is modified from ISiCLE
def bolt_ccs(ccs_sample, energy_sample):
    g = energy_sample * 627.503
    mn = g.min()
    relG = g - mn
    b = np.exp(-relG / 0.5924847535)
    w = (b / b.sum()) * len(b)

    ws = DescrStatsW(ccs_sample, weights=w, ddof=0)
    return ws.mean

def threshold_ave(ccs_sample, energy_sample, threshold=5):
    """Returns average CCS of conformations associated with energies below a threshold. 
       Note the threshold is a relative to the lowest energy.
    Args:
      ccs_sample (list): CCS corresponding to the energy sample.
      energy_sample (list): units hartree/mol, which are converted to kcal/mol
                             and taken relative to the minimum energy which is set to 0.
      threshold (float): units kcal/mol. Only energies below this threshold will be averaged.
    """
    # Convert Hartree/mol to kcal/mol
    kcal = np.array(energy_sample * 627.509469)

    # Make everything relative to the lowest energy
    kmin = np.min(kcal)
    kcal -= kmin

    # Grab indices for energies less than threshold
    below = [i for i, x in enumerate(kcal[:]) if x < threshold]
    
    # Grab corresponding CCS
    t = np.array(ccs_sample)[below]

    # t = np.mean([i for i, x in enumerate((energy_sample/(9.597*10**20))[:]) if x < threshold]) # One liner
    return np.mean(t)

def ccs_converge(dfarr, nlst, iter=1000):
    """
    Args:
      dfarr (pandas.DataFrame): array of ccs and energy dataframe split up by cycle (in numerical order)
      nlst (list, int): array of sample sizes, e.g. [x for x in range(50, 50001, 50)]
      iter (int): Number of Monte Carlo iterations

    Returns:
    """
    np.random.seed(0) #fix seed for reproducibility
    
    # Initialize dataframe
    df = pd.DataFrame(columns=['Boltzmann Weighted', 'bw std', 'Lowest Energy', 'le std', 
                               'Lowest Energy CCS', 'lec std', 'Simple Ave', 'sa std']) # ,
                               #'Ave Low Threshold', 'alt std'])

    for j, n in enumerate(nlst):
        if n % 100 == 0:
            print(n)

        b = []
        e = []
        m = []
        s = []
        #t = []
        for i in range(0, iter):
            idx = []
            cyc = np.random.choice(range(0, 1000), n, replace=False)


            ccs_sample = []
            energy_sample = []
            for c in cyc:
                ccs_sample.extend(dfarr[c]['ccs'].values)
                energy_sample.extend(dfarr[c]['dft_energy'].values)

            ccs_sample = np.array(ccs_sample)
            energy_sample = np.array(energy_sample)
            

            b.append(bolt_ccs(ccs_sample, energy_sample))
            e.append(np.min(energy_sample))
            m.append(ccs_sample[np.argmin(energy_sample)])
            s.append(np.mean(ccs_sample))
            #t.append(threshold_ave(ccs_sample, energy_sample))

        row = []
        for method in [b, e, m, s]: #, t]:
            row.append(np.mean(method))
            row.append(np.std(method))
        df.loc[j] = row

    return df
    

def file_handle(molid, dfpath, writepath='', ssize=50, iter=1000):
    """
    """
    df = pd.read_csv(dfpath)
    ccs = df['ccs'].values
    energy = df['dft_energy'].values
    
    # Build sample array for per cycle. Split the dataframe up by cycles and put these sub-dataframes into an array
    nlst = [x for x in range(1, 1001)]
    dfarray = []
    for c in nlst:
        dfarray.append(df.loc[df['cycle'] == c])

    wdf = ccs_converge(dfarray, nlst, iter)

    # Write to .csv
    wdf['nlst'] = nlst
    wdf.to_csv(join(writepath, f'{molid}_MC_ccs_iter{iter}_percycle.csv'), index=False)
                      


if __name__ == '__main__':

    import time
    starttime = time.time()

    parser = argparse.ArgumentParser()
    parser.add_argument('dfpath', type=str, help='path to .csv containing conformer ccs and energy values')
    parser.add_argument('molid', type=str, help='Molecular identifier, e.g. \'molid03\'')
    parser.add_argument('-w', '--writepath', type=str, default='.' help='write path')
    parser.add_argument('-s', '--ss', type=int, default='50', help='sample size step for building nlst, e.g. 50')
    parser.add_argument('-i', '--itr', type=int, default=1000, help='Number of Monte Carlo iterations')

    args = parser.parse_args()

    file_handle(args.molid, args.dfpath, args.writepath, iter=args.itr)

    print((time.time()-starttime)/3600, ' hours')
