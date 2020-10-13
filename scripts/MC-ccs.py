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

def ccs_converge(ccs, energy, nlst, itr=1000):
    """
    Args:
      ccs (list, float): array of ccs values for a single molecule, multiple conformers 
      energy (list, float): array of energy values mapping exactly to the ccs array
      nlst (list, int): array of sample sizes, e.g. [x for x in range(50, 50001, 50)]
      method (str): Selection and/or averaging method. Options are 
                    `bolt` (Boltzmann weighted average)
                    `minen` (minimum energy)
                    `minenccs` (ccs corresponding to the minimum energy)
                    `simple` (simple average)
      itr (int): Number of Monte Carlo iterations

    Returns:
      List of boltzmann weighted ccs averaged for 10,000 iterations for every n in n_list.
      List of standard deviations for each average
    """
    np.random.seed(1) #fix seed for reproducibility
    
    # Initialize dataframe
    df = pd.DataFrame(columns=['Boltzmann Weighted', 'bw std', 'Lowest Energy', 'le std', 
                               'Lowest Energy CCS', 'lec std', 'Simple Ave', 'sa std']) # ,
                               #'Ave Low Threshold', 'alt std'])

    for j, n in enumerate(nlst):
        if n % 10000 == 0:
            print(n)

        b = []
        e = []
        m = []
        s = []
        #t = []
        for i in range(0, itr):
            idx = np.random.choice(range(len(ccs)), n, replace=False)
            ccs_sample = ccs[idx]
            energy_sample = energy[idx]

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
    

def file_handle(molid, dfpath, writepath='', ssize=50, itr=1000):
    """
    """
    df = pd.read_csv(dfpath)
    ccs = df['ccs'].values
    energy = df['dft_energy'].values
    
    # Build sample array
    popsize = len(ccs)
    nlst = [x for x in range(ssize, popsize, ssize)]
    nlst.extend([popsize])
    print('tail nlst: ', nlst[-2:])

    df = ccs_converge(ccs, energy, nlst, itr)

    # Write to .csv
    df['nlst'] = nlst
    df.to_csv(join(writepath, f'{molid}_MC_ccs_iter{itr}.csv'), index=False)
                      


if __name__ == '__main__':
                  
    import time
    starttime = time.time()

    parser = argparse.ArgumentParser()
    parser.add_argument('molid', type=str, help='Molecular identifier, e.g. \'molid03\'')
    parser.add_argument('-s', '--ss', type=int, default='50', help='sample size step for building nlst, e.g. 50')
    parser.add_argument('-i', '--itr', type=int, default=1000, help='Number of Monte Carlo iterations')
    parser.add_argument('-d', '--dfpath', type=str, help='path to .csv containing conformer ccs and energy values')
    parser.add_argument('-w', '--writepath', type=str, default='.', help='write path')

    args = parser.parse_args()


    file_handle(args.molid, args.dfpath, args.writepath, itr=args.itr)

    print((time.time()-starttime)/3600, ' hours')
