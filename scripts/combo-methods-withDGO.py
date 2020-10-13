import pandas as pd
import numpy as np
import argparse
import glob
from os.path import *
from statsmodels.stats.weightstats import DescrStatsW

####################################### 1

def BW(ccs_sample, energy_sample):
    """Calculates the Boltzmann weighted average CCS 
    based on given energies. Energies are expected to be in kcal/mol, relative to some minimum... Full population minimum?
    or individual sample minimum? I'm not sure... I need to look at how DescrStatsW calculates weights.
    """
    if len(energy_sample) == 1:
        return ccs_sample[0]
    
    relG = energy_sample
    b = np.exp(-relG / 0.5924847535)
    w = (b / b.sum()) * len(b)
    
    ws = DescrStatsW(ccs_sample, weights=w, ddof=0)
    return ws.mean

def LE(ccs_sample, energy_sample):
    """Finds the lowest energy conformer of a sample.
    """
    idx = np.argmin(energy_sample)
    return ccs_sample[idx]

def SA(ccs_sample):
    return np.mean(ccs_sample)


####################################### 2


def nDis(mtrx, n=3):
    """Finds the most dissimilar set of size n. Logsum friendly, np.nan friendly
    Args:
      mtrx (np.array): contains pairwise rmsd matrix where each row and colum correponds to a conformer.
      n (int): The n most dissimilar geometries
    Returns:
      indices of n most dissimilar conformers.
    """
    # Check df matrix is square
    N = len(mtrx)                
    assert N == len(mtrx[0])

    if n == 0:
        return np.array([])                                                                                        

    row_mx = np.nanmax(mtrx, axis=1)
    indices = (-row_mx).argsort()[:2]
    ind1 = indices[0]
    if n == 1:
        return np.array([ind1])
    else:
        ind2 = indices[1]


    # Initialize dissimilar RMSD-array with the two most dissimilar geometries
    disarray = [np.array(mtrx[ind1]), np.array(mtrx[ind2])]

    # Initialize array for log summing
    logsum = [0 for x in range(N)]
    logsum += np.log(disarray[0])

    for i in range(n-2):
        logsum += np.log(disarray[-1])
        indn = np.nanargmax(logsum)
        #indices.append(indn)
        indices = np.append(indices, indn)
        disarray.append(np.array(mtrx[indn]))    #(df.loc[indn]))

    return indices



####################################### 3


def combo_methods(df, times, pwRMSD): #, dropallnan=False): might replace dropallnan with something automatic ***
    '''Performs a parameter sweep of different combinations of below energy threshold (BET), 
    similarity down-selection (SDS), random selection, DFT energy calculation, and DFT geometry optimization methods,
    to down select and optimize from a larger set of conformers, and then perform Boltzmann weighting (BW),
    lowest energy (LE), and simple average (SA) at the end. The actual energy and optimizations are not calculated 
    here. This function merely counts up the time it would take to do each task, and outputs the corresponding
    collisional cross section (CCS) that would be obtained in MOBCAL if the down-selected and generated conformers
    were chosen.
    
    The goal is to provide an estimation on the time vs accuracy tradeoff.
    
    Args:
      df (pandas.DataFrame): Rows are 50,000 conformers. Columns are... 
        'MD CCS', 'MD Energy (EPtot)', 'MD Energy (Etot)', 'DFT Energy', 'DFT GO CCS', 'DFT GO Energy', '50k SDS Rank'
        
      times (pandas.DataFrame): Contains times for the following steps (these are columns in the dataframe)...
        'SDS time per pair', 'AMBER ave time', 'DFT Energy ave time', 'DFT Geom Opt ave time', 'MOBCAL ave time'
        
      pwRMSD (numpy.array, matrix): 50,000 x 50,000  matrix containing all 
                                      pwRSMD relations between the 50k conformers.
     # dropallnan (bool): Default False. If running a molecule with DFT geom opt info, you'll want to drop
        # all molecules with any nan, as opposeed to only dropping rows with nan. I might replace this ***
    '''
    
    print('length of df ', len(df))
    print('length of times ', len(times))
    print('length of pwRMSD ', len(pwRMSD))

    # Initialize dataframe for capturing combo parameters, CCS results, and time
    df_results = pd.DataFrame(columns=['Time (min)', 'BW', 'BW std', 'LE', 'LE std', 'SA', 'SA std'])

    df_methods = pd.DataFrame(columns=['AMBER cycles', 'num conformers 1',
                                       'AMBER energy BET', 'num conformers 2',
                                       'SDS or random', 'num conformers 3',
                                       'DFT energy BET', 'num conformers 4',
                                       'SDS or random 2', 'num conformers 5',
                                       'Final method', 'final num conformers'])


    
    # Initialize parameters to sweep over. 
      #They will be modified after the initial downselection (number "generated" in AMBER)

    # Below energy threshold, units are kcal/mol
    BET_params_init = [1, 2, 5, 10, 20]
    #BET_params_init = [0.5, 1, 2, 5, 10, 15, 20, 25, 35]
   
    # Similarity down selection
    # Really, these are the "dissimilar" params. We will also choose the most similar every time, so we start with 0.
    SDS_params_init = [0, 3, 5, 10, 50, 100, 1000, 25000]
    #SDS_params_init = [0, 1, 3, 5, 10, 20, 30, 40, 50, 100, 200, 300, 500,
    #                   1000, 2000, 5000, 10000, 20000, 40000]
    
    # Random
    R_params_init = [1, 3, 5, 10, 50, 100, 1000, 25000]    
    #R_params_init = [1, 3, 5, 10, 20, 30, 40, 50, 100, 200, 300, 500,
    #                 1000, 2000, 5000, 10000, 20000, 40000]
    
    # if len(df.dropna()) != 0:
    #     df = df.dropna()
    #     DGO = True # I could use a bool like this for later, 
    #                 #as opposed to dropping nan again when I append DGO results ***         
    #     print('DGO flagged as True')
            
    # Start by simulating the initial number of AMBER cycles
    # When I do cycles, use the maximum number of cycles I get to for molid25 as the halfway point for comparison ***
    # Use the length of df if it has rows that have DGO elements (len of df, as modified by the if statement above)    
    for c in [1, 5, 10, 100, 500, 1000]:
        subdf = df.loc[:c*50-1].copy()
        print('cycles: ', c)
        print('length of subdf: ', len(subdf))
        print('subdf: \n', subdf)

        # Make the energy columns relative so the minimum energy is 0.
        #subdf['DFT Energy'] -= subdf['DFT Energy'].min() # dont do this here because Ill do in Tier2
        subdf['MD Energy (EPtot)'] -= subdf['MD Energy (EPtot)'].min()
        subdf['MD Energy (Etot)'] -= subdf['MD Energy (Etot)'].min()
        subdf['DFT GO Energy'] -= subdf['DFT GO Energy'].min()

        #print('MD Energy EPtot column \n', subdf['MD Energy (EPtot)'])
        #print('DFT GO column \n', subdf['DFT GO Energy'])

        # Initialize subdf as the full df without nan's
          # This is important so we don't "over capture" during SDS, or select non-existent conformers
          # If doing combo methods with a molecule that has DFT geometry optimization info, remove molecules
          # that do not have this information from consideration.    
        subdf = subdf.dropna(how='all')
        N = len(subdf)
        print('len of subdf after dropping nan', N)

        # Modify parameters to sweep over, based off the max energy or number of conformers still availiable

        # BET (AMBER energy, MD)    
        maxE = subdf['MD Energy (EPtot)'].max() #max energy, to capture the full energy range, and subsequently, all conformers    
        BET_params = [BET_params_init[x] for x in range(len(BET_params_init)) if BET_params_init[x] < maxE]
        BET_params.append(maxE) #add maximum, to capture all conformers
        # I could do a simple check here to make sure I am not appending maxE (essentially) twice (e.g. [..., 34.98, 35] ) ***
        # Could fix by making the ` < maxE` some range like  ` < (maxE - 5)` Except would be difficult if maxE <= 5... ***
        #print('max amber energy', maxE)
        #print('BET_params ', BET_params)

        # Initialize time
        # Let N be the number of conformers, or the cycles?.... ***
        totaltime = times['AMBER ave time'].values*N
        #print('total time after AMBER: ', totaltime)


        #### That concludes the initial MC set up. Here is were I begin the BET

        # Begin MD (AMBER) BET param sweep. I'm going to use Potential Energy, but I could change it or make it an input
         # variable.
        for bp in BET_params:
            subdf_b = subdf.loc[subdf['MD Energy (EPtot)'] <= bp].copy()
            sub_b_len = len(subdf_b)
            print('bp', bp)
            print('sub_b_len', sub_b_len)


            #### SDS     

            # Adjust the size of SDS params according to the size of the sub-dataframe after BET
            SDS_params = [SDS_params_init[x] \
                          for x in range(len(SDS_params_init)) if SDS_params_init[x] < (sub_b_len - 1)]
            SDS_params.append(sub_b_len - 1) # minus 1 to account for the most similar which will be appended as well         
#            print('SDS_params after adjustment', SDS_params)

            # Construct a new pwRSMD reduced to only the conformers which have been selected
            #pwrmsd = pwRMSD.loc[idx_b].copy().loc[idx_b].copy()
            idx_b = subdf_b.index
            pwrmsd = np.array([pwRMSD[i][idx_b] for i in idx_b]) #Check thisis faster than using dataframe 
                                                                    #method on cascade ***
#            print('len idx_b ', len(idx_b))
#            print('len pwrmsd ', len(pwrmsd))

            # Rank by pairwise RMSD dissimilarity
#            print('checking sub_b_len again: ',sub_b_len)
            idx_rank = nDis(pwrmsd, sub_b_len)
            rank = [x for x in range(1,sub_b_len+1)]
            rankdf = pd.DataFrame([idx_rank, rank], index=['index', 'rank']).T
#            print('rankdf before sorting \n', rankdf)
            rankdf = rankdf.sort_values(by='index')
            subdf_b['Dissimilarity Rank'] = rankdf['rank'].values #test this whole process I did ***
            
#            print('check the rank sorting worked ok')
#            print('Also check any changes I made in SDS tier 1 match SDS in tier 2')
#            print('idx_rank ', idx_rank)
#            print('rank ', rank)
#            print('rankdf after sorting \n', rankdf)
#            print('subdf_b[dissimilarity rank]', subdf_b['Dissimilarity Rank'])

            # Begin the sweep
            for sp in SDS_params:
                subdf_interm = subdf_b.loc[subdf_b['Dissimilarity Rank'] <= sp].copy()
#                print('subdf_interm ', subdf_interm)
#                print('subdf_b[Dissimilarity Rank].idxmax() ', subdf_b['Dissimilarity Rank'].idxmax())
                subdf_s = subdf_interm.append(subdf_b.loc[subdf_b['Dissimilarity Rank'].idxmax()].copy())
                sub_s_len = len(subdf_s)            
    #             subdf_interm = subdf_b.iloc[:sp]
    #             subdf_s = subdf_interm.append(subdf_b.iloc[-1])


                # If SDS is grabbing all, that's the same as not doing SDS
                if sp == (sub_b_len - 1):
                    totaltime_s = totaltime
                else:
                    number_of_pairs = sub_b_len * (sub_b_len-1) / 2 #sub_s_len * (sub_s_len-1) / 2
                    totaltime_s = totaltime + times['SDS time per pair'].values * number_of_pairs
#                    print('adjust this')
                    # I will need to adjust this to estimate how long it would take to 
                    # create the pwRMSD at this population size ***

                df_s, df_sm = Tier2(subdf_s.copy(), 1, sub_s_len, [c, N, bp, sub_b_len, 'SDS', sp+1]) #, totaltime_s)

                # Modify with current time and combo *** reeword wht I just said heere, can't brain
                df_s['Time (min)'] += totaltime_s
                #df_s['Methods'] = [c, N, bp, sub_b_len, 'SDS', sp+1] + df_s['Methods']
                df_results = df_results.append(df_s)
                df_methods = df_methods.append(df_sm)


            #### Random

            # Begin random selection param sweep
            # Adjust the size of the param sweep according to how many conformers are left after BET
            R_params = [R_params_init[x] for x in range(len(R_params_init)) if R_params_init[x] < int(sub_b_len/2)]
            R_params.append(int(sub_b_len/2)) # I've changed this to only grab a half, since I'm letting the
                                            # full range be grabbed through the SDS branch. 

            # Make sure to grab at least one conformer. This will lead to redundant information,
            # but it will at least allow the algorithm to continue. Redudancies can be filtered out later.
            # OR I could tell it to quit right here. Not sure what I would have it return... ***
            if R_params == [0]:
                R_params = [1]

#            print('since I am grabbing up to half the full range, \
#                  make sure R params make sense to compare against SDS still')

            # Launch Monte Carlo simulation
            #for rp2 in R2_params:
            #    if rp2 == sub_b_len:
            #        # Only bother with one iteration if you're grabbing the full sub-population everytime anyway
            #        MC2_iter = 1
            #    else:
            #        MC2_iter = 10 # could make this a variable
            MC_iter = 10

            for rp in R_params:

                df_r, df_rm = Tier2(subdf_b.copy(), MC_iter, rp, [c, N, bp, sub_b_len, 'random', rp])
                df_r['Time (min)'] += totaltime # might move these by turning into params ***
                # might move these by turning into params ***
                #df_r['Methods'] = [c, N, bp, sub_b_len, 'random', rp] + df_s['Methods']
                df_results = df_results.append(df_r)
                df_methods = df_methods.append(df_rm)




            
    return df_results, df_methods


####################################### 4


def Tier2(df, MC_iter, num_confs, precombo): #, totaltime): #=0, combo_so_far=[], combos=[]): 
    # Make energy relative to minimum in the sub-sample
    df['DFT Energy'] -= df['DFT Energy'].min()
    #print('df[\'DFT Energy\'] \n', df['DFT Energy'])
    #print('halfmaxE', df['DFT Energy'].max()/2)
    # Initialize parameters
    # Below energy threshold, units are kcal/mol
    BET2_params_init = [0.5, 1, 2, 5, 15]   #10, 15, 20, 25, 35
    
    # Modify BET2 params according to max energy available
    halfmaxE = df['DFT Energy'].max()/2 #half max energy, to capture the half energy range as opposed to whole,
                                        #so we're not doing SDS on SDS, random on random, random on SDS, etc.
                                        #because the first tier already accounts for grabbing everything
    BET2_params = [BET2_params_init[x] for x in range(len(BET2_params_init)) if BET2_params_init[x] < halfmaxE]
    BET2_params.append(halfmaxE)
    #print('BET2_params', BET2_params)
    
    # Similarity down selection, "dissimilar" parameters
    SDS2_params_init = [0, 3, 5, 10, 50, 100, 1000, 25000]
    
    # Random
    R2_params_init = [1, 3, 5, 10, 50, 100, 1000, 25000]
#    print('make sure these parameters are consistent, or else make sense with, the params in tier 1')
#    print('For BET, SDS, and random')
    
    
    # The results have to be averaged over the iterations, so let's dump everything into a matrix
    # Each row will contain the results of every iteration for each of the x final method paths
    # We will then take the average and stdev of each row.
    # The combination of methods only needs to be declared once for each path, as does time
    
    # Just set up the maximum number possible, which will leave a bunch of blanks at the end
    # and then drop the rows of nan's it creates when turned into a pandas dataframe at the end of the simulation
    # The +2 for the two only done once (not under BET), and the other +2 incase 
    # the params after adjustment increase one each (such as the case when the length of the BET
    # sub-dataframe is greater than all of the parameters in the <method>_params_init    
    num_cells = int(3 + 3*len(BET2_params)*(len(SDS2_params_init) + len(R2_params_init) + 2)) #I changed 2+3* to 1+2*  ***

    timeM = [[] for x in range(num_cells)]
    comboM = [[] for x in range(num_cells)]
    bwM = [[] for x in range(num_cells)]
    leM = [[] for x in range(num_cells)]
    saM = [[] for x in range(num_cells)]
    
    totaltimes = []
    combos = []
    bw = []
    bw_stdev = []
    le = []
    le_stdev = []
    sa = []
    sa_stdev = []

    # Count up time under tier 2. Tier 1 will add tier 2's time to its time, and append combos to its combos.
    totaltime = 0 # ***, or I could not. I'll know at the end what the better method is
    #combos = []

    
    for i in range(MC_iter):
        idx = np.random.choice(df.index, num_confs, replace=False)
        subdf = df.loc[idx].copy()
        sub_len = len(subdf)


        # First, do without a DFT energy threshold 
        
        # AMBER Energies
        timeM[0] = totaltime + times['MOBCAL ave time'].values * sub_len
        comboM[0] = precombo + [np.nan, np.nan, 'nan', np.nan, 'AMBER', sub_len]
        bwM[0].append(BW(subdf['MD CCS'].values, subdf['MD Energy (EPtot)'].values))
        leM[0].append(LE(subdf['MD CCS'].values, subdf['MD Energy (EPtot)'].values))
        saM[0].append(SA(subdf['MD CCS'].values))                    

        # DFT Energies
        timeM[1] = totaltime + (times['DFT Energy ave time'].values + times['MOBCAL ave time'].values) * sub_len
        comboM[1] = precombo + [np.nan, np.nan, 'nan', np.nan, 'DFT E', sub_len]
        bwM[1].append(BW(subdf['MD CCS'].values, subdf['DFT Energy'].values))
        leM[1].append(LE(subdf['MD CCS'].values, subdf['DFT Energy'].values))
        saM[1].append(SA(subdf['MD CCS'].values)) # Could make this np.nan since SA is independent of energy

        # DFT GO
        timeM[2] = totaltime  \
                   + (times['DFT Geom Opt ave time'].values  \
                   + times['MOBCAL ave time'].values)  \
                   * sub_len

        comboM[2] = precombo + [np.nan, np.nan, 'nan', np.nan, 'DFT GO', sub_len]
        bwM[2].append(BW(subdf['DFT GO CCS'].values, subdf['DFT GO Energy'].values))
        leM[2].append(LE(subdf['DFT GO CCS'].values, subdf['DFT GO Energy'].values))
        saM[2].append(SA(subdf['DFT GO CCS'].values))


        #### Begin second tier BET (DFT energy)

        # Set DFT energies relative to the minimum
        subdf['DFT Energy'] -= subdf['DFT Energy'].min()        

        k = 3

        for bp2 in BET2_params:
            #print('bp2', bp2)
            #print('subdf \n', subdf)
            subdf_b2 = subdf.loc[subdf['DFT Energy'] <= bp2].copy()
            #print('subdf_b2 \n', subdf_b2)
            sub_b2_len = len(subdf_b2)
            #print('sub_b2_len', sub_b2_len)
            totaltime_b = totaltime + times['DFT Energy ave time'].values * sub_len
            
            # Modify SDS and random parameters
#            print('I do not need grab the full range for both SDS and random, so I should only do so for one')
#            print('This ^ also applies for tier 1, so I should modify that as well.')
#            print('I should leave out the full range for random, since I do not want to do MC on full range anyway.')
            SDS2_params = [SDS2_params_init[x] for x in range(len(SDS2_params_init)) \
                           if SDS2_params_init[x] < (sub_b2_len - 1)]
            SDS2_params.append(sub_b2_len - 1) # minus 1 to leave room for the most similar
            #print('SDS2_params', SDS2_params)
            
            
            
            #### SDS, 2nd tier
            

            # Construct a new pwRSMD reduced to only the conformers which have been selected
            idx_b = subdf_b2.index
            pwrmsd = np.array([pwRMSD[x][idx_b] for x in idx_b]) #Check thisis faster 
                                                                  #than using dataframe method on cascade ***
            print('idx_b ', idx_b)
            # Rank by pairwise RMSD dissimilarity
            idx_rank = nDis(pwrmsd, sub_b2_len)
            rank = [x for x in range(1,sub_b2_len+1)]
            rankdf = pd.DataFrame([idx_rank, rank], index=['index', 'rank']).T
            rankdf = rankdf.sort_values(by='index')
            subdf_b2['Dissimilarity Rank'] = rankdf['rank'].values

            # Begin SDS sweep
            for sp2 in SDS2_params:
                subdf_interm = subdf_b2.loc[subdf_b2['Dissimilarity Rank'] <= sp2].copy()
                subdf_s2 = subdf_interm.append(subdf_b2.loc[subdf_b2['Dissimilarity Rank'].idxmax()].copy())
                sub_s2_len = len(subdf_s2)            

                # If SDS is grabbing all, that's the same as not doing SDS
                if sp2 == (sub_b2_len - 1):
                    totaltime_s = totaltime_b
                else:
                    number_of_pairs = sub_b2_len * (sub_b2_len-1) / 2
                    totaltime_s = totaltime_b + times['SDS time per pair'].values * number_of_pairs            
            

                # Collect times, method combinations, and BW, LE, SA CCS's

                # DFT Energies
                timeM[k] = totaltime_s  \
                           + (times['DFT Energy ave time'].values  \
                           + times['MOBCAL ave time'].values)  \
                           * sub_s2_len

                comboM[k] = precombo + [bp2, sub_b2_len, 'SDS', sp2+1, 'DFT E', sub_s2_len] #[bp, sp, rp, bp2, 'DFT Energy', sub_b2_len]
                bwM[k].append(BW(subdf['MD CCS'].values, subdf['DFT Energy'].values))
                leM[k].append(LE(subdf['MD CCS'].values, subdf['DFT Energy'].values))
                saM[k].append(SA(subdf['MD CCS'].values))

                # AMBER Energies
                timeM[k+1] = totaltime_s  \
                           + (times['DFT Energy ave time'].values  \
                           + times['MOBCAL ave time'].values)  \
                           * sub_s2_len

                comboM[k+1] = precombo + [bp2, sub_b2_len, 'SDS', sp2+1, 'AMBER', sub_s2_len]
                bwM[k+1].append(BW(subdf['MD CCS'].values, subdf['MD Energy (EPtot)'].values))
                leM[k+1].append(LE(subdf['MD CCS'].values, subdf['MD Energy (EPtot)'].values))
                saM[k+1].append(SA(subdf['MD CCS'].values))

                # DFT GO
                timeM[k+2] = totaltime_s  \
                             + (times['DFT Energy ave time'].values  \
                             + times['DFT Geom Opt ave time'].values  \
                             + times['MOBCAL ave time'].values)  \
                             * sub_s2_len

                comboM[k+2] = precombo + [bp2, sub_b2_len, 'SDS', sp2+1, 'DFT GO', sub_s2_len]
                bwM[k+2].append(BW(subdf['DFT GO CCS'].values, subdf['DFT GO Energy'].values))
                leM[k+2].append(LE(subdf['DFT GO CCS'].values, subdf['DFT GO Energy'].values))
                saM[k+2].append(SA(subdf['DFT GO CCS'].values))

                k += 3

#            print('k: ', k)
#            print('the otehr side:', 2 + 3*len(BET2_params)*len(SDS2_params))
            #assert k == 2 + 3*len(BET2_params)*len(SDS2_params) # I think this is bad. Idk why it appeared to be passin before
            
            
            
            #### Random, 2nd tier

            R2_params = [R2_params_init[x] for x in range(len(R2_params_init)) if R2_params_init[x] < int(sub_b2_len/2)]
            R2_params.append(int(sub_b2_len/2)) # Again, make sure doing 1/2 b2 len makes sense for these parameters ***             
#            print('R2_params: ', R2_params)

            # Make sure to grab at least one conformer. This will lead to redundant information,
            # but it will at least allow the algorithm to continue. Redudancies can be filtered out later.
            # OR I could tell it to quit right here. Not sure what I would have it return... ***
            if R2_params == [0]:
                R2_params = [1]

            # Initalize a matrix for the second MC
            num_cells2 = int(3 * len(R2_params)) # here I dont need to add an extra +1 because the params have already been modified.
            timeM2 = [[] for x in range(num_cells2)]
            comboM2 = [[] for x in range(num_cells2)]
            bwM2 = [[] for x in range(num_cells2)]
            leM2 = [[] for x in range(num_cells2)]
            saM2 = [[] for x in range(num_cells2)]
            
            MC_iter2 = 10 # Could change this ***
            
            k2 = 0            
            
            for rp2 in R2_params:
#                print('rp2', rp2)
                for j in range(MC_iter2):
                    idx = np.random.choice(subdf_b2.index, rp2, replace=False)
                    subdf_r2 = subdf_b2.loc[idx].copy()
                    sub_r2_len = len(subdf_r2)

                    # DFT Energies
                    timeM2[k2] = totaltime_b  \
                               + (times['DFT Energy ave time'].values  \
                               + times['MOBCAL ave time'].values)  \
                               * sub_r2_len

                    comboM2[k2] = precombo + [bp2, sub_b2_len, 'random', rp2, 'DFT E', sub_r2_len] 
                    bwM2[k2].append(BW(subdf['MD CCS'].values, subdf['DFT Energy'].values))
                    leM2[k2].append(LE(subdf['MD CCS'].values, subdf['DFT Energy'].values))
                    saM2[k2].append(SA(subdf['MD CCS'].values))

                    
                    # AMBER Energies
                    timeM2[k2+1] = totaltime_b  \
                               + (times['DFT Energy ave time'].values  \
                               + times['MOBCAL ave time'].values)  \
                               * sub_r2_len

                    comboM2[k2+1] = precombo + [bp2, sub_b2_len, 'random', rp2, 'AMBER', sub_r2_len]
                    bwM2[k2+1].append(BW(subdf['MD CCS'].values, subdf['MD Energy (EPtot)'].values))
                    leM2[k2+1].append(LE(subdf['MD CCS'].values, subdf['MD Energy (EPtot)'].values))
                    saM2[k2+1].append(SA(subdf['MD CCS'].values))


                    # DFT GO
                    timeM2[k2+2] = totaltime_b  \
                                   + (times['DFT Energy ave time'].values  \
                                   + times['DFT Geom Opt ave time'].values  \
                                   + times['MOBCAL ave time'].values)  \
                                   * sub_r2_len

                    comboM2[k2+2] = precombo + [bp2, sub_b2_len, 'random', rp2, 'DFT GO', sub_r2_len]
                    bwM2[k2+2].append(BW(subdf['DFT GO CCS'].values, subdf['DFT GO Energy'].values))
                    leM2[k2+2].append(LE(subdf['DFT GO CCS'].values, subdf['DFT GO Energy'].values))
                    saM2[k2+2].append(SA(subdf['DFT GO CCS'].values))
                
                k2 += 3 # On to next rp
            
#            print('num_cells2', num_cells2)
#            print('k2', k2)

            for cell in range(num_cells2): # check this is right, should it be k2-2? ***
#                print('k ', k)
#                print('cell ', cell)
#                print('timeM2', timeM2)
#                print('timeM2[cell] ', timeM2[cell])
                #timeM[k+cell].append(timeM2[cell])
                #comboM[k+cell].append(comboM2[cell])
                timeM[k+cell] = timeM2[cell]  
                comboM[k+cell] = comboM2[cell]

                #take an average (and sstdev?)
                bwM[k+cell].append(np.nanmean(bwM2[cell])) # include the stdev, modify with an appropriate equation from google search ***
                leM[k+cell].append(np.nanmean(leM2[cell]))
                saM[k+cell].append(np.nanmean(saM2[cell]))

            k += cell + 1
            

    for cell in range(k): # I thinnk this could be chnage to range(k) not range(num_cells)
        totaltimes.append(timeM[cell])                  
        combos.append(comboM[cell])

        bw.append(np.nanmean(bwM[cell]))
        bw_stdev.append(np.nanstd(bwM[cell]))

        le.append(np.nanmean(leM[cell]))
        le_stdev.append(np.nanstd(leM[cell]))  

        sa.append(np.nanmean(saM[cell]))
        sa_stdev.append(np.nanstd(saM[cell]))

    
    
    df_results = pd.DataFrame([totaltimes, bw, bw_stdev, le, le_stdev, sa, sa_stdev],
                              index=['Time (min)', 'BW', 'BW std', 'LE', 'LE std', 'SA', 'SA std']).T
    df_results = df_results.loc[df_results['Time (min)'].str.len() != 0]
    df_results['Time (min)'] = df_results['Time (min)'].astype(float)

    df_methods = pd.DataFrame(combos, 
                          columns=['AMBER cycles', 'num conformers 1', 
                                   'AMBER energy BET', 'num conformers 2',
                                   'SDS or random', 'num conformers 3',
                                   'DFT energy BET', 'num conformers 4',
                                   'SDS or random 2', 'num conformers 5',
                                   'Final method', 'final num conformers'])
    df_methods = df_methods.loc[df_results.index]

    return df_results, df_methods


####################################### 5


if __name__ == '__main__':
    from time import time
    start = time()    
    
    parser = argparse.ArgumentParser()
    parser.add_argument('molid', type=str, help='Molecular identifier, i.e. the \'02\' in \'molid02\'')
    args = parser.parse_args()    
    ID = args.molid
                        
    directory = abspath('.') # Assume the files are in the working directory
    
    # Load dataframe
    df = pd.read_csv(join(directory, f'molid{ID}_combo_methods.csv'))
    
    # Load time information
    times = pd.read_csv(join(directory, f'molid{ID}_combo_times.csv'))
    
    # Load pwRMSD dataframe
    mtrx_file = glob.glob('*_50k_50k_reflect_logsum.pkl')[0] # could change this and the others ^ to arguments in argparse
    pwRMSD = pd.read_pickle(mtrx_file) # this will necessitate I work on the supercomputer. Need 20+ GB RAM
    pwRMSD = pwRMSD.values

    
    # Run combo methods paramter sweep simulation
    results, methods = combo_methods(df, times, pwRMSD)
    results.to_csv(join(directory, f'molid{ID}_combo_methods_DGO_results-v2.csv'), index=False)
    methods.to_csv(join(directory, f'molid{ID}_combo_methods_DGO_resultsM-v2.csv'), index=False)
    
    print((time()-start)/3600, ' hours')

