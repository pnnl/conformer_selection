import pandas as pd
import glob
import numpy as np

inchidf = pd.read_csv('inchi-ID-table.csv')

paths = glob.glob('*.tsv')
paths.sort()


for mol in paths:
    df = pd.read_csv(mol, delim_whitespace=True)

    length = len(df)
    placeholder = [np.nan for x in range(length)]
    csvdf = pd.DataFrame()
    csvdf['cycles'] = placeholder
    csvdf['geometry'] = placeholder
    csvdf['ccs'] = df['ccs']
    csvdf['dft_energy'] = df['dft_energy']

    dot = mol.split('.')[0]
    molid = dot.split('_')[0]
    row = inchidf.loc[inchidf['Inchi Key'] == molid]
    ID = row['ID'].values[0]

    csvdf.to_csv(f'{ID}_30dgo_ccs.csv', index=False)
