import glob
import pandas as pd

paths = glob.glob('*-N_+*_*_geom*.tsv')
paths.sort()

df = pd.read_csv(paths[0], delimiter='\t')
for p in paths[1:]:
    df2 = pd.read_csv(p, delimiter='\t')
    df.append(df2)

df.to_csv('combined.tsv', index=False)
