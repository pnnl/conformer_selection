'''Generates RDKit conformers. Saves conformers into 'cycles' to allow use of the same Monte Carlo simulation
scripts used on AMBER generated conformers without modification.  

Option to optimize using UFF uncomment below
'''

from rdkit import Chem
from rdkit.Chem import AllChem
from os.path import abspath, join, exists
from os import makedirs
from openbabel import pybel #obabel v3
#import pybel # obabel v2.4.1
from time import time

if __name__ == '__main__':
    START = time()

    smi = # e.g. '[OH2+][C@@H](C#N)C1=CC=CC=C1'
    id_add = # e.g. 'NNICRUQPODTGRU-SVGMAFHSNA-N_+H'
    m = Chem.MolFromSmiles(smi)
    m2 = Chem.AddHs(m)
    cnm = 1000
    gnm = 50

    # Make output and mol directory
    dirmol = 'output/mol'
    if not exists(dirmol):
        makedirs(dirmol)

    ids = list(AllChem.EmbedMultipleConfs(m2, numConfs=cnm*gnm)) #, randomSeed=r))

    for c in range(cnm+1): # Cycles

        # Fake the cycles and geometries. 
         # This will allow using previous conformer selection scripts without the need to modify them.
        start = c * gnm
        stop = (c + 1) * gnm

        for g, i in enumerate(ids[start:stop]):
            cycle = '%04d' % (c+1)
            geom = '%02d' % (g+1)

            # Optimize using UFF
            #AllChem.UFFOptimizeMolecule(m2, confId=i)

            # Write a temporary molfile
            molfile = join(dirmol, f'{id_add}_{cycle}_geom{geom}.mol')
            Chem.SDWriter(molfile).write(m2, confId=i)

            # Make specific dft directory for the .xyz, which corresponds to how the directories are set up
              # for AMBER SA conformers. 
            dirxyz = join('output', 'dft', f'{id_add}', f'cycle_{cycle}_geom{geom}')
            if not exists(dirxyz):
                makedirs(dirxyz)

            # Convert to .xyz
            mol = list(pybel.readfile('mol', molfile))
            xyzfile = join(dirxyz, f'{id_add}_{cycle}_geom{geom}.xyz')
            mol[0].write('xyz', xyzfile, True)

    print('Finished at: ', (time()-START)/60)









