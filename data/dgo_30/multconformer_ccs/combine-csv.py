import glob
import pandas as pd

#inchikeys = ['BHQCQFFYRZLCQQ-OELDTZBJSA-N_+Na','CZMRCDWAGMRECN-UGDNZRGBSA-N_-H','DDRJAANPRJIHGJ-UHFFFAOYSA-N_+Na',
#             'HSMNQINEKMPTIC-UHFFFAOYSA-N_-H','NNICRUQPODTGRU-SVGMAFHSNA-N_+H','QWIZNVHXZXRPDR-WSCXOGSTSA-N_+H',
#             'RWZYAGGXGHYGMB-UHFFFAOYSA-N_-H','UTAIYTHAJQNQDW-KQYNXXCUSA-N_+H','UYARPHAXAJAZLU-KQYNXXCUSA-N_-H',
#             'VQAYFKKCNSOZKM-IOSLPCCCSA-N_+Na','ZROGCCBNZBKLEL-MPRHSVQHSA-N_+Na','ZYEMGPIYFIJGTP-UHFFFAOYSA-N_+Na']
#molid = ['16','25','11','22','10','06','21','03','24','14','17','12']
inchikeys = ['AUNGANRZJHBGPY-SCRDCRAPSA-N_+H', 'BXNJHAXVSOCGBA-UHFFFAOYSA-N_+H', 'DFPMSGMNTNDNHN-ZPHOTFPESA-N_-H',
             'QBUVFDKTZJNUPP-BBROENKCNA-N_+Na', 'UVLWLKCNNYTXDT-XDTORHTBNA-N_+Na', 'WWUZIQQURGPMPG-CCEZHUSRSA-N_+H']
molid = ['05','02','28',
         '19','18','04']


for ID in molid:
    ccs_file = f'~/Documents/excel-and-data/conformer_ccs/dgo_30/multconformer_ccs/molid{ID}_30dgo_ccs.csv'
    energies_file = f'~/Documents/excel-and-data/conformer_ccs/dgo_30/multconformer_ccs/molid{ID}_30dgo_dft_energy_map.csv'

    ccs_df = pd.read_csv(ccs_file)
    energies_df = pd.read_csv(energies_file)
    ccs_df['cycles'] = energies_df['cycle']
    ccs_df['geometry'] = energies_df['geometry']
    ccs_df['first_dft_energy'] = energies_df['first_dft_energy']
    ccs_df = ccs_df.rename(columns={'cycles': 'cycle'})
    ccs_df.to_csv(f'~/Documents/excel-and-data/conformer_ccs/dgo_30/multconformer_ccs/molid{ID}_ccs_30dgo.csv', index=False)
