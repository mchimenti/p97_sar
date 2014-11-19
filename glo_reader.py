from pandas import read_table, merge
from numpy import mean, std

cmpd = read_table('359_SAR.csv', sep=',',header=0)
cmpd.columns = ['structure', 'SMILES', 'SMDC', 'GLO', 'LOT', 'NOTES']

cmpd.GLO = cmpd.GLO.replace('> 50.0', 50)
cmpd.GLO = cmpd.GLO.astype('float64', copy=False)
cmpd.SMDC = cmpd.SMDC.astype('int32', copy=False)

cmpd_mean = cmpd.groupby('SMDC').agg([mean, std])
cmpd_mean['SMDC'] = cmpd_mean.index  ## create new SMDC column from SMDC index
cmpd_mean = cmpd_mean.reset_index(drop=True)  

cmpd_smiles = cmpd[['SMILES','SMDC']]

cmpd_final = merge(cmpd_mean, cmpd_smiles, on='SMDC')
cmpd_final = cmpd_final[cmpd_final["SMDC"] != cmpd_final["SMDC"].shift()]  ##delete identical rows

##clean up final dataframe

cmpd_final.drop((u'SMDC', u''), axis=1, inplace=True)
cmpd_final.reset_index(drop=True, inplace=True)
cmpd_final.columns = ['SMDC','GLO_average', 'GLO_std', 'SMILES']

##write to file

cmpd_final.to_csv('./cmpd_out.csv', index=False)