###script to read an ADP GLO excel report and a docking report
###clean and merge the tables together on smdc number

###input requires two csv files, one with ADP GLO data in formatted columns
###with order "structure (empty)", "SMILES", "SMDC", "GLO(IC50)", "LOT" and "NOTES"
###the other file contains docking data with columns "SMDC", "GLIDEscore"

###NOTE: be sure to strip the MS-DOS style "^M" characters from the Excel created csv files
###Use the vi command: "%s/control-V, ^M/control-V, ^M/g"  this works for some reason

import numpy as np
import matplotlib.pyplot as plt
import statsmodels.api as sm
import pandas as pd
from scipy.stats import linregress

cmpd = pd.read_table('raw_sar_input.csv', sep=',',header=0)
cmpd.columns = ['structure', 'SMILES', 'SMDC', 'GLO', 'LOT', 'NOTES']

cmpd.GLO = cmpd.GLO.replace('> 50.0', 50)
cmpd.GLO = cmpd.GLO.replace('> 50', 50)
cmpd.GLO = cmpd.GLO.astype('float64', copy=False)
cmpd.SMDC = cmpd.SMDC.astype('int32', copy=False)

##calculate mean IC50s

cmpd_mean = cmpd.groupby('SMDC').agg([np.mean, np.std])
cmpd_mean['SMDC'] = cmpd_mean.index  ## create new SMDC column from SMDC index
cmpd_mean = cmpd_mean.reset_index(drop=True)  

cmpd_smiles = cmpd[['SMILES','SMDC']]

##merge mean IC50s with raw data

cmpd_final = pd.merge(cmpd_mean, cmpd_smiles, on='SMDC')
cmpd_final = cmpd_final[cmpd_final["SMDC"] != cmpd_final["SMDC"].shift()]  ##delete identical rows

##clean up final sar dataframe

cmpd_final.drop((u'SMDC', u''), axis=1, inplace=True)
cmpd_final.reset_index(drop=True, inplace=True)
cmpd_final.columns = ['SMDC','GLO_ave', 'GLO_std', 'SMILES']

##read in docking data and merge with sar dataframe

cmpd_dock = pd.read_table('docking_scores.csv', sep=',',header=0)
cmpd_dock.columns= ['SMDC','g_score']
cmpd_dock_GLO = pd.merge(cmpd_final, cmpd_dock, on='SMDC')

##plot

fig=plt.figure()
ax=fig.add_subplot(3,1,1)
ax.plot(np.log(cmpd_dock_GLO.GLO_ave),cmpd_dock_GLO.g_score,'bo')

ax.set_title('GLIDEXP Docking Score vs IC50 for Phenyl Indole Cmpds')
ax.set_xlabel('Log[IC50] (uM)')
ax.set_ylabel('GLIDE XP score')

ax2 = fig.add_subplot(3,1,2)
ax2.plot(np.log(cmpd_dock_GLO.GLO_ave[cmpd_dock_GLO.GLO_ave < 4]),cmpd_dock_GLO.g_score[cmpd_dock_GLO.GLO_ave < 4],'bo')
ax2.set_title('GLIDEXP Docking Score vs IC50 **IC50 < 4 uM only**')
ax2.set_xlabel('Log [IC50] (uM)')
ax2.set_ylabel('GLIDE XP score')

plt.subplots_adjust(wspace=0,hspace=0.5)

x = cmpd_dock_GLO.GLO_ave[cmpd_dock_GLO.GLO_ave < 4]
y = cmpd_dock_GLO.g_score[cmpd_dock_GLO.GLO_ave < 4]

result = sm.OLS(y,x).fit()
rsquared = result.rsquared

slope, intercept, r_value, p_value, std_err = linregress(x, y)
fit_fn = np.poly1d([slope, intercept])

ax2.plot(x, fit_fn(x), 'k--')
ax2.annotate('r^2 =' + str(rsquared), xy=(1, fit_fn(1)), xytext=(2.5,-9.3),
            arrowprops=dict(facecolor='black', shrink=0.05),
            )


ax3 = fig.add_subplot(3,1,3)
ax3.plot(result.resid, "k")
ax3.set_title('Linear Least Squares Residuals')
ax3.set_ylabel('Residual')



fig.show()
