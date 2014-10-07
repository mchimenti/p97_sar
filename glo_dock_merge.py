from pandas import read_table, merge
from numpy import mean, std, polyfit
import matplotlib.pyplot as pyp

cmpd_out = read_table('cmpd_out.csv', sep=',',header=0)
cmpd_dock = read_table('359_dock_oct2014.csv', sep=',',header=0)

cmpd_dock.columns= ['SMDC','GLIDE_score', 'GLO average', '3','4','5','6']
cmpd_dock.drop(['GLO average','3', '4', '5', '6'], axis=1, inplace=True)

cmpd_dock_GLO = merge(cmpd_out, cmpd_dock, on='SMDC')

##plot

fit = polyfit(cmpd_dock_GLO.GLO_average, cmpd_dock_GLO.GLIDE_score, 1)

pyp.scatter(cmpd_dock_GLO.GLO_average, cmpd_dock_GLO.GLIDE_score)
pyp.show()