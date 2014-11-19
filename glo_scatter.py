import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import statsmodels.api as sm
from scipy.stats import linregress

colnames = ("smdc","dock","glo","3","4","5","6")
sar = pd.read_csv('359series_oct2014.csv', names = colnames, header=0)
sar.drop(["3","4","5","6"],axis=1, inplace='true')

fig=plt.figure()
ax=fig.add_subplot(3,1,1)
ax.plot(sar.glo,sar.dock,'bo')

ax.set_title('GLIDEXP Docking Score vs IC50 for 62 Phenyl Indole Cmpds')
ax.set_xlabel('IC50 (uM)')
ax.set_ylabel('GLIDE XP score')

ax2 = fig.add_subplot(3,1,2)
ax2.plot(sar.glo[sar.glo < 4],sar.dock[sar.glo < 4],'bo')
ax2.set_title('GLIDEXP Docking Score vs IC50 **IC50 < 3 uM only**')
ax2.set_xlabel('IC50 (uM)')
ax2.set_ylabel('GLIDE XP score')

plt.subplots_adjust(wspace=0,hspace=0.5)

x = sar.glo[sar.glo < 4]
y = sar.dock[sar.glo < 4]

result = sm.OLS(y,x).fit()
rsquared = result.rsquared

slope, intercept, r_value, p_value, std_err = linregress(x, y)
fit_fn = np.poly1d([.999, -8.5])

ax2.plot(x, fit_fn(x), 'k--')
ax2.annotate('r^2 = 0.50', xy=(2.1, -6.5), xytext=(2.5,-9.3),
            arrowprops=dict(facecolor='black', shrink=0.05),
            )


ax3 = fig.add_subplot(3,1,3)
ax3.plot(result.resid, "k")
ax3.set_title('Linear Least Squares Residuals')
ax3.set_ylabel('Residual')



fig.show()

