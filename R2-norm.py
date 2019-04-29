import matplotlib
import matplotlib.pyplot as plt

import numpy as np
from scipy import stats

x=np.arange(10)     #x values range from 0 to 9
y=np.zeros(10)      #y values range from 0 to 9 with small error
for i in range(10):
    y[i]=i+np.random.randn()*1
y_mean=np.mean(y)   #mean of y values

#fitting y=f(x)
slope, intercept, r_value, p_value, std_err = stats.linregress(x,y)
y_line = slope*x+intercept             #y_line=f(x)
plt.plot(x,y,linestyle='', marker='o') #plotting data
plt.plot(x,y_line)                     #plotting fitted line

#Calculating the R-value
SSRES  = 0.0 ##Residual   sum of squares
SSTOT  = 0.0 ##Total      sum of squares
SSREG  = 0.0 ##Regression sum of squares
for i in range(len(x)) :
    SSRES = SSRES + (y[i]      - y_line[i])**2
    SSTOT = SSTOT + (y[i]      - y_mean   )**2
    SSREG = SSREG + (y_line[i] - y_mean   )**2
r_comp=np.sqrt(1-SSRES/SSTOT)

print("r_value from algorithm     =",r_value,"\nr_value from 1-SSRES/SSTOT =",r_comp)
print("SSTOT-SSRES-SSREG = 0")
print("SSTOT-SSRES-SSREG =",SSTOT-SSRES-SSREG, "  (Should be approximately zero)")
