##http://datagenetics.com/blog/november22017/index.html
import numpy as np
E=np.arange(0,101)

## function to recursively calculate 
## mean and standard deviation
def recursiveStats(mu_old,sigma_old,xn,n):
    mu_new     = ( mu_old*(n-1) + xn )/n
    sigma_new  = np.sqrt( ( (n-1)*sigma_old**2 + (xn-mu_old)*(xn-mu_new) )/n )
    return(mu_new,sigma_new)

## calculating mean and standard deviation
## by calling a funciton
mu=0
sigma=0
for it,val in enumerate(E):
    mu,sigma=recursiveStats(mu,sigma,val,it+1)
    if ((it+1)%10==1): print(mu,sigma)
    
## calculating mean and standard deviation
## in line
tot=0   #sum of values
tot2=0  #sum of squared vales
for it,val in enumerate(E):
    tot    = tot  + val                       #sum of values
    tot2   = tot2 + val**2                    #sum of squared vales
    mean   = tot/(it+1)                       #mean of values ; since "it" starts from zero
    #mean2  = mean**2                         #squared mean
    #mtot2  = tot2/(it+1)                     #mean of squared sum
    #sigma2 = mtot2-mean2                     #sigma squared = mean of squared sum - squared mean
    sigma  = np.sqrt( tot2/(it+1) - mean**2 ) #since "it" starts from zero
    if ((it+1)%10==1): print((it+1),mean,sigma)
