IT=0
from scipy.optimize import minimize
def getError(temp):
    global IT
    x=temp[0]
    y=temp[1]
    if (0==0): print(IT)
    IT=IT+1
    return(x**2+4*x+4+y**2)

res = minimize(getError,[2,-1],method='SLSQP', bounds=[[-3,3],[0,2]])
print(res)
