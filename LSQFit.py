import numpy as np
from numpy.linalg import inv
from matplotlib import pyplot as plt
from fit_fcn_to_data import fit_fcn_to_data

xmin=1.0
xmax=20.0
npoints=100
sigma=0.2
lx=np.zeros(npoints)
ly=np.zeros(npoints)
ley=np.zeros(npoints)
pars=[0.5,1.3,0.5]


from math import log
def f(x,par):
    return par[0]+par[1]*log(x)+par[2]*log(x)*log(x)

from random import gauss
def getX(x):  # x = array-like
    step=(xmax-xmin)/npoints
    for i in range(npoints):
        x[i]=xmin+i*step
        
def getY(x,y,ey):  # x,y,ey = array-like
    for i in range(npoints):
        y[i]=f(x[i],pars)+gauss(0,sigma)
        ey[i]=sigma


# the function corresponding to each term 
fit_terms = [lambda x: 1., lambda x: np.log(x), lambda x: np.log(x) ** 2]  

# *** modify and add your code here ***
nexperiments = 10000  # for exampxle

# perform many least squares fits on different pseudo experiments here
par_a = np.empty(nexperiments)
par_b = np.empty(nexperiments)
par_c = np.empty(nexperiments) 
chi2_reduced = np.empty(nexperiments) 

for i in range(0, nexperiments):

    getX(lx)
    getY(lx, ly, ley)
    coeffs, chi2 = fit_fcn_to_data(fit_terms, lx, ly, ley)

    par_a[i] = coeffs[0]
    par_b[i] = coeffs[1]
    par_c[i] = coeffs[2]

    chi2_reduced[i] = chi2


chi2_reduced *= 1./nexperiments

fig, axs = plt.subplots(2, 2)
plt.tight_layout()

# careful, the automated binning may not be optimal for displaying your results!
axs[0, 0].hist2d(par_a, par_b, bins=[100,100])
axs[0, 0].set_title('Parameter b vs a')

axs[0, 1].hist2d(par_a, par_c, bins=[100,100])
axs[0, 1].set_title('Parameter c vs a')

axs[1, 0].hist2d(par_b, par_c, bins=[100,100])
axs[1, 0].set_title('Parameter c vs b')

axs[1, 1].hist(chi2_reduced, bins=100, range=[np.min(chi2_reduced),np.max(chi2_reduced)])
axs[1, 1].set_title('Reduced chi^2 distribution')

fig.show()

plt.savefig("my-plot.png")

# **************************************
  
input("hit Enter to exit")
