### Flattened script from MIMOSCA ###
# Original can be found on gitgub at asncd/MIMOSCA


import sys
import pandas as pd
import numpy.matlib
import numpy as np
import scipy
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import stats
from numpy import unravel_index

n = pd.read_csv(sys.argv[1], header=None, index_col = None, sep = '\t')[1]
num_virus = int(sys.argv[2]) # Number of viruses in pool
nums = int(sys.argv[3]) # gridsize for performing lambda and alpha search

min_moi = 0.05
max_moi = 2

# Maximum number of possible viruses per cell
maxk = 10 

#specify start and finishing MOI to search over, it is set to 0.1 and 3 here
mois=np.linspace(min_moi, max_moi, nums)

#specify start and finishing detection probability to search over, it is set to 0.1 and 0.99 here
detects=np.linspace(0.1, 0.99, nums)

#initialize search array
LL=np.zeros((nums,nums))

#loop through square grid of different poission parameters and detection probabilities
for i in range(nums):
    for m in range(nums):
        #current parameter guesses
        moi_guess=mois[i]
        detect_guess=detects[m]
          
        #initialize possion distribution with current guess    
        pdf=scipy.stats.poisson.pmf(k=range(maxk),mu=moi_guess)
        
        #Zero truncation and renormalization
        pdf[0]=0
        pdf=np.divide(pdf,np.sum(pdf))

        
        #get probabilities after convolving with binomial distribution
        zibpdf=np.zeros((maxk,1))
        for k in range(maxk):
            pf=0
            for j in np.arange(k,maxk):
                pf+=pdf[j]*scipy.stats.binom.pmf(k,j,detect_guess)
            zibpdf[k]=pf
        
        #evaluate log likelihood after multiplying with observed values
        ll=1.0
        for k in range(len(n)):
            ll+=n[k]*np.log(zibpdf[k])
        LL[i,m]=ll       

#Log likelihood vs. paramter space
#plt.contour(np.round(detects,2),np.round(mois,2),LL,400,cmap='magma')
#plt.xlabel('Detection Probability')
#plt.ylabel('MOI')
#plt.savefig(pathout+'/_ph_moi_LL.pdf')

#Find paramters that maximize the log likelihood
final_tuple=unravel_index(LL.argmax(), LL.shape)
moi_guess=mois[final_tuple[0]]
detect_guess=detects[final_tuple[1]]
print(moi_guess,detect_guess)

#Create expected probability distribution given these paramters
pdf=scipy.stats.poisson.pmf(range(maxk),moi_guess)
pdf[0]=0
pdf=np.divide(pdf,np.sum(pdf))

zibpdf=np.zeros((maxk,1))

for k in range(maxk):
    pf=0
    for m in np.arange(k,maxk):
        pf+=pdf[m]*scipy.stats.binom.pmf(k,m,detect_guess)
    zibpdf[k]=pf
    
zibpdf_nocorrect=zibpdf.copy()

# Plot expected vs observed probabilities
#birthday problem correction, calculate the probability of drawing two of the same virus copies...
#in this example the contribution is negligible since the majoirty of cells only our observed to have one guide
top2bottom=np.arange(2,maxk,1)[::-1]
for k in top2bottom:
    delta_prob=zibpdf[k]*(1-scipy.stats.poisson.pmf(0,np.divide(scipy.misc.comb(k,2),num_virus)))
    zibpdf[k-1]+=delta_prob[0]
    zibpdf[k]-=delta_prob[0]

for i in range(0,len(n)):
    print str(i) + '\t' + str(float(n[i])/np.sum(n)) + '\t' + str(zibpdf[i][0])

plt.plot(range(maxk),np.cumsum(zibpdf),label='expected',alpha=0.75)
plt.plot(range(len(n)),np.cumsum(np.divide(1.0*n,np.sum(n))),label='observed',alpha=0.75)
plt.xlabel('MOI')
plt.ylabel('Cumulative Probability')
plt.legend()
plt.ylim([0,1.0])
plt.xlim([0,8])
plt.savefig('expected_vs_observed_prob.png')
plt.clf()



#plt.plot(range(maxk),np.log10(zibpdf),label='birthday corrected',alpha=0.75)
#plt.xlim([0,8])
#plt.ylim([-10,0])
#plt.xlabel('Number of sgRNAs')
#plt.ylabel('Log10(probability)')
