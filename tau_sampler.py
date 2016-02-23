import tau_calculator as tau
import emcee
import numpy as np

#Emcee sampler
nthreads = 12


ndim = 5

nwalkers = 250
p0=[]
for walker in xrange(nwalkers):
    p0_here = []
    p0_here.append(np.random.rand()) #f_6
    p0_here.append(np.random.rand()*4.0) #\alpha
    p0_here.append(np.random.rand()*(5.0-1.0)+1.0) #C_HII
    p0_here.append(np.random.rand()*(-9.5+11.0) - 11.0) #M_SF
    p0_here.append(np.random.rand()*(-0.031 + 0.035) - 0.035) #M_SF'
    p0.append(p0_here)
p0=np.array(p0)


sampler = emcee.EnsembleSampler(nwalkers, ndim, tau.logpost, threads=nthreads)
pos, prob, state = sampler.run_mcmc(p0, 1000)

np.save(tau.directory + 'output_flatchain.npy', sampler.flatchain)
