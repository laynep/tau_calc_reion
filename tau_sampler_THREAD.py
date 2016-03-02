#Emcee sampler
import tau_calculator as tau
import emcee
import numpy as np

ndim = 5
nwalkers = 200
niterations = 1000

nthreads = 12

save_fname = "chain.dat"

p0=[]
for walker in xrange(nwalkers):
    p0_here = []
    p0_here.append(np.random.rand()) #f_6
    p0_here.append(np.random.rand()*4.0) #\alpha
    p0_here.append(np.random.rand()*(5.0-1.0)+1.0) #C_HII
    p0_here.append(np.random.rand()*(-9.5+11.0) - 11.0) #M_SF
    p0_here.append(np.random.rand()*(-0.031 + 0.035) - 0.035) #M_SF'

    if tau.data_type == "marg_cosmo":
        p0_here.append(np.random.normal(tau.globe.ombh2,tau.globe.ombh2*0.01)) #ombh2

        p0_here.append(np.random.normal(tau.globe.ommh2,tau.globe.ommh2*0.01)) #ommh2

    p0.append(p0_here)
p0=np.array(p0)

sampler = emcee.EnsembleSampler(nwalkers, ndim, tau.logpost, threads=nthreads)

fname = open(save_fname, "w")
fname.close()

for result in sampler.sample(p0, iterations=niterations, storechain=False):
    position = result[0]
    fname = open(save_fname, "a")
    for k in range(position.shape[0]):
        fname.write("{0:4d} {1:s}\n".format(k, " ".join(str(posk) for posk in position[k])))
    fname.close()
