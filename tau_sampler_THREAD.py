#Emcee sampler
import tau_calculator as tau
import emcee
import numpy as np

nwalkers = 300
niterations = 2000

nthreads = 16

save_fname = "chain_power2.dat"
#save_fname = "chain_polint.dat"
#save_fname = "chain_polint_margcosmo.dat"
#save_fname = "chain_margcosmo.dat"

p0_random = False
p0_file = "chain_power.dat"

if tau.f_esc_flag == "Power":
    ndim = 6
elif tau.f_esc_flag == "Polint":
    ndim = 8

if tau.data_type == "marg_cosmo":
    ndim += 2

if p0_random:
    p0=[]
    for walker in xrange(nwalkers):
        p0_here = []
        if tau.f_esc_flag == "Power":
            p0_here.append(np.random.rand()) #f_6
            p0_here.append(np.random.rand()*4.0) #\alpha
        elif tau.f_esc_flag == "Polint":
            p0_here.append(np.random.rand()) #f_3
            p0_here.append(np.random.rand()) #f_6
            p0_here.append(np.random.rand()) #f_9
            p0_here.append(np.random.rand()) #f_12
        p0_here.append(np.random.rand()*(5.0-1.0)+1.0) #C_HII
        p0_here.append(np.random.rand()*(-9.5+11.0) - 11.0) #M_SF
        p0_here.append(np.random.rand()*(-0.30 + 0.4) - 0.4) #M_SF'
        p0_here.append(np.random.rand()*(26.0 - 24.0) + 24.0) #Photon norm
    
        if tau.data_type == "marg_cosmo":
            p0_here.append(np.random.normal(tau.globe.ombh2,tau.globe.ombh2*0.01)) #ombh2
    
            p0_here.append(np.random.normal(tau.globe.ommh2,tau.globe.ommh2*0.01)) #ommh2
    
        p0.append(p0_here)
    p0=np.array(p0)
else:
    p0 = np.loadtxt(p0_file)
    p0 = p0[-nwalkers:,1:]

sampler = emcee.EnsembleSampler(nwalkers, ndim, tau.logpost, threads=nthreads)

fname = open(save_fname, "w")
fname.close()

for result in sampler.sample(p0, iterations=niterations, storechain=False):
    position = result[0]
    fname = open(save_fname, "a")
    for k in range(position.shape[0]):
        fname.write("{0:4d} {1:s}\n".format(k, " ".join(str(posk) for posk in position[k])))
    fname.close()
