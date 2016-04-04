"""Emcee sampler for importance sampling tau with pre-defined galaxy luminosity function parameters."""
import tau_calculator as tau
import emcee
import numpy as np

import importlib

#Get params file from command line
import sys
sys.path.append('./parameters/')
cline_args = sys.argv
modname = str(sys.argv[1]).strip()
if modname[-3:]==".py":
    modname=modname[:-3]

params = importlib.import_module(modname)


nwalkers = 300
niterations = 2000
nthreads = 16

save_fname = params.save_fname
p0_random = params.p0_random
p0_file = params.p0_file


if p0_random:
    p0=[]
    for walker in xrange(nwalkers):
        p0_here = []

        if params.ion_model=="Standard":

            if tau.f_esc_flag == "Power":
                p0_here.append(np.random.rand()) #f_8
                p0_here.append(np.random.rand()*4.0) #\alpha
            elif tau.f_esc_flag == "Linear":
                p0_here.append(np.random.rand()) #f_8
                p0_here.append(np.random.rand()) #slope
            elif tau.f_esc_flag == "Polint":

                p0_here.append(np.random.rand()) #f_3
                p0_here.append(np.random.rand()) #f_6
                p0_here.append(np.random.rand()) #f_9
                p0_here.append(np.random.rand()) #f_12

        elif params.ion_model=="Nonparametric":

            p0_here.append(np.random.rand()*1e-19) #emiss_3
            p0_here.append(np.random.rand()*1e-19) #emiss_6
            p0_here.append(np.random.rand()*1e-19) #emiss_9
            p0_here.append(np.random.rand()*1e-19) #emiss_12

        else:
            raise Exception('This ion model not implemented.')

        if 'C_HII' in params.nuisance:
            p0_here.append(np.random.rand()*(5.0-1.0)+1.0) #C_HII
        if 'xi_ion' in params.nuisance:
            p0_here.append(np.random.rand()*(26.0 - 24.0) + 24.0) #Photon norm

        if tau.data_type == "marg_cosmo":
            p0_here.append(np.random.normal(tau.globe.ombh2,tau.globe.ombh2*0.01)) #ombh2

            p0_here.append(np.random.normal(tau.globe.ommh2,tau.globe.ommh2*0.01)) #ommh2

        p0.append(p0_here)
    p0=np.array(p0)
else:
    p0 = np.loadtxt(p0_file)
    p0 = p0[-nwalkers:,1:]

ndim = len(p0[0])

sampler = emcee.EnsembleSampler(nwalkers, ndim, tau.logpost, threads=nthreads)

fname = open(save_fname, "w")
fname.close()

for result in sampler.sample(p0, iterations=niterations, storechain=False):
    position = result[0]
    fname = open(save_fname, "a")
    for k in range(position.shape[0]):
        fname.write("{0:4d} {1:s}\n".format(k, " ".join(str(posk) for posk in position[k])))
    fname.close()
