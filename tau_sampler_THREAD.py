"""Emcee sampler for importance sampling tau with pre-defined galaxy luminosity function parameters."""
import tau_calculator as tau
import emcee
import numpy as np
import scipy.optimize as opt

import importlib

#Get params file from command line
import sys
sys.path.append('./parameters/')
cline_args = sys.argv
modname = str(sys.argv[1]).strip()
if modname[-3:]==".py":
    modname=modname[:-3]

params = importlib.import_module(modname)

save_fname = params.save_fname

nwalkers = 300
niterations = 700
nthreads = 16

p0_random = params.p0_random
p0_file = params.p0_file

def get_IC():

    p0_here = []
    
    if params.ion_model=="Standard":
    
        if tau.f_esc_flag == "Power":
            p0_here.append(np.random.rand()*0.1) #f_8
            p0_here.append(np.random.rand()*4.0) #\alpha
        elif tau.f_esc_flag == "Linear":
            p0_here.append(np.random.rand()) #f_8
            p0_here.append(np.random.rand()) #slope
        elif tau.f_esc_flag == "Polint":
    
            if 'planck' in params.data_file:
                factor = 0.25
            else:
                factor = 1.0
    
            p0_here.append(factor*(np.random.rand()*0.1)) #f_3
            p0_here.append(factor*(np.random.rand()*(0.2-0.1)+0.1)) #f_6
            p0_here.append(factor*(np.random.rand()*(0.3-0.2)+0.2)) #f_9
            p0_here.append(factor*(np.random.rand()*(0.4-0.3)+0.3)) #f_12
    
    
    elif params.ion_model=="Nonparametric":
    
        p0_here.append(np.random.rand()*(-17.2+21.2)-21.2) #emiss_3
        p0_here.append(np.random.rand()*(-17.2+21.2)-21.2) #emiss_6
        p0_here.append(np.random.rand()*(-17.2+21.2)-21.2) #emiss_9
        p0_here.append(np.random.rand()*(-17.2+21.2)-21.2) #emiss_12
        p0_here.append(np.random.rand()*(-17.2+21.2)-21.2) #emiss_15
        p0_here.append(np.random.rand()*(-17.2+21.2)-21.2) #emiss_18
    
        #p0_here.append(np.random.rand()*1e-19) #emiss_3
        #p0_here.append(np.random.rand()*1e-19) #emiss_6
        #p0_here.append(np.random.rand()*1e-19) #emiss_9
        #p0_here.append(np.random.rand()*1e-19) #emiss_12
        #p0_here.append(np.random.rand()*1e-19) #emiss_15
        #p0_here.append(np.random.rand()*1e-19) #emiss_18
    
    else:
        raise Exception('This ion model not implemented.')
    
    if 'C_HII' in params.nuisance:
        p0_here.append(np.random.rand()*(5.0-1.0)+1.0) #C_HII
    if 'xi_ion' in params.nuisance:
        p0_here.append(np.random.rand()*(26.0 - 24.0) + 24.0) #Photon norm

    if 'dMSFdz' in params.nuisance:
        p0_here.append(np.random.rand()*(0.0 + 0.5) - 0.5) #MSF evolve with z
    
    if tau.data_type == "marg_cosmo":
        p0_here.append(np.random.normal(tau.globe.ombh2,tau.globe.ombh2*0.01)) #ombh2
    
        p0_here.append(np.random.normal(tau.globe.ommh2,tau.globe.ommh2*0.01)) #ommh2

    return np.array(p0_here)


if p0_random:
    p0=[]
    for walker in xrange(nwalkers):

        cycling = True
        while cycling:
            p0_here = get_IC()
            post =  tau.logpost(p0_here)
            if post>-1e60:
                #print "IC set", walker
                p0.append(p0_here)
                cycling = False

    p0=np.array(p0)
else:
    p0 = np.loadtxt(p0_file)
    p0 = p0[-nwalkers:,1:]

print "------------------"
print "ICs have been set."
print "------------------"

ndim = len(p0[0])

sampler = emcee.EnsembleSampler(nwalkers, ndim, tau.logpost, threads=nthreads)

fname = open(save_fname, "w")
fname.close()

for result in sampler.sample(p0, iterations=niterations, storechain=False):
    fname = open(save_fname, "a")

    #result[0] = chain
    #result[1] = like

    for elmn in zip(result[1],result[0]):
        fname.write("%s " % str(elmn[0]))
        for k in list(elmn[1]):
            fname.write("%s " % str(k))
        fname.write("\n")
