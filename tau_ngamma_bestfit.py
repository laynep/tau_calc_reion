"""Calculates the escaped photon production rate density
for a given point in parameter space, usually the best-fit."""
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

z_list = np.linspace(6.0,15.0,10)
Mpc3_to_cm3 = (3.086e24)**3

x0 = params.bestfit[1:]

f_esc_params, c_hii, photon_norm_factor, dMSFdz, ombh2, ommh2 = tau.unpack(x0)
fesc = tau.f_esc_funct(params.f_esc_flag, f_esc_params)

print "this is fesc"
for z in z_list:
    print z, fesc.f_esc(z)

print "this is ngamma"
for z in z_list:
    print z, tau.ioniz_emiss(z, fesc, photon_norm_factor,dMSFdz)*Mpc3_to_cm3

print "This is tau",
tau, Q, bad = tau.tau_calculator(x0)
print tau
