#This evaluates the limiting magnitude of the galaxy luminosity function of http://arxiv.org/pdf/1506.01035v2.pdf that matches the number density of galaxies of http://iopscience.iop.org/article/10.1088/0004-637X/813/1/54/pdf at redshift z=6.
import tau_calculator as tau
import numpy as np
import tau_parameters as params

import scipy.optimize as opt

schecter_fname1 = params.directory + 'schecter_params/schecter_params_om=0.30_h0=0.70_sig8=0.82.txt'
schecter = np.loadtxt(schecter_fname1)
z_list1 = schecter[:,0]
phi_list1 = schecter[:,1]
m_list1 = schecter[:,2]
alpha_list1 = schecter[:,3]
muv_list1 = schecter[:,4]

schecter_fname2 = params.directory + 'schecter_params/schecter_params_Bouwens.txt'
schecter = np.loadtxt(schecter_fname2)
z_list2 = schecter[:,0]
phi_list2 = schecter[:,1]
m_list2 = schecter[:,2]
alpha_list2 = schecter[:,3]
muv_list2 = schecter[:,4]

f_esc = tau.f_esc_funct('Power', [1.0,0.0])

z_match = 6.0

def emiss1():
    return tau.ioniz_emiss(z=z_match, f_esc=f_esc, photon_norm_factor=25.2,z_list=z_list1,phi_list=phi_list1,m_list=m_list1,alpha_list=alpha_list1,muv_list=muv_list1)

def emiss2(Muv):
    muv_list = Muv*np.ones(len(z_list2)) #Assume constant
    return tau.ioniz_emiss(z=z_match, f_esc=f_esc, photon_norm_factor=25.2,z_list=z_list2,phi_list=phi_list2,m_list=m_list2,alpha_list=alpha_list2,muv_list=muv_list)

Muv = opt.newton(lambda Muv: emiss1() - emiss2(Muv),-10.0)

print "schecter_fname1:", schecter_fname1
print "schecter_fname2:", schecter_fname2
print "ionizing emissitivity:", emiss1()
print "matching when M_{UV,SF}=", Muv
