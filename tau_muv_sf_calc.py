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

def emiss1(z_match):
    return tau.ioniz_emiss(z=z_match, f_esc=f_esc, photon_norm_factor=25.2,z_list=z_list1,phi_list=phi_list1,m_list=m_list1,alpha_list=alpha_list1,muv_list=muv_list1)

def emiss2(z_match,Muv):
    muv_list = Muv*np.ones(len(z_list2)) #Assume constant
    return tau.ioniz_emiss(z=z_match, f_esc=f_esc, photon_norm_factor=25.2,z_list=z_list2,phi_list=phi_list2,m_list=m_list2,alpha_list=alpha_list2,muv_list=muv_list)

z_list = [5.0, 5.625, 6.25, 6.875, 7.5, 8.125, 8.75, 9.375, 10.0, 10.625, 11.25, 11.875, 12.5, 13.125, 13.75, 14.375, 15.0]
for z in z_list:
    Muv = opt.newton(lambda Muv: emiss1(z) - emiss2(z,Muv),-10.0)

    #print "schecter_fname1:", schecter_fname1
    #print "schecter_fname2:", schecter_fname2
    #print "ionizing emissitivity:", emiss1(z), emiss2(z,Muv)
    #print "matching when M_{UV,SF}=", Muv
    #print "z, Muv", z, Muv
    print "Muv", Muv

