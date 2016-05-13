#Parameters file

#Define f_esc parameters
#f_esc_flag ='Power'
f_esc_flag ='Polint'
#f_esc_flag ='Linear'

#Define the data type
data_type = "tau_only"
#data_type = "marg_cosmo"

#data_file = 'data/total_TTTEEE_lowTEB.csv'
#data_file = 'data/total_TT_lowl_lensing.csv'
#data_file = 'data/total_TT_lowl_lensing_BAO.csv'
data_file = 'data/planck16_tau.csv'

use_lowfesc_const = True

#Nuisance params
nuisance = ['C_HII','xi_ion','dMSFdz']

#Where the main reionization directory is
#directory = '/Users/laynep/work/reionization/importance_sampler/python_implementation/'
directory = '/home/laynep/reion_importance/'

#Which set of Schecter params to use for the GLF
#schecter_fname = 'schecter_params/schecter_params.txt'
#schecter_fname = 'schecter_params/schecter_params_om=0.30_h0=0.70_sig8=0.80.txt'
#
schecter_fname = 'schecter_params/schecter_params_Bouwens.txt'
#schecter_fname = 'schecter_params/schecter_params_om=0.30_h0=0.70_sig8=0.82.txt'

#Save the chains to this file
save_fname = "chain_polint_planck16_Bouwens_xiion_dMSFdz.dat"

#Use random ICs for chains or load from file
p0_random = True
p0_file = "chain_power.dat"

#Model for the ionizing emissitivity
#ion_model = "Nonparametric"
ion_model = "Standard"

#Require f_esc to be monotonically increasing with z for the polint model
f_esc_monotonic = True


#Best fit value
bestfit = [4.79775034e+00,2.29563003e-02,3.71701191e-02,4.73402690e-02,8.64285554e-02,1.95375749e+00,2.51625988e+01,-3.99293448e-01]
