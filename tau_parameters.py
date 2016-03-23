#Parameters file

#Define f_esc parameters
f_esc_flag ='Power'
#f_esc_flag ='Polint'
#f_esc_flag ='Linear'

#Define the data type
#data_type = "tau_only"
data_type = "marg_cosmo"

data_file = 'data/total_TTTEEE_lowTEB.csv'
#data_file = 'data/total_margcosmo.txt'
#data_file = 'data/posterior_tauonly.txt'

use_lowfesc_const = True

#Nuisance params
#nuisance = ['C_HII','M_SF','xi_ion','dMdz']
nuisance = ['xi_ion']

#Where the main reionization directory is
directory = '/Users/laynep/work/reionization/importance_sampler/python_implementation/'
#directory = '/home/laynep/reion_importance/'

#Which set of Schecter params to use for the GLF
schecter_fname = 'schecter_params.txt'
#schecter_fname = 'schecter_params_Bouwens.txt'

#Save the chains to this file
save_fname = "chain_test.dat"

#Use random ICs for chains or load from file
p0_random = True
p0_file = "chain_power.dat"
