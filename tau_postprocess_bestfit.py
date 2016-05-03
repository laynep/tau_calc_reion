#Post processor that adds the likelihood to the data file

import numpy as np
import tau_calculator as tau

import importlib

#Get params file from command line
import sys
sys.path.append('./parameters/')
cline_args = sys.argv
modname = str(sys.argv[1]).strip()
if modname[-3:]==".py":
    modname=modname[:-3]

params = importlib.import_module(modname)

directory = '/home/laynep/reion_importance/'

name = params.save_fname

data = np.loadtxt(directory+name)


burnin = 0

new_data = []

#Just loop through it like an idiot
for dat in data[burnin:]:
    x0 = dat[1:]

    x0=np.array(x0)

    post = tau.logpost(x0)
    #post_list.append(post)

    x0 = np.insert(x0,0,post)

    new_data.append(x0)


np.savetxt(directory+name[:-4]+'_likelihoods.dat',np.array(new_data))
