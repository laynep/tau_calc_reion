#Stupid post processor

import numpy as np
import tau_calculator as tau
from sklearn import tree

import warnings
warnings.filterwarnings("ignore",category=DeprecationWarning)

directory = '/home/laynep/reion_importance/'
#name = 'chain_power_lensing_BAO_Bouwens_tauonly.dat'
#name = 'chain_power_TTTEEE_Bouwens_tauonly.dat'
#name = 'chain_power_TTTEEE_margcosmo.dat'
name = 'chain_ngamma_lensing_BAO_tauonly2.dat'
burnin = 100000


trainsize = 500

data = np.loadtxt(directory+name)

good_data = []

post_list = []

for dat in data[burnin:burnin+trainsize]:
    x0 = dat[1:]

    x0=np.array(x0)

    post = tau.logpost(x0)
    #post_list.append(post)


    #tau0= tau.tau_calculator(x0)
    #print "this is tau:", tau0[0]
    print "this is x0:", x0
    print "this is post:", post
    #bad = tau0[2]


    if post > 1e-30:
        post_list.append(1)
        #good_data.append(dat)
    else:
        post_list.append(0)

post_list = np.array(post_list)

#Train a really stupid tree
clf = tree.DecisionTreeRegressor(max_depth=2)     
clf = clf.fit(data[burnin:burnin+trainsize,1:],post_list)

for dat in  data:
    x0 = dat[1:]

    if clf.predict(x0)>0.25:
        good_data.append(dat)

    #print "testing here:",x0
    #print "learned post:", clf.predict(x0)
    #print "real post:", tau.logpost(x0)

np.savetxt(directory+name[:-4]+'_postprocess.dat',np.array(good_data))
