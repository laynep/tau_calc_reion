#Tool to find the reionization parameters that provide the closest fit to a given \tau
#which we provide, assuming that it's the best-fit value
import tau_calculator as tau
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

if params.data_file == 'data/total_TTTEEE_lowTEB.csv':
    tau_ref = 0.079
elif params.data_file == 'data/total_TT_lowl_lensing_BAO.csv':
    tau_ref = 0.067
elif params.data_file == 'data/planck16_tau.csv':
    tau_ref = 0.055
else:
    raise Exception('Set target tau here.')

print "this is tau_ref:", tau_ref



def closest_tau(x):

    tau0= tau.tau_calculator(x)

    #print "----------------------"
    #print "This is x:", x
    print "This is tau:", tau0[0]

    return (tau_ref - tau0[0])**2

#Optimize
#x0= [-19.0, -19.0, -19.0, -19.0, -19.0, -19.0, 3.0]
#x0 = [0.039,0.041,0.51,0.793,25.05]
#x0=np.array(x0)

#print opt.minimize(closest_tau,x0)
#print opt.basinhopping(closest_tau,x0,T=10.0)


#Loop through chain and find closest
data = np.loadtxt(params.save_fname)

data = data[data[:,0].argsort()[::-1]] #Sort by likelihood (1st colmn)
maxlike = data[0,0]

#Calculate \tau for degenerate ones
diff = np.inf
bestfit = data[0]
for dat in data:
    print dat[0], maxlike
    if dat[0] < maxlike:
        break
    else:
        newdiff = closest_tau(dat[1:])
        if newdiff < diff:
            diff = newdiff
            bestfit = dat
            print "new best fit:", bestfit


print "This is the likelihood:", params.data_file
print "This is the best fit value:", bestfit


save_fname = params.data_file[:-4]+'_bestfit.dat'
np.savetxt(save_fname, bestfit)
#fname = open(save_fname, "a")
#
#for k in bestfit:
#    fname.write("%s " % str(k))
