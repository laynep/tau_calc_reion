import numpy as np

mean = 0.055
std = 0.009

tau = np.random.normal(mean,std,size=100000)

np.savetxt('planck16_tau.dat',tau)
