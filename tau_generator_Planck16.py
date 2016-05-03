import numpy as np

mean = 0.05
std = 0.01

tau = np.random.normal(mean,std,size=100000)

np.savetxt('planck16_tau.dat',tau)
