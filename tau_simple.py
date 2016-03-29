#Simple tool to calculate \tau
import tau_calculator as tau
import numpy as np

x0=[]
x0.append(2e-1)
x0.append(0.0)
x0.append(25.2)
x0.append(0.022)
x0.append(0.15)
x0=np.array(x0)

tau0= tau.tau_calculator(x0)

z_list = np.linspace(0,25.0,100)

print "This is tau"
print tau0[0]

print "This is loglike"
print tau.loglike(tau0[0],tau0[1],x0)

#print "this is Q"
#for q,z in zip(map(tau0[1], z_list), z_list):
#    print z,q
