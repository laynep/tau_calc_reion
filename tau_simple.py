#Simple tool to calculate \tau
import tau_calculator as tau
import numpy as np

x0=[]
x0.append(3.73604206e-01)
x0.append(2.13735353e+00)
x0.append(2.59957181e+01)
x0.append(2.11533863e-02)
x0.append(1.40088349e-01)
x0=np.array(x0)

tau0= tau.tau_calculator(x0)

z_list = np.linspace(0,25.0,100)

print "This is tau"
print tau0[0]

print "this is Q"
for q,z in zip(map(tau0[1], z_list), z_list):
    print z,q
