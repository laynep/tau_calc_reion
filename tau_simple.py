#Simple tool to calculate \tau
import tau_calculator as tau
import numpy as np

x0=[]
x0.append(0.2) #f_6
x0.append(0.0) #\alpha
x0.append(3.0) #C_HII
x0.append(-10.0) #M_SF
x0.append(-0.35) #M_SF'
x0=np.array(x0)

tau0= tau.tau_calculator(x0)

z_list = np.linspace(0,25.0,100)

print "This is tau"
print tau0[0]

print "this is Q"
for q,z in zip(map(tau0[1], z_list), z_list):
    print z,q
