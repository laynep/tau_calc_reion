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

print tau.tau_calculator(x0)
