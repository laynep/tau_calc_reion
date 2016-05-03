#Simple tool to calculate \tau
import tau_calculator as tau
import numpy as np


x0=[]
#x0.append(1.58160173e-11)
#x0.append(8.30644385e-12)
#x0.append(4.65758855e-11)
#x0.append(2.61772389e-12)
#x0.append(5.47992470e-16)
#x0.append(4.18676655e+00)

x0= [0.05, 0.1, 0.15, 0.2, 0.25, 27.0]
x0= [0.5*0.0833101816631, 0.5*0.17017923816, 0.5*0.294128625116, 0.5*0.337454428298, 24.000 ]

x0=np.array(x0)


tau0= tau.tau_calculator(x0)

print "This is tau"
print tau0[0]

print "This is loglike"
print tau.loglike(tau0[0],tau0[1],x0)

print "This is logpost"
print tau.logpost(x0)

print "Is params bad?", tau0[2]

#print "This is n_gamma:"
#for z in np.linspace(0,15,100):
#    f_esc = tau.f_esc_funct("Power", x0[0:2])
#    print z, tau.ioniz_emiss(z, f_esc, 25.2)


#print "This is n_H:", tau.comov_h_density(ombh2=0.022)

#print "this is Q"
#z_list = np.linspace(0,25.0,100)
#for q,z in zip(map(tau0[1], z_list), z_list):
#    print z,q
