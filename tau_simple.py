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

x0=[0.2, 0.00, 3.0]

#x0= [1e-15, 1e-20, 5e-22, 1e-21, 1e-25, 1e-25, 3.0]
#x0= [1e-15, 1e-15, 1e-20, 1e-20, 1e-20, 1e-20, 3.0]

x0= [4.17250233273e-16, 1.77888086156e-17, 2.75708251745e-19, 2.79253311941e-20, 5.50774139926e-20, 1.68754049981e-30, 1.41738459482]

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

print "this is Q"
z_list = np.linspace(0,25.0,100)
for q,z in zip(map(tau0[1], z_list), z_list):
    print z,q
