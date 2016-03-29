import numpy as np
from scipy.integrate import ode
from scipy.integrate import quad
from scipy.interpolate import interp1d
from scipy.interpolate import NearestNDInterpolator
import pandas as pd

import mpmath as mpm
import sys

import tau_parameters as params


#Define the redshift evolution of f_esc, etc
f_esc_flag = params.f_esc_flag
data_type = params.data_type
directory = params.directory
schecter_fname = params.schecter_fname
data_file = params.data_file

class f_esc_funct():

    def __init__(self,f_esc_flag, f_esc_params):

        if f_esc_flag == "Power":
            self.f_esc = lambda z: self.f_esc_POWER(z, f_esc_params)
        elif f_esc_flag == "Linear":
            self.f_esc = lambda z: self.f_esc_LINEAR(z, f_esc_params)
        elif f_esc_flag == "Tanh":
            self.f_esc = lambda z: self.f_esc_TANH(z, f_esc_params)
        elif f_esc_flag == "Polint":

            self.z_pos = np.linspace(3,12,4)
            self.f_polint = interp1d(self.z_pos, f_esc_params[0:4], kind='cubic')
            self.f_esc = lambda z: self.f_esc_POLINT(z, f_esc_params)
        else:
            raise Exception('This f_esc_flag not implemented.')


    def f_esc_POWER(self, z, f_esc):
        f8 = f_esc[0]
        alpha = f_esc[1]

        out = f8*((1.0+z)/9.0)**alpha

        return np.min([np.max([0.0, out]), 1.0])

    def f_esc_LINEAR(self, z, f_esc):
        f8 = f_esc[0]
        slope = f_esc[1]

        out = f8 + slope*(z-8.0)

        return np.min([np.max([0.0, out]), 1.0])

    def f_esc_POLINT(self, z, f_esc):

        if z< np.min(self.z_pos):
            return f_esc[0]
        elif z>np.max(self.z_pos):
            return f_esc[-1]
        else:
            return np.min([np.max([0.0, self.f_polint(z)]), 1.0])

    def f_esc_TANH(self, z,f_esc):
        f_inf = f_esc[0]
        f_0 = f_esc[1]
        z_half = f_esc[2]
        dz = f_esc[3]

        f_out = 0.5*(f_inf-f_0)*np.tanh((z-z_half)/dz) + 0.5*(f_inf+f_0)
        return np.min([np.max([0.0, f_out]), 1.0])



#Fixed params
class global_params(object):

    def __init__(self,schecter_fname):
        """Some general set up required."""

        self.m_bright = -35.0
        #self.photon_norm_factor = 25.2

        self.T0_IGM = 2.0e4
        self.X_h = 0.747
        self.HeII_z = 3.0

        self.schecter_fname = directory + schecter_fname
        schecter = np.loadtxt(self.schecter_fname)
        self.z_list = schecter[:,0]
        self.phi_list = schecter[:,1]
        self.m_list = schecter[:,2]
        self.alpha_list = schecter[:,3]
        self.muv_list = schecter[:,4]

        self.mpc_to_cm = 3.086e24
        self.mass_h = 1.674e-24
        self.G_newton = 6.674e-8
        self.sigT = 6.6524e-25
        self.c = 2.998e10

        #Planck2015 best fit
        self.h = 0.6727
        self.ommh2 = 0.02225+0.1198
        self.ombh2 = 0.02225

        #Scorch I
        ##Can also be iterated over
        #self.h = 0.7
        ##self.ommh2 = 0.27*self.h**2
        #self.ommh2 = 0.3*self.h**2
        #self.ombh2 = 0.045*self.h**2

        if data_type == "tau_only":
            self.tau_post_fname =directory + data_file
            try:
                tau_pdf = np.loadtxt(self.tau_post_fname)
            except:
                raise Exception('The data file is not available.')
            self.tau_list = tau_pdf[:,1]
            self.taupdf_list = tau_pdf[:,0]
        elif data_type == "marg_cosmo":
            try:
                self.data_margcosmo = pd.read_csv(directory + data_file)
            except:
                raise Exception('The data file is not available.')


globe = global_params(schecter_fname)

def unpack(x):

    y = np.array(x)

    if f_esc_flag=="Power":
        offset=2
    elif f_esc_flag=="Linear":
        offset=2
    elif f_esc_flag=="Polint":
        offset=4
    elif f_esc_flag=="Tanh":
        offset=4
    else:
        raise Exception('This f_esc_flag not supported.')

    f_esc_params = x[0:offset]

    index = offset

    if 'C_HII' in params.nuisance:
        c_hii = x[index]
        index += 1
    else:
        c_hii = 3.0
    if 'xi_ion' in params.nuisance:
        photon_norm_factor = x[index]
        index += 1
    else:
        photon_norm_factor = 25.2

    if data_type == "tau_only":
        #Can get these either from x or from fixed
        ombh2 = globe.ombh2
        ommh2 = globe.ommh2
    elif data_type == "marg_cosmo":
        ombh2 = x[index]
        index += 1
        ommh2 = x[index]
    else:
        raise Exception('This data_type not supported.')

    return f_esc_params, c_hii, photon_norm_factor, ombh2, ommh2

def mag_to_lumin(mag):
    """Converts AB magnitude to luminosity."""
    conv_factor = 4.345e20
    return float(conv_factor*10.0**(-mag/2.5))

def schecter_params(z,z_list=False,phi_list=False,m_list=False,alpha_list=False,muv_list=False):
    """Get the Schecter parameters as a function of redshift from file."""

    logical = type(False)
    if (type(z_list) == logical):
        z_list = globe.z_list
    if (type(phi_list)  == logical):
        phi_list = globe.phi_list
    if (type(m_list)  == logical):
        m_list = globe.m_list
    if (type(alpha_list)  == logical):
        alpha_list = globe.alpha_list
    if (type(muv_list)  == logical):
        muv_list = globe.muv_list

    #Use range limits if z extends beyond the interpolation list
    if z<np.min(z_list):
        phi_star, m_star, alpha_star, muv_star = phi_list[0],m_list[0],alpha_list[0],muv_list[0]
    elif z>np.max(globe.z_list):
        phi_star, m_star, alpha_star, muv_star = phi_list[-1], m_list[-1], alpha_list[-1],muv_list[-1]

    else:
        #Interpolate
        phi, m, alpha, muv = map(lambda x: interp1d(z_list,x,kind='cubic'),
                            [phi_list, m_list, alpha_list, muv_list])

        phi_star, m_star, alpha_star, muv_star = phi(z), m(z), alpha(z), muv(z)

    L_star = mag_to_lumin(m_star)
    L_SF = mag_to_lumin(muv_star)

    #Convert phi into cgs --- currently Mpc^-3
    phi_star *= globe.mpc_to_cm**(-3.0)

    return float(phi_star), float(L_star), float(alpha_star), float(L_SF)

def ioniz_emiss(z, f_esc, photon_norm_factor,z_list=False,phi_list=False,m_list=False,alpha_list=False,muv_list=False):

    L_bright = mag_to_lumin(globe.m_bright)

    phi_star, L_star, alpha, L_faint = schecter_params(z,z_list,phi_list,m_list,alpha_list,muv_list)

    prefactor = -2.5/np.log(10.0)
    prefactor *= 10.0**(photon_norm_factor)
    prefactor *= f_esc.f_esc(z)*phi_star*L_star

    if alpha<-2.0:
        #Close approximation
        return prefactor*(alpha+1.0)**(-1.0)*(L_faint/L_star)**(alpha+1.0)
    else:
        return prefactor*(mpm.gammainc(1.0+alpha,L_bright/L_star) - mpm.gammainc(1.0+alpha,L_faint/L_star))

def t_recomb(z, c_hii):
    prefactor = 0.93*1e9*365.25*24.0*60.0*60.0 #Gyr to sec
    c_term = (c_hii/3.0)**(-1.0)
    t_term = (globe.T0_IGM/2.0e4)**0.7
    z_term = ((1.0+z)/7.0)**(-3.0)
    return prefactor*c_term*t_term*z_term

def comov_h_density(ombh2):
    prefactor = (3.0/8.0/np.pi/globe.G_newton)*1e4*(1e5/globe.mpc_to_cm)**2

    return prefactor*(globe.X_h/globe.mass_h)*ombh2

def hubble(z, ommh2):
    h2 = 1e2**2 *(ommh2*(1.0+z)**3.0)
    return np.sqrt(h2)*1e5/globe.mpc_to_cm

def electron_from_helium(z):
    if z > globe.HeII_z:
        return 1.0
    else:
        return 2.0


def tau_calculator(x):
    """Calculate the Thomson optical depth from the galaxy reionization parameters."""
    f_esc_params, c_hii, photon_norm_factor, ombh2, ommh2 = unpack(x)

    #Initiate f_esc object
    f_esc = f_esc_funct(f_esc_flag, f_esc_params)

    def derivs(z, Q):
        """The Q'(z0) = derivs(z, Q)"""

        hub = hubble(z, ommh2)
        trec = t_recomb(z, c_hii)
        nh = comov_h_density(ombh2)
        ion = ioniz_emiss(z, f_esc, photon_norm_factor)

        return -1.0*(1.0/hub/(1.0+z))*(ion/nh/(1.0+z)**3.0 - Q/trec)

    #Solve the ODE for the volume filling factor Q
    solver = ode(derivs).set_integrator("vode")

    #ICs
    Q0 = 2.1979099400481504e-004 #From CAMB, connects to residual ionization from recombination
    z0= 25.0
    solver.set_initial_value(Q0,z0)

    nsteps = 100.0
    dz = -z0/nsteps
    z_ode = []
    Q_ode = []
    while solver.successful() and solver.t > 0.0 and solver.y<1.0 and solver.y>0.0:
        solver.integrate(solver.t + dz)
        z_ode.append(solver.t)
        Q_ode.append(solver.y)
        #print solver.t, solver.y


    #Set up interpolating function
    if len(z_ode) > 3 and len(Q_ode)>3:
        naive_Q = interp1d(np.array(z_ode),np.array(Q_ode).flatten(),kind='cubic')
    elif len(z_ode)>=2 and len(Q_ode)>=2:
        naive_Q = interp1d(np.array(z_ode),np.array(Q_ode).flatten(),kind='linear')
    else:
        raise Exception('The Q integrator took only one step.')

    def Q_of_z(z):
        #Volume filling factor
        if z<np.min(z_ode):
            Q=1.0
        elif z > np.max(z_ode):
            Q=Q_ode[0]
        else:
            Q = np.max([np.min([1.0,naive_Q(z)]),0.0])

        return Q

    #Integrate Q to get \tau
    def tau_integrand(z):
        prefactor = globe.sigT*globe.c*1e-2*comov_h_density(ombh2)*globe.mpc_to_cm/1e5/np.sqrt(ommh2)

        integrand = Q_of_z(z)*np.sqrt(1.0+z)*(1.0+(1.0-globe.X_h)/(4.0*globe.X_h)*electron_from_helium(z))

        return prefactor*integrand

    #Calculate \tau
    (tau, dtau) = quad(tau_integrand, 0.0, z0)

    if dtau > 1e-4:
        raise Exception('tau integrator has large error.')

    return tau, Q_of_z


def logprior(x):
    f_esc_params, c_hii, photon_norm_factor, ombh2, ommh2 = unpack(x)

    bad = False
    const = 1.0  #Unnormalized prior

    #Prior ranges
    if f_esc_flag=='Power':
        if f_esc_params[0]<0.0 or f_esc_params[0]>1.0:
            bad = True
        elif f_esc_params[1]<0.0 or f_esc_params[1]>4.0:
            bad = True
    elif f_esc_flag=='Linear':
        if f_esc_params[0]<0.0 or f_esc_params[0]>1.0:
            bad = True
        elif f_esc_params[1]<0.0 or f_esc_params[1]>1.0:
            bad = True
    elif f_esc_flag=="Polint":
        if any(f_esc_params<0.0) or any(f_esc_params>1.0):
            bad = True
    elif f_esc_flag == "Tanh":
        if f_esc_params[0]<0.0 or f_esc_params[0]>1.0:
            bad = True
        if f_esc_params[1]<0.0 or f_esc_params[1]>1.0:
            bad = True
        if f_esc_params[2]<0.01 or f_esc_params[2]>1.0:
            bad = True
        if f_esc_params[3]<0.0 or f_esc_params[3]>25.0:
            bad = True
    else:
        raise Exception('This f_esc_flag not implemented')

    if c_hii<1.0 or c_hii>5.0:
        bad = True

    if photon_norm_factor <24.0 or photon_norm_factor>26.0:
        bad = True

    if bad:
        return -np.inf, False
    else:
        return const, True

def loglike(tau,x):
    f_esc_params, c_hii, photon_norm_factor, ombh2, ommh2 = unpack(x)

    if params.use_lowfesc_const:
        #Initiate f_esc object
        f_esc = f_esc_funct(f_esc_flag, f_esc_params)
        #Implement 2\sigma constraint from Boutsia et al 2011
        if f_esc.f_esc(3.3)>0.10:
            return -np.inf

    if data_type=="tau_only":

        if tau<np.min(globe.tau_list) or tau > np.max(globe.tau_list):
            return -np.inf
        else:

            #Likelihood just from the LCDM tau values with cosmo params fixed
            post_funct = interp1d(globe.tau_list, globe.taupdf_list,kind='cubic')

            #print "This is post_funct", post_funct(tau)
            post = np.max([0.0,post_funct(tau)])

            if post<=0.0:
                return -np.inf
            else:
                return np.log(post)

    elif data_type == "marg_cosmo":

        if tau<np.min(globe.data_margcosmo['tau']) or tau > np.max(globe.data_margcosmo['tau']):
            #print "Caught by tau"
            return -np.inf
        elif ombh2<np.min(globe.data_margcosmo['omegabh2']) or ombh2 > np.max(globe.data_margcosmo['omegabh2']):
            #print "Caught by ombh2"
            return -np.inf
        elif ommh2<np.min(globe.data_margcosmo['omegamh2']) or ommh2 > np.max(globe.data_margcosmo['omegamh2']):
            #print "Caught by ommh2"
            return -np.inf
        else:

            #Only really feasible to do nearest neighbors interpolation here
            #Here we assume flat priors in LCDM so that like=post
            post_funct = NearestNDInterpolator(
                np.array(globe.data_margcosmo[['omegabh2','omegamh2','tau']]),
                np.array(globe.data_margcosmo['like']))

            post = np.max([0.0,post_funct([ombh2,ommh2,tau])])

            #print "this is post:", post

            if post<=0.0:
                return -np.inf
            else:
                #return np.log(post)
                return post #Like is already log-like

    else:
        raise Exception('This data type not supported.')


def logpost(x):

    f_esc_params, c_hii, photon_norm_factor, ombh2, ommh2 = unpack(x)

    prior, success = logprior(x)
    if not success:
        return prior

    #Likelihood

    tau, Q = tau_calculator(x)

    #print "this is tau", tau, Q(5.9)

    #Implement 2sig constraint from 1411.5375 as a step function
    if Q(z=5.9)<0.84:
        like = -np.inf

    else:
        like = loglike(tau, x)

    return prior + like
