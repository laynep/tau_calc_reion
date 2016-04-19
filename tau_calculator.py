import numpy as np
from scipy.integrate import ode
from scipy.integrate import quad
from scipy.interpolate import interp1d
from scipy.interpolate import NearestNDInterpolator
import pandas as pd

import mpmath as mpm
import sys

#Get params file from command line
import importlib
import sys
sys.path.append('./parameters/')
cline_args = sys.argv
modname = str(sys.argv[1]).strip()
if modname[-3:]==".py":
    modname=modname[:-3]
elif modname[-4:]==".pyc":
    modname=modname[:-4]

params = importlib.import_module(modname)

#Define the redshift evolution of f_esc, etc
f_esc_flag = params.f_esc_flag
data_type = params.data_type
directory = params.directory
schecter_fname = params.schecter_fname
data_file = params.data_file


class f_esc_funct():

    def __init__(self,f_esc_flag, f_esc_params):

        self.params = f_esc_params

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

        try:
            self.data = pd.read_csv(directory + data_file)
        except:
            raise Exception('The data file is not available.')

        if data_type == "tau_only":

            #Weighted 3D histogram normed to the posterior value
            nbins=40
            try:
                hist, bins = np.histogramdd(
                        np.array(self.data['tau']),
                        bins=nbins,
                        normed=True,
                        weights=np.array(self.data['weight']))
            except:
                try:
                    hist, bins = np.histogramdd(
                            np.array(self.data['tau']),
                            bins=nbins,
                            normed=True)
                except:
                    raise Exception('Could not build density estimator.')

        elif data_type == "marg_cosmo":

            #Weighted 3D histogram normed to the posterior value
            nbins=15
            try:
                hist, bins = np.histogramdd(
                        np.array(self.data[['omegabh2','omegamh2','tau']]),
                        bins=nbins,
                        normed=True,
                        weights=np.array(self.data['weight']))
            except:
                try:
                    hist, bins = np.histogramdd(
                            np.array(self.data[['omegabh2','omegamh2','tau']]),
                            bins=nbins,
                            normed=True)
                except:
                    raise Exception('Could not build density estimator.')

        else:
            raise Exception('This data type not supported.')

        hist = np.array(hist)

        def post_funct(x):
            """Log-like of posterior."""
            #ombh2,ommh2,tau = x[0], x[1], x[2] marg_cosmo
            #tau = x[0] tau_only

            loc=[]
            for index in xrange(len(x)):
                pos = np.searchsorted(bins[index],x[index])
                if pos <1 or pos == len(bins[index]):
                    return -np.inf
                else:
                    loc.append(pos)

            loc = [l-1 for l in loc]

            if data_type == "tau_only":
                post = hist[loc[0]]
            elif data_type == "marg_cosmo":
                post= hist[loc[0],loc[1],loc[2]]

            if post<=0.0:
                return -np.inf
            else:
                return np.log(post)

        self.post_funct = post_funct

        #DEBUG
        #if data_type == "marg_cosmo":
        #    print "This is post funct", post_funct([0.0224,0.144,0.040])
        #    print "This is post funct", post_funct([0.0224,0.144,0.050])
        #    print "This is post funct", post_funct([0.0224,0.144,0.066])
        #    print "This is post funct", post_funct([0.0224,0.144,0.079])
        #    print "This is post funct", post_funct([0.0224,0.144,0.09])
        #else:
        #    print "This is post funct", post_funct([0.040])
        #    print "This is post funct", post_funct([0.050])
        #    print "This is post funct", post_funct([0.066])
        #    print "This is post funct", post_funct([0.079])
        #    print "This is post funct", post_funct([0.09])
        #sys.exit()

        self.print_counter = 0


globe = global_params(schecter_fname)

def unpack(x):

    y = np.array(x)

    index = 0

    if params.ion_model=="Standard":

        if f_esc_flag=="Power":
            index+=2
        elif f_esc_flag=="Linear":
            index+=2
        elif f_esc_flag=="Polint":
            index+=4
        elif f_esc_flag=="Tanh":
            index+=4
        else:
            raise Exception('This f_esc_flag not supported.')

    elif params.ion_model=="Nonparametric":

        index +=6

    else:
        raise Exception('This ion model not supported.')

    f_esc_params = x[0:index]

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

    if params.ion_model=="Nonparametric":
        #This overrides the f_esc params and uses them directly as knot-spline approach for the ionizing emissitivity.


        if z<3.0:
            return f_esc.params[0]
        elif z >18.0:
            return f_esc.params[-1]
        else:
            #I know it's bad to interpolate at every step.  Sue me.
            z_list = np.linspace(3.0,18.0,6)
            n_gamma = interp1d(z_list,f_esc.params,kind='linear')

            return np.max([0.0,n_gamma(z)])

    elif params.ion_model=="Standard":

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

    else:
        raise Exception('This model for the ionizing emissitivity is not implemented.')

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

    bad = False

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
    if params.ion_model=="Standard":
        z0= 25.0
    elif params.ion_model=="Nonparametric":
        z0= 17.9
    else:
        raise Exception('This ion model not supported.')
    solver.set_initial_value(Q0,z0)

    nsteps = 100.0
    dz = -z0/nsteps
    z_ode = [z0]
    Q_ode = [Q0]
    while solver.successful() and solver.t > 0.0 and solver.y<1.0 and solver.y>0.0:
        solver.integrate(solver.t + dz)
        z_ode.append(solver.t)
        Q_ode.append(solver.y)
        #print "z and Q:", solver.t, solver.y

    if z_ode[-1]>17.0:
        #Reionization happens way too early
        bad = True

    #Set up interpolating function
    if len(z_ode) > 3 and len(Q_ode)>3:
        naive_Q = interp1d(np.array(z_ode),np.array(Q_ode).flatten(),kind='cubic')
    elif len(z_ode)>=2 and len(Q_ode)>=2:
        #Acts like instantaneous reioniz at z0, but could be at much higher z.
        bad = True
        naive_Q = interp1d(np.array(z_ode),np.array(Q_ode).flatten(),kind='linear')
    else:
        bad = True
        print "z_ode", z_ode
        print "Q_ode", Q_ode
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

    if np.abs(tau > 1e-6):
        if np.abs(dtau/tau) > 1e-2: #percent-level error
            if tau > 0.15:
                bad=True #Doesn't matter, too large
            else:
                print "This is tau:", tau
                print "This is dtau:", dtau
                raise Exception('tau integrator has large error.')
    else:
        if np.abs(dtau)>1e-3:
            if tau > 0.15:
                bad=True #Doesn't matter, too large
            else:
                print "This is tau:", tau
                print "This is dtau:", dtau
                raise Exception('tau integrator has large error.')

    return tau, Q_of_z, bad


def logprior(x):
    f_esc_params, c_hii, photon_norm_factor, ombh2, ommh2 = unpack(x)

    bad = False
    const = 1.0  #Unnormalized prior

    #Prior ranges
    if params.ion_model=="Standard":
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
            if params.f_esc_monotonic:
                if  any([f_sort != f_nosort for f_sort, f_nosort in zip(np.sort(f_esc_params),f_esc_params)]):
                #if np.sort(f_esc_params)!=f_esc_params:
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


    elif params.ion_model=="Nonparametric":

        if any(f_esc_params<0.0):
            bad = True

        elif np.mean(f_esc_params[-4]>1e-18):
            bad = True

    else:
        raise Exception('This ion model not implemented.')

    if c_hii<1.0 or c_hii>5.0:
        bad = True

    if photon_norm_factor <24.0 or photon_norm_factor>27.0:
        bad = True

    if bad:
        return -np.inf, False
    else:
        return const, True

def loglike(tau,Q,x):
    f_esc_params, c_hii, photon_norm_factor, ombh2, ommh2 = unpack(x)

    #Non-CMB constraints

    #2sig constraint from McGreer+, 1411.5375 as a step function
    try:
        #if Q(z=5.9)<0.84:
        if any(map(lambda z: Q(z)<0.84,np.linspace(0,5.9,10))):
            #print "Caught by Q"
            return -np.inf
    except:
        raise TypeError('Call to Q did not work.')

    #2\sigma constraint from Boutsia et al 2011
    if params.ion_model=="Standard":
        if params.use_lowfesc_const:
            #Initiate f_esc object
            f_esc = f_esc_funct(f_esc_flag, f_esc_params)
            if any(map(lambda z: f_esc.f_esc(z)>0.10,np.linspace(0.0,3.3,10))):
                #print "Caught by f_esc"
                return -np.inf

    if data_type=="tau_only":

        if tau<np.min(globe.data['tau']) or tau > np.max(globe.data['tau']):
            return -np.inf

        else:

            #print "This is post_funct", globe.post_funct(tau)
            post = globe.post_funct([tau])

            return post #Like is already log-like
            #return np.log(post)

    elif data_type == "marg_cosmo":

        if tau<np.min(globe.data['tau']) or tau > np.max(globe.data['tau']):
            #print "Caught by tau"
            return -np.inf
        elif ombh2<np.min(globe.data['omegabh2']) or ombh2 > np.max(globe.data['omegabh2']):
            #print "Caught by ombh2"
            return -np.inf
        elif ommh2<np.min(globe.data['omegamh2']) or ommh2 > np.max(globe.data['omegamh2']):
            #print "Caught by ommh2"
            return -np.inf
        else:

            post = globe.post_funct([ombh2,ommh2,tau])

            #print "this is post:", post

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

    tau, Q, bad = tau_calculator(x)

    if bad:
        like = -np.inf
    else:
        like = loglike(tau, Q, x)


        #print "this is tau", tau, Q(5.9)
        #print "this is loglike", like
        #print "this is logprior", prior
        #print "this is bad", bad
        if np.mod(globe.print_counter,100)==0:
            print "----------------------", globe.print_counter
            print "this is tau", tau, Q(5.9)
            print "this is loglike", like
            print "this is logprior", prior
            print "this is bad", bad

    globe.print_counter +=1

    return prior + like
