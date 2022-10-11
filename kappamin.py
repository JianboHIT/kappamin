import configparser
import numpy as np
from scipy.integrate import quad
import os
import sys
import re


def kernel(x):
    '''
    x^2*exp(x)/(exp(x)-1)^2
    '''
    p = np.power(x, 2) * np.exp(-x)
    q = np.power(1-np.exp(-x), 2)
    v = np.divide(p, q, 
                  out=np.ones_like(p+q), 
                  where=(np.absolute(q)>1E-4))
    return v

def Cahill(vT, vL, Natom, Vcell, T):
    '''
    Calculate Cahill-Debye minimum limit to thermal conductivity
    '''
    
    # define constants for matching common units
    kB = 13.80649
    hb = 105.457182
    
    # define some factors and integrals
    Vatom = Vcell/Natom     # [A^3]
    Kc = np.power(6*np.pi*np.pi/Vatom, 1/3)  # [1/A]
    WcT = vT * Kc * 10       # [km/s * 1/A] = [10 rad/ps] to [rad/ps]
    WcL = vL * Kc * 10
    TaT = hb/kB * WcT   # [K]
    TaL = hb/kB * WcL
    factor_cv = 3       # [kB]
    factor_kmT = (kB*vT*vT)/Vatom * np.pi/WcT  # [W/m.K]
    factor_kmL = (kB*vL*vL)/Vatom * np.pi/WcL
    factor_tmT = np.pi/WcT      # [ps]
    factor_tmL = np.pi/WcL
    
    f_cv = lambda t, u: 3*t*t*kernel(u*t)
    f_km = lambda t, u: 3*t*kernel(u*t)
    
    # calculate
    out = dict()
    if isinstance(T, float) and (T == float('inf')):
        # kernel --> 1
        out['T'] = T
        out['Cv'] = factor_cv
        out['Kappa_min'] = (2*factor_kmT+factor_kmL)/3 * 3/2
        out['Tau_min'] = (2*factor_tmT+factor_tmL)/3 * 3/2
        out['Omega_a_T'] = WcT
        out['Omega_a_L'] = WcL
        out['T_a_T'] = TaT
        out['T_a_L'] = TaL
    else:
        T = np.array(T)
        CrT = quad(f_cv, 0, 1, args=(TaT/T,))[0]
        CrL = quad(f_cv, 0, 1, args=(TaL/T,))[0]
        KMrT = quad(f_km, 0, 1, args=(TaT/T,))[0]
        KMrL = quad(f_km, 0, 1, args=(TaL/T,))[0]
        TMrT = KMrT/CrT
        TMrL = KMrL/CrL
        out['T'] = T
        out['Cv'] = factor_cv * (2*CrT+CrL)/3
        out['Kappa_min'] = (2*factor_kmT*TMrT+factor_kmL*TMrL)/3
        out['Tau_min'] = (2*factor_tmT*KMrT+factor_tmL*KMrL)/(2*CrT+CrL)
        out['Omega_a_T'] = WcT*np.ones_like(T)
        out['Omega_a_L'] = WcL*np.ones_like(T)
        out['T_a_T'] = TaT*np.ones_like(T)
        out['T_a_L'] = TaL*np.ones_like(T)
    return out

def BvK(vt, vl, natom, vcell, T):
    Vatom = vcell/natom
    v = (2*vt+vl)/3
    if T == float('inf'):
        pass
    else:
        pass

def BvK_Pei(vt, vl, natom, vcell, T):
    raise NotImplementedError

def fileparser(filename, ktypes=None):
    # read config file
    config = configparser.ConfigParser()
    config.SECTCRE = re.compile(r"\[ *(?P<header>[^]]+?) *\]")
    config.read(filename)
    
    # optional: refine config file
    if ktypes is None:
        return config
    else:
        config2 = configparser.ConfigParser()
        for sect in config.sections():
            for ktype in ktypes:
                if sect.lower() == ktype.lower():
                    config2[ktype] = config[sect]
        return config2


if __name__ == '__main__':
    # access filename of config
    ktypes = ['Cahill', 'BvK', 'BvK-Pei']
    if len(sys.argv) > 1:
        filename = sys.argv[1]
    else:
        prefix = 'KAPPAMIN'
        suffix = ['cfg', 'ini', 'txt']
        for fn in ['{}.{}'.format(prefix, suf) for suf in suffix]:
            if os.path.exists(fn):
                filename = fn
                break
        else:
            raise FileNotFoundError('Failed to find configuration file')
    config = fileparser(filename, ktypes)
    
    # parser config file
    if 'Cahill' in config.sections():
        keys = ['vT', 'vL', 'Natom', 'Vcell', 'T']
        paras = dict()
        for key in keys:
            paras[key] = config.getfloat('Cahill', key)
        Kmin = Cahill(**paras)
        if Kmin.ndim == 0:
            Kmin = np.array([Kmin])
        fileout = 'Kappamin_Cahill.dat'
        np.savetxt(fileout, Kmin, fmt='%.6f')
    elif 'BvK' in config.sections():
        keys = ['vT', 'vL', 'Natom', 'Vcell', 'T']
        paras = dict()
        for key in keys:
            paras[key] = config.getfloat('BvK', key)
        Kmin = BvK(**paras)
        if Kmin.ndim == 0:
            Kmin = np.array([Kmin])
        fileout = 'Kappamin_BvK.dat'
        np.savetxt(fileout, Kmin, fmt='%.6f')
    elif 'BvK-Pei' in config.sections():
        keys = ['vT', 'vL', 'Natom', 'Vcell', 'T']
        paras = dict()
        for key in keys:
            paras[key] = config.getfloat('BvK-Pei', key)
        Kmin = BvK_Pei(**paras)
        if Kmin.ndim == 0:
            Kmin = np.array([Kmin])
        fileout = 'Kappamin_BvK-Pei.dat'
        np.savetxt(fileout, Kmin, fmt='%.6f')
    else:
        raise ValueError('Unknown method.(Valid: %s)', ', '.join(ktypes))
    