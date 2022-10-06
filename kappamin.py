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
    kB = 1.380649
    hb = 105.457182
    
    # calculate
    Vatom = Vcell/Natom     # in A^3
    if isinstance(T, float) and (T == float('inf')):
        # kernel --> 1
        itg = 1/2
        vitgs = (2*vT+vL) * itg    # in km/s
    else:
        fx = lambda x, u: x * kernel(x*u)        # u = Td/T
        kc = np.power(6*np.pi*np.pi/Vatom, 1/3)  # in 1/A
        TdT = hb/kB * kc * vT
        TdL = hb/kB * kc * vL
        itgT = quad(fx,0,1,args=(TdT/T,))[0]
        itgL = quad(fx,0,1,args=(TdL/T,))[0]
        vitgs = 2*vT*itgT + vL*itgL
    factor = np.power(np.pi/6, 1/3) * kB
    n_23 = np.power(1/Vatom, 2/3)
    return factor * n_23 * vitgs

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
    