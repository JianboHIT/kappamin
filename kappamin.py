import configparser
from collections import OrderedDict
import numpy as np
from scipy.integrate import quad
import os
import sys
import io
import re


UNITS = OrderedDict()
UNITS['T'] = 'K'                    # temperature
UNITS['Cv'] = 'kB'                  # heat capacity in kB
UNITS['Kappa_min'] = 'W/(m*K)'      # minimum limit kappa under tau=pi/omega
UNITS['Kappa_min_A'] = 'W/(m*K)'    # contribution of acoustic branches
UNITS['Kappa_min_O'] = 'W/(m*K)'    # contribution of optical branches
UNITS['vs'] = 'km/s'                # average sound velocity, vs = (3/(2/vT**3+1/vL**3))**(1/3)
UNITS['Tau_min'] = 'ps'             # Kappa_min/(1/3*Cv*vs^2)
UNITS['MFP_min'] = 'A'              # Kappa_min/(1/3*Cv*vs)
UNITS['Omega_a_T'] = 'rad/ps'       # Cut-off angular frequency of TA
UNITS['Omega_a_L'] = 'rad/ps'       # Cut-off angular frequency of LA
UNITS['T_a_T'] = 'K'                # Debye temperature of TA
UNITS['T_a_L'] = 'K'                # Debye temperature of LA

def _kernel(x):
    '''
    x^2*exp(x)/(exp(x)-1)^2
    '''
    p = np.power(x, 2) * np.exp(-x)
    q = np.power(1-np.exp(-x), 2)
    v = np.divide(p, q, 
                  out=np.ones_like(p+q), 
                  where=(np.absolute(q)>1E-4))
    return v

def _core_cv(x):
    '''
    3*t^2
    '''
    v = 3 * x**2
    return v

def _core_debye(x):
    '''
    3*t
    '''
    v = 3 * x
    return v

def _core_bvk(x):
    '''
    3 * x^2 * (cos(pi/2 * x))^2 / sin(pi/2 * t)
    '''
    p = 3 * x**2 * np.cos(np.pi/2 * x)**2
    q = np.sin(np.pi/2 * x)
    v = np.divide(p, q,
                  out=np.zeros_like(p+q),
                  where=(np.absolute(q)>1E-4))
    return v

def _quad_t(func, Tr):
    '''
    Calculate integral in reduced k-space with an additional parameter
    '''
    def itg(Trx):
        return quad(func, 0, 1, args=(Trx,))[0]
    itg_v = np.vectorize(itg)
    return itg_v(Tr)

def Debye(vT, vL, Natom, Vcell, T):
    '''
    Calculate Debye-Cahill minimum limit to thermal conductivity.
    All  parameters in common units.
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
    vs = (3/(2/vT**3+1/vL**3))**(1/3)
    factor_cv = 3       # [kB]
    factor_kmT = (kB*vT*vT)/Vatom * np.pi/WcT  # [W/m.K]
    factor_kmL = (kB*vL*vL)/Vatom * np.pi/WcL
    # factor_tmT = np.pi/WcT      # [ps]
    # factor_tmL = np.pi/WcL
    # factor_mfp = np.pi/Kc       # [A]
    
    # calculate
    out = dict()
    if isinstance(T, float) and (T == float('inf')):
        # kernel --> 1
        # factor_itg = quad(_core_debye, 0, 1)[0]
        factor_itg = 3/2
        out['T'] = T
        out['Cv'] = factor_cv
        out['Kappa_min'] = (2*factor_kmT+factor_kmL)/3 * factor_itg
        out['Omega_a_T'] = WcT
        out['Omega_a_L'] = WcL
        out['T_a_T'] = TaT
        out['T_a_L'] = TaL
        out['vs'] = vs
    else:
        T = np.array(T)
        f_cv = lambda t, u: _core_cv(t)*_kernel(u*t)
        f_km = lambda t, u: _core_debye(t)*_kernel(u*t)
        CrT = _quad_t(f_cv, TaT/T)
        CrL = _quad_t(f_cv, TaL/T)
        KMrT = _quad_t(f_km, TaT/T)
        KMrL = _quad_t(f_km, TaL/T)
        out['T'] = T
        out['Cv'] = factor_cv * (2*CrT+CrL)/3
        out['Kappa_min'] = (2*factor_kmT*KMrT+factor_kmL*KMrL)/3
        out['Omega_a_T'] = WcT*np.ones_like(T)
        out['Omega_a_L'] = WcL*np.ones_like(T)
        out['T_a_T'] = TaT*np.ones_like(T)
        out['T_a_L'] = TaL*np.ones_like(T)
        out['vs'] = vs*np.ones_like(T)
        
    out['Kappa_min_A'] = out['Kappa_min']
    out['Kappa_min_O'] = 0 * out['Kappa_min']
    out['Tau_min'] = out['Kappa_min']/(1/3 * out['Cv']*kB/Vatom * out['vs']**2)     # [ps]
    out['MFP_min'] = out['Tau_min'] * out['vs'] * 10 # [ps*km/s] = [ps*nm/ps] = [nm] to [A]
    return out

def BvK(vT, vL, Natom, Vcell, T):
    '''
    Calculate BvK-Cahill minimum limit to thermal conductivity.
    All  parameters in common units.
    '''
    
    # define constants for matching common units
    kB = 13.80649
    hb = 105.457182
    
    # define some factors and integrals
    Vatom = Vcell/Natom     # [A^3]
    Kc = np.power(6*np.pi*np.pi/Vatom, 1/3)  # [1/A]
    WcT = 2/np.pi * vT * Kc * 10    # [km/s * 1/A] = [10 rad/ps] to [rad/ps]
    WcL = 2/np.pi * vL * Kc * 10
    TaT = hb/kB * WcT   # [K]
    TaL = hb/kB * WcL
    factor_cv = 3       # [kB]
    vs = (3/(2/vT**3+1/vL**3))**(1/3)
    factor_kmT = (kB*vT*vT)/Vatom * np.pi/WcT  # [W/m.K]
    factor_kmL = (kB*vL*vL)/Vatom * np.pi/WcL
    # factor_tmT = np.pi/WcT      # [ps]
    # factor_tmL = np.pi/WcL
    # factor_mfp = (np.pi/Kc) * np.pi/2   # [A]
    
    # calculate
    out = dict()
    if isinstance(T, float) and (T == float('inf')):
        # kernel --> 1
        # factor_itg = quad(_core_bvk, 0, 1)[0]
        factor_itg = 0.31456063126172384
        out['T'] = T
        out['Cv'] = factor_cv
        out['Kappa_min'] = (2*factor_kmT+factor_kmL)/3 * factor_itg
        out['Omega_a_T'] = WcT
        out['Omega_a_L'] = WcL
        out['T_a_T'] = TaT
        out['T_a_L'] = TaL
        out['vs'] = vs
    else:
        T = np.array(T)
        f_cv = lambda t, u: _core_cv(t)*_kernel(u*np.sin(np.pi/2 * t))
        f_km = lambda t, u: _core_bvk(t)*_kernel(u*np.sin(np.pi/2 * t))
        CrT = _quad_t(f_cv, TaT/T)
        CrL = _quad_t(f_cv, TaL/T)
        KMrT = _quad_t(f_km, TaT/T)
        KMrL = _quad_t(f_km, TaL/T)
        out['T'] = T
        out['Cv'] = factor_cv * (2*CrT+CrL)/3
        out['Kappa_min'] = (2*factor_kmT*KMrT+factor_kmL*KMrL)/3
        out['Omega_a_T'] = WcT*np.ones_like(T)
        out['Omega_a_L'] = WcL*np.ones_like(T)
        out['T_a_T'] = TaT*np.ones_like(T)
        out['T_a_L'] = TaL*np.ones_like(T)
        out['vs'] = vs*np.ones_like(T)
    
    out['Kappa_min_A'] = out['Kappa_min']
    out['Kappa_min_O'] = 0 * out['Kappa_min']
    out['Tau_min'] = out['Kappa_min']/(1/3 * out['Cv']*kB/Vatom * out['vs']**2)     # [ps]
    out['MFP_min'] = out['Tau_min'] * out['vs'] * 10 # [ps*km/s] = [ps*nm/ps] = [nm] to [A]
    return out

def Pei(vT, vL, Natom, Vcell, T):
    '''
    Calculate Pei-Cahill minimum limit to thermal conductivity.
    All  parameters in common units.
    '''
    
    vs = (3/(2/vT**3+1/vL**3))**(1/3)
    out = BvK(vT=vs,
              vL=vs,
              Natom=1,
              Vcell=Vcell,
              T=T)
    
    # same as BvK model when Natom = 1
    if Natom == 1:
        return out
    # else:
    #     raise NotImplementedError('Only support Natom = 1')
    
    # for Natom != 1
    # define constants for matching common units
    kB = 13.80649
    hb = 105.457182
    
    # define some factors and integrals
    Vatom = Vcell/Natom     # [A^3]
    Kc = np.power(6*np.pi*np.pi/Vatom, 1/3)  # [1/A], Kc ~ Vatom, Kcc ~ Vcell
    Wc = 2/np.pi * vs * Kc * 10     # [km/s * 1/A] = [10 rad/ps] to [rad/ps]
    To = hb/kB * Wc     # [K]
    factor_cv = 3       # [kB]
    factor_km = (kB*vs*vs)/Vcell * np.pi/Wc  # [W/m.K], here Vcell=Natom*Natom
    # factor_tm = np.pi/Wc      # [ps]
    # factor_mfp = (np.pi/Kc) * np.pi/2   # [A]
    
    # calculate
    kir = [(np.power(i/Natom, 1/3)+np.power((i+1)/Natom, 1/3))/2 for i in range(1,round(Natom))]
    if isinstance(T, float) and (T == float('inf')):
        kir = np.array(kir)
        wr = np.sin(np.pi/2 * kir)
        vr = np.cos(np.pi/2 * kir)
        km_O = factor_km * vr**2 / wr
        out['Kappa_min_O'] = np.sum(km_O)
        out['Kappa_min'] = out['Kappa_min_A'] + out['Kappa_min_O']
    else:
        T = out['T']
        shp = [-1,]+[1 for _ in range(T.ndim)]
        kir = np.array(kir).reshape(shp)
        wr = np.sin(np.pi/2 * kir)
        vr = np.cos(np.pi/2 * kir)
        km_O = factor_km * vr**2 / wr * _kernel(To/T*wr)
        cv_O = factor_cv * _kernel(To/T*wr)
        out['Kappa_min_O'] = np.sum(km_O, axis=0)
        out['Kappa_min'] = out['Kappa_min_A'] + out['Kappa_min_O']
        out['Cv'] = (out['Cv']+np.sum(cv_O, axis=0))/Natom
    
    out['Tau_min'] = out['Kappa_min']/(1/3 * out['Cv']*kB/Vatom * out['vs']**2)     # [ps]
    out['MFP_min'] = out['Tau_min'] * out['vs'] * 10 # [ps*km/s] = [ps*nm/ps] = [nm] to [A]
    return out

def fileparser(filename, ktypes=None):
    '''
    Parser the configuration file and return a ConfigParser() object.
    '''
    # read config file
    config = configparser.ConfigParser(inline_comment_prefixes=('#',))
    config.SECTCRE = re.compile(r"\[ *(?P<header>[^]]+?) *\]")  # ignore blanks in section name
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
    
def _section_parser(config, modeltype):
    '''
    Parser parameters from a specific model-section.
    '''
    
    # TODO: config = fileparser(config) if is string else config
    # keys = ['vT', 'vL', 'Natom', 'Vcell', 'T']
    
    paras = dict()
    
    # parser float
    for key in ['vT', 'vL', 'Natom', 'Vcell']:
        paras[key] = config.getfloat(modeltype, key)
    
    # parser temperature    
    valueT = config[modeltype]['T']
    isStep, valueT = _sequence_parser(valueT, func=float, defaultstep=1)
    if isStep:
        start, step, end = valueT
        valueT = np.arange(start, end+0.01, step)
    
    # check type and return
    if len(valueT) == 0:
        raise ValueError('Failed to read temperature [T].')
    elif len(valueT) == 1:
        isSingle = True
        paras['T'] = valueT[0]
    else:
        isSingle = False
        paras['T'] = valueT
    return isSingle, paras

def _sequence_parser(value, func=float, defaultstep=1):
    '''
    Parser a sequence or slice from a string.
    '''
    
    if ':' in value:
        # sequence with step
        values = value.strip().split(':')
        if func is not None:
            values = list(map(func, values))
        if len(values) == 2:
            start, end = values
            step = defaultstep
        elif len(values) == 3:
            start, step, end = values
        else:
            raise ValueError('Failed to parser the sequence data.')
        isStep = True
        return isStep, (start, step, end)
    else:
        # any sequence
        values = value.strip().split()
        if func is not None:
            values = map(func, values)
        isStep = False
        return isStep, list(values)

def _savedat_to_file(filename, datas, keys=None, 
                     fmt='%.6f', header='auto', 
                     isSingle=None):
    '''
    Save output of model calculation to a file.
    Default header(='auto') is Properties & Units, while is ignored if it is None. 
    '''
    
    # TODO: fix isSingle check bug if datas['Kappa_min'] not a ndarray
    
    # check default
    if keys is None:
        keys = UNITS.keys()
    if isSingle is None:
        if datas['Kappa_min'].ndim == 0:
            isSingle = True
        else:
            isSingle = False
    if header.lower().startswith('auto'):
        props = []
        units = []
        for prop in keys:
            unit = UNITS[prop]
            Nmax = max(map(len, [prop, unit]))
            props.append(prop.rjust(Nmax))
            units.append(unit.rjust(Nmax))
        propsline = ', '.join(props)
        unitsline = ', '.join(units)
        header = '\n'.join([propsline, unitsline])
    
    # access output datas
    out = [datas[key] for key in keys]
    if isSingle:
        out = np.atleast_2d(out)
    else:
        out = np.vstack(out).T
    s = io.StringIO()
    np.savetxt(s, out, fmt=fmt, header=header)
    rst = s.getvalue().replace('\ninf', '\ninf       ')
    with open(filename, 'w') as f:
        f.write(rst)

def execute(filename=None, toFile=True, hasReturn=False):
    '''
    Read configuration file and do calculation.
    '''
    
    # TODO: allow muilt-output result when muilt-model
    
    # access filename of config
    ktypes = ['Debye', 'BvK', 'Pei']
    if filename is None:
        # auto-detect filename
        prefix = 'KAPPAMIN'
        suffix = ['cfg', 'ini', 'txt']
        for fn in ['{}.{}'.format(prefix, suf) for suf in suffix]:
            if os.path.exists(fn):
                filename = fn
                break
        else:
            raise FileNotFoundError('Failed to find configuration file')
    config = fileparser(filename, ktypes)
    
    # parser config file and calculate
    if 'Debye' in config.sections():
        fileout = 'Kappamin_Debye.dat'
        isSingle, paras = _section_parser(config, 'Debye')        
        out = Debye(**paras)
    elif 'BvK' in config.sections():
        fileout = 'Kappamin_BvK.dat'
        isSingle, paras = _section_parser(config, 'BvK')        
        out = BvK(**paras)
    elif 'Pei' in config.sections():
        fileout = 'Kappamin_Pei.dat'
        isSingle, paras = _section_parser(config, 'Pei')
        out = Pei(**paras)
    else:
        raise ValueError('Unknown method.(Valid: %s)', ', '.join(ktypes))
    
    # output and(or) save to file
    if hasReturn:
        return out
    if toFile:
        _savedat_to_file(fileout, out, isSingle=isSingle)
    
    
if __name__ == '__main__':
    if len(sys.argv) > 1:
        filename = sys.argv[1]
    else: 
        filename = None
    execute(filename)

