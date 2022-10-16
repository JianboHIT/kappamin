#!/usr/bin/env python3

from kappamin import Debye, BvK, Pei

# Parameters
fileout = 'KAPPAMIN_models.dat'
paras = {
    'T':  float('inf'), 
    'vT': 4.37, 
    'vL': 7.36,
    'Natom': 2,
    'Vcell': 160.1,
}

# Calculate
out_deb = Debye(**paras)
out_bvk = BvK(**paras)
out_pei = Pei(**paras)

# Output
hdsp = '# {:<7s}{:<17s}{:<17s}\n'
ddsp = '{:<9s}{:<17.6f}{:<17.6f}\n'
with open(fileout, 'w') as f:
    f.write(hdsp.format('model', 'kL_min[W/(m.K)]', 'MFP[A]'))
    f.write(ddsp.format('Debye', out_deb['Kappa_min'], out_deb['MFP_min']))
    f.write(ddsp.format('BvK', out_bvk['Kappa_min'], out_bvk['MFP_min']))
    f.write(ddsp.format('Pei', out_pei['Kappa_min'], out_pei['MFP_min']))
