import numpy as np
import matplotlib.pyplot as plt
import matplotlib._cntr as cntr

from scipy import interpolate
from math import isnan

def table_from_tab(filename):
    """
    Generic function to read a file containing floats and commented lines
    into a 2D numpy array.

    Commented lines are prefixed by the characters # or %.
    """
    f=open(filename, 'r')
    data = []
    datastream = f.read()
    f.close()
    datalines = [ line.strip().split() for line in datastream.split('\n') if line.strip() ]
    
    minp = float(datalines[4][0])
    stepp = float(datalines[5][0])
    nps = int(datalines[6][0])
    maxp = minp + (nps-1.)*stepp
    mint = float(datalines[8][0])
    stept = float(datalines[9][0])
    nts = int(datalines[10][0])
    maxt = mint + (nts-1.)*stept
    nproperties = int(datalines[11][0]) - 2
    properties = datalines[12][2:]

    pressures = np.linspace(minp*1.e5, maxp*1.e5, nps)
    temperatures = np.linspace(mint, maxt, nts)
    property_array = np.array([[[0. for p in pressures] for t in temperatures] for prop in properties])
    
    for c, prop in enumerate(properties):
        col = c+2
        n=13
        for i, temperature in enumerate(temperatures):
            for j, pressure in enumerate(pressures):
                data.append(map(float, line))
                value=float(datalines[n][col])
                if isnan(value):
                    value=0. # Crude replacement of nans for zeros to avoid interp throwing the toys out of the pram.
                property_array[c][i][j] = value # property, then T, then P
                n = n+1
                
    return pressures, temperatures, property_array

def interpolate_isentropic_T_rho_Vp_Vs(T_ref, filename):
    pressures, temperatures, property_array = table_from_tab(filename)
    s, rho, vp, vs = property_array
    entropy = interpolate.interp2d(pressures, temperatures, s, kind='linear')
    density = interpolate.interp2d(pressures, temperatures, rho, kind='linear')
    v_p = interpolate.interp2d(pressures, temperatures, vp, kind='linear')
    v_s = interpolate.interp2d(pressures, temperatures, vs, kind='linear')
    S_ref = entropy(pressures.min(), T_ref)
    if isnan(S_ref):
        print('Oh no! For some reason the entropy interpolation on your perplex file returns NaN at the minimum pressure.')
        print('Conditions: P='+str(pressures.min()/1.e5)+' bar, T='+str(T_ref)+' K.')
        exit()
    adiabatic_temperatures = np.empty_like(pressures)
    adiabatic_densities = np.empty_like(pressures)
    adiabatic_Vps = np.empty_like(pressures)
    adiabatic_Vss = np.empty_like(pressures)
    for i, s_p in enumerate(zip(*s)):
        entropy_slice = interpolate.interp1d(s_p, temperatures, kind='linear')
        adiabatic_temperatures[i] = entropy_slice(S_ref)
        adiabatic_densities[i] = density(pressures[i], adiabatic_temperatures[i])
        adiabatic_Vps[i] = v_p(pressures[i], adiabatic_temperatures[i])
        adiabatic_Vss[i] = v_s(pressures[i], adiabatic_temperatures[i])
    
    adiabat_density = interpolate.interp1d(pressures, adiabatic_densities, kind='linear')
    adiabat_Vp = interpolate.interp1d(pressures, adiabatic_Vps, kind='linear')
    adiabat_Vs = interpolate.interp1d(pressures, adiabatic_Vss, kind='linear')
    adiabat_temperature = interpolate.interp1d(pressures, adiabatic_temperatures, kind='linear')

    return adiabat_temperature, adiabat_density, adiabat_Vp, adiabat_Vs
