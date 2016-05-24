# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for the Earth and Planetary Sciences
# Copyright (C) 2012 - 2015 by the BurnMan team, released under the GNU GPL v2 or later.



'''
This python code has been scrounged from example_build_planet, which was written for
the burnman code by Ian Rose.
'''


from __future__ import absolute_import
from __future__ import print_function

import os, sys

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.mlab import griddata

from scipy.integrate import odeint
from scipy.integrate import quad
from scipy.interpolate import interp1d
from scipy.interpolate import UnivariateSpline

from scipy.optimize import minimize

from process_perplex_output import *

from planet import Planet

if __name__ == "__main__":

    
    infile='data/mantle/grid_compositions/NCFMAS_grid_compositions_fine_2.txt'
    outfile='output_properties_fine_grid_2_1773.txt'
    T_pot = 1773.15
    moho_depth = 100.e3
    
    #gravitational constant
    G = 6.67408e-11
    
    GM_Mars_system = 4.282837453e13
    GM_Phobos = 7.17e5
    GM_Deimos = 9.43e4
    GM_Mars_system_uncertainty = 5.3e5
    
    M_Mars = (GM_Mars_system - GM_Phobos - GM_Deimos)/G
    # These are the actual observables
    observed_mass = M_Mars # 6.4171e23 # kg (Konopliv et al., 2011)
    mass_uncertainty = GM_Mars_system_uncertainty/G
    observed_mean_radius = 3389.5e3 # km (Seidelmann et al., 2007)
    observed_moment = 0.3635  # +/- 0.0012 (Sohl et al., 2005)

    core_Vp_fraction = 1.0 # only mineral physics constraints, set at 1.0 for now.

    def misfit(args, n_slices, n_iterations, moho_depth, T_pot, perplex_mantle_file):
        cmb_depth, core_density_deficit = args
        core_density_deficit = core_density_deficit / 1000.
        if core_density_deficit < -3000. or core_density_deficit > 5000. or cmb_depth < 100.e3:
            misfit = 1.e12
        else:
            mars = Planet(n_slices, moho_depth, cmb_depth, observed_mean_radius, core_density_deficit, core_Vp_fraction, T_pot, perplex_mantle_file)
            mars.generate_profiles( n_iterations )
            misfit = np.sqrt((1. - mars.mass/observed_mass)*(1. - mars.mass/observed_mass)
                            + (1. - mars.moment_of_inertia_factor/observed_moment)*(1. - mars.moment_of_inertia_factor/observed_moment))
        print(args, misfit)
        return misfit
    
    # Here we actually do the iteration.  We make an instance
    # of our Planet, then call generate_profiles.
    # Emprically, 300 slices and 5 iterations seem to do
    # a good job of converging on the correct profiles.
    n_slices = 1000
    n_iterations = 10

    perplex_mantle_files = []
    NCFMAS_compositions=np.loadtxt(fname=infile)
    for composition in NCFMAS_compositions:
        cstring='{0:.1f}_{1:.1f}_{2:.1f}_{3:.1f}_{4:.1f}_{5:.1f}'.format(composition[0], composition[1], composition[2], composition[3], composition[4], composition[5]).replace(".", ",")
        filename='data/mantle/grid_compositions/mars_mantle_'+cstring+'_1.tab'
        if os.path.isfile(filename):
            perplex_mantle_files.append(filename)

    NCFMAS_compositions = []
    cmb_depths = np.zeros(len(perplex_mantle_files))
    core_density_deficits = np.zeros(len(perplex_mantle_files))
    cmb_pressures = np.zeros(len(perplex_mantle_files))
    cmb_temperatures = np.zeros(len(perplex_mantle_files))
    misfits = np.zeros(len(perplex_mantle_files))
    Mg_numbers = np.zeros(len(perplex_mantle_files))
    Si_contents = np.zeros(len(perplex_mantle_files))
    output = []
    for i, perplex_mantle_file in enumerate(perplex_mantle_files):
        NCFMAS_compositions.append(map(float, perplex_mantle_file.replace(',', '.').split("_")[3:9]))
        
        mFeO = 71.844
        mMgO = 40.3044
        Mg_numbers[i] = 100.*(NCFMAS_compositions[i][3]/mMgO) / ((NCFMAS_compositions[i][2]/mFeO) + (NCFMAS_compositions[i][3]/mMgO))
        Si_contents[i] = NCFMAS_compositions[i][5]
        
        print('Processing composition', str(i+1)+'/'+str(len(perplex_mantle_files))+':', NCFMAS_compositions[i], Mg_numbers[i])
    
        
        res = minimize(misfit, [1700.e3, 2000.e3], method='Nelder-Mead', tol=1.e4,
                       args=(n_slices, n_iterations, moho_depth, T_pot, perplex_mantle_file)) # tol is depth in m
                      
        cmb_depths[i] = res.x[0]
        core_density_deficits[i] = res.x[1]/1000.
        mars = Planet(n_slices, moho_depth, cmb_depths[i], observed_mean_radius, core_density_deficits[i], core_Vp_fraction, T_pot, perplex_mantle_file)
        mars.generate_profiles( n_iterations )
        cmb_pressures[i] = mars.cmb_pressure
        cmb_temperatures[i] = mars.cmb_temperature
        misfits[i] = res.fun

        data_to_print = NCFMAS_compositions[i]
        data_to_print.extend([Mg_numbers[i], cmb_depths[i]/1000., cmb_pressures[i]/1.e9, cmb_temperatures[i], core_density_deficits[i], misfits[i]])

        print(data_to_print)
        output.append(data_to_print)
        
    
    np.savetxt(fname=outfile, header='NCFMAS, Mg number, CMB depth, CMB pressure, CMB temperature, Core density deficit, misfit', fmt="%.4f", X=output)

    '''
    extent = (min(Mg_numbers),max(Mg_numbers),min(Si_contents),max(Si_contents))
    N = 30
    xs,ys = np.mgrid[extent[0]:extent[1]:N, extent[2]:extent[3]:N]

    resampled = griddata(Mg_numbers, Si_contents, cmb_pressures/1.e9, xs, ys)

    plt.imshow(resampled.T, extent=extent)
    plt.plot(xs0, ys0, "r.")
    plt.plot(xs, ys, "b.")
    plt.title("CMB pressures")
    plt.show()
    
    '''
