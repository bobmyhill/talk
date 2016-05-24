
import numpy as np
from scipy.integrate import odeint
from scipy.integrate import quad
from scipy.interpolate import interp1d
from scipy.interpolate import UnivariateSpline

from process_perplex_output import *


#gravitational constant
G = 6.67408e-11

class Planet(object):
    def __init__(self, n_slices, moho_depth, cmb_depth, planet_radius, core_density_deficit, core_Vp_fraction, T_pot, perplex_mantle_file):
        #The constructor takes the number of depth slices which will
        #be used for the calculation.  More slices will generate more
        #accurate profiles, but it will take longer.
        
        #The planet will be represented by a three layer model, crust, mantle and core.
        self.mantle_potential_temperature = T_pot
        self.moho = planet_radius - moho_depth
        self.cmb = planet_radius - cmb_depth  # Guess for the radius of the core-mantle-boundary
        self.outer_radius = planet_radius # Outer radius of the planet
        self.core_density_deficit = core_density_deficit
        self.core_Vp_fraction = core_Vp_fraction
            
        self.radii = np.linspace(0.e3, self.outer_radius, n_slices) # Radius list
        self.pressures = np.linspace(40.0e9, 0.0, n_slices) # initial guess at pressure profile
        self.temperatures = np.ones_like(self.pressures)*1000. #Assume an isothermal interior (not a great assumption, but we do this for simplicity's sake).

        core_pressures, core_temperatures, core_rhos, core_Vps = (np.loadtxt('data/core/isentropic_PTrhoVp.dat')).T
        self.core_density_function = interp1d(core_pressures*1.e9, core_rhos, kind='linear')
        self.core_Vp_function = interp1d(core_pressures*1.e9, core_Vps/1000., kind='linear')
            
        self.mantle_temperature_function, self.mantle_density_function, self.mantle_Vp_function, self.mantle_Vs_function = interpolate_isentropic_T_rho_Vp_Vs(self.mantle_potential_temperature, perplex_mantle_file)
    
    def generate_profiles( self, n_iterations ):
        #Generate the density, gravity, pressure, vphi, and vs profiles for the planet.
        #The n_iterations parameter sets how many times to iterate over density, pressure, 
        #and gravity.  Empirically, five-ish iterations is adequate.  After the iterations,
        #this also calculates mass and moment of inertia of the planet.

        for i in range(n_iterations):
            self.densities, self.V_p, self.V_s = self._populate_model(self.pressures, self.radii)
            self.gravity = self._compute_gravity(self.densities, self.radii)
            self.pressures = self._compute_pressure(self.densities, self.gravity, self.radii)

        self.mass = self._compute_mass(self.densities, self.radii)
        self.moment_of_inertia = self._compute_moment_of_inertia(self.densities, self.radii)
        self.moment_of_inertia_factor = self.moment_of_inertia / self.mass / self.outer_radius / self.outer_radius

        pfn = interp1d(self.radii, self.pressures, kind='linear')
        self.cmb_pressure = pfn(self.cmb)
        self.cmb_temperature = self.mantle_temperature_function(self.cmb_pressure)
            
    def _populate_model(self, pressures, radii):
        #Evaluates the equation of state for each radius slice of the model.
        #Returns density, V_p, V_s
        rho = np.empty_like(radii)
        Vp = np.empty_like(radii) 
        Vs = np.empty_like(radii) 
    
        for i, r in enumerate(radii):
            if r < self.cmb:
                rho[i] = self.core_density_function(pressures[i]) - self.core_density_deficit
                Vp[i] = self.core_Vp_function(pressures[i])*self.core_Vp_fraction
                Vs[i] = 0.00
            elif r < self.moho:
                rho[i] = self.mantle_density_function(pressures[i])
                Vp[i] = self.mantle_Vp_function(pressures[i])
                Vs[i] = self.mantle_Vs_function(pressures[i])
            else:
                rho[i]=3000.
                Vp[i] = 6.5
                Vs[i] = 3.75
        pfunc = interpolate.interp1d(radii, pressures)
        self.cmb_pressure = pfunc(self.cmb)
        self.moho_pressure = pfunc(self.moho)
            
        ps = [self.cmb_pressure, self.moho_pressure]
        rhos=self.mantle_density_function(ps)
        Vps=self.mantle_Vp_function(ps)
        Vss=self.mantle_Vs_function(ps)
                
        self.cmb_rhoVpVs = [[self.core_density_function(self.cmb_pressure) - self.core_density_deficit,
                             self.core_Vp_function(self.cmb_pressure)*self.core_Vp_fraction, 0.00],
                             [rhos[0], Vps[0], Vss[0]]]
        self.moho_rhoVpVs = [[rhos[1], Vps[1], Vss[1]], [3000., 6.5, 3.75]]
                
                
        return rho, Vp, Vs


    def _compute_gravity(self, density, radii):
        #Calculate the gravity of the planet, based on a density profile.  This integrates
        #Poisson's equation in radius, under the assumption that the planet is laterally
        #homogeneous. 
     
        #Create a spline fit of density as a function of radius
        rhofunc = UnivariateSpline(radii, density )

        #Numerically integrate Poisson's equation
        poisson = lambda p, x : 4.0 * np.pi * G * rhofunc(x) * x * x
        grav = np.ravel(odeint( poisson, 0.0, radii ))
        grav[1:] = grav[1:]/radii[1:]/radii[1:]
        grav[0] = 0.0 #Set it to zero a the center, since radius = 0 there we cannot divide by r^2
        return grav

    def _compute_pressure(self, density, gravity, radii):
        #Calculate the pressure profile based on density and gravity.  This integrates
        #the equation for hydrostatic equilibrium  P = rho g z.

        #convert radii to depths
        depth = radii[-1]-radii
        
        #Make a spline fit of density as a function of depth
        rhofunc = UnivariateSpline( depth[::-1], density[::-1] )
        #Make a spline fit of gravity as a function of depth
        gfunc = UnivariateSpline( depth[::-1], gravity[::-1] )

        #integrate the hydrostatic equation
        pressure = np.ravel(odeint( (lambda p, x : gfunc(x)* rhofunc(x)), 0.0,depth[::-1]))
        return pressure[::-1]

    def _compute_mass( self, density, radii):
        #Returns a list of moments of inertia of the planet [kg m^2]
        rhofunc = UnivariateSpline(radii, density )
        mass = quad( lambda r : 4*np.pi*rhofunc(r)*r*r, 
                                 radii[0], radii[-1] )[0]
        return mass
  
   
    def _compute_moment_of_inertia( self, density, radii):
        #Returns the moment of inertia of the planet [kg m^2]

        rhofunc = UnivariateSpline(radii, density )
        moment = quad( lambda r : 8.0/3.0*np.pi*rhofunc(r)*r*r*r*r,
                       radii[0], radii[-1] )[0]
        return moment
