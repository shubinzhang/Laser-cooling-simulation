# -*- coding: utf-8 -*-
setting = {
    "dE": 0.024,                #energy difference between excitation and emission  (eV)
    "Eg": 1.417,                 #emission peak energy (eV)
    "tn": 2*10**8,               #time step uplimit
    "tolerance_charge": 0.005,   #charge distribution converge criteria
    "tolerance_T": 10**(-7), #temperature distribution converge criteria
    "mpnt": 100,                 #calculoate tolerance every m time steps
    "Tcool": 300,                #boundary temperature (K)

    "L": 9*10**(-4),             #cooling material radius (cm)
    "xn": 512,                   #number of finite element in radius
    "thickness":10**(-5),        #thickness(cm)

    "eps": 13., 	                #GaAs relative dielectric constant
    "mue": 8500.,                #electron mobility cm2/(V·s) 
    "muh": 400.,		            #hole mobility cm^2/(V*s)
    "Dce": 200,                  #electron diffusivity cm^2/s
    "Dch":10,		            #hole diffusivity cm^2/s 
    "Dt": 0.31,		            #heat diffusivity
    "M": 144.645,                #molar mass g/mol
    "rho": 5.32, 	            #volume density g/cm^3
    "Cp": 47.02,		            #heat capacity J/(mol*K)

    "kr": 4*10**(-10),           #bimocular recombination rate constant   cm^3/s
    "extract_eta": 0.024,        #photon extraction efficiency
    "kt": 5*10**(3),             #(non-radiative) electron trapping rate constant /s  From Mansor's paper
    "ka": 4*10**(-30),           #Auger recombination rate constant  cm^6/s
    

    "sigma": 0.5*10**(-4),       #laser spot size (cm)
    "laser_power": 10**(-3),     #laser power (W)
    "showprogress": False,       #show progress
    "save": False               #save data



}