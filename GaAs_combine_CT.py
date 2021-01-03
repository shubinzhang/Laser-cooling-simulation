'''
Python script to measure the emission quantum yields with differrnt excitation wavelength.

Original code by  Shubin Zhang 2017 
Adapted by Yurii Morozov 2018
Adapted by Shubin Zhang 2020

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation
'''

import sys
import os
import copy
import atexit

import time
from time import localtime, strftime

import numpy as np
from scipy import integrate
from scipy.constants import e, epsilon_0
from scipy.integrate import odeint

import matplotlib
import matplotlib.pyplot as plt

from config import setting
#from PyQt4 import QtCore, QtGui

class Small_GaAs():
	def __init__(self):
		#load material parameters from config.py, more detail can be found in config.py

		self.dE = setting["dE"]
		self.Eg = setting["Eg"]
		self.tn = setting["tn"]
		self.tolerance_charge = setting["tolerance_charge"]
		self.tolerance_T = setting["tolerance_T"]
		self.mpnt = setting["mpnt"]
		self.Tcool = setting["Tcool"]
		self.L = setting["L"]	
		self.T = setting["thickness"]
		self.xn =  setting["xn"]
		eps = setting["eps"]
		self.mue = setting["mue"]
		self.muh = setting["muh"]
		self.Dce = setting["Dce"]
		self.Dch = setting["Dch"]
		Dt = setting["Dt"]  						
		self.kr = setting["kr"] * setting["extract_eta"]   
		self.kt = setting["kt"]
		self.ka = setting["ka"]
		self.sigma = setting["sigma"] 
		laser_power = setting["laser_power"]
		self.M = setting["M"]
		self.rho = setting["rho"]
		self.Cp = setting["Cp"]

		self.dx = self.L/self.xn   								#finite element size (cm)
		self.xv = np.linspace(0.,self.L-self.dx, self.xn)   	#radius array 
		self.dx2 = self.dx**2	

		xv00 = self.xv[1:-1]
		xim = self.xv - self.dx
		xip = self.xv + self.dx
		self.eoveps = e/(eps*epsilon_0)	
		self.gampvalues =  laser_power/(self.Eg*e) 		#charge carrier generation rate 
		self.cdt =  1.*self.dx2/(4*self.Dce)   			#charge calculation time step


		self.De_plus =  self.Dce*(xv00+self.dx/2)/(xv00*self.dx2)
		self.De_minus = self.Dce*(xv00-self.dx/2)/(xv00*self.dx2)
		self.Dh_plus =  self.Dch*(xv00+self.dx/2)/(xv00*self.dx2)
		self.Dh_minus = self.Dch*(xv00-self.dx/2)/(xv00*self.dx2)

		self.mue_plus =  self.mue*(xv00+self.dx)/(2*xv00*self.dx)
		self.mue_minus = self.mue*(xv00-self.dx)/(2*xv00*self.dx)
		self.muh_plus =  self.muh*(xv00+self.dx)/(2*xv00*self.dx)
		self.muh_minus = self.muh*(xv00-self.dx)/(2*xv00*self.dx)

		self.last_factor = (self.xn - 1.5)/(self.xn - 1)
		self.last_factor_mu = (self.xn - 2)/(self.xn - 1)

		self.Ascld = self.kt
		self.Bscld = self.kr/self.T
		self.Cscld = 0.5*self.ka/self.T**2

		self.Tdt =  1*self.dx2/(4*Dt) 					#Temperature calculation time step
		self.Dt1coef = Dt*self.Tdt/self.dx2   			
		self.Dt2coef = Dt*self.Tdt/(2*self.dx*self.xv[1:-1])
		self.D0coef = 2*Dt*(self.Tdt/self.dx2) 

		self.showprogress = setting["showprogress"]
		self.save = setting["save"]

		if self.showprogress:
			app = QtGui.QApplication(sys.argv)
			line_plot = pg.plot()
			line_plot2.showGrid(x=True, y=True, alpha=1.)
			line_plot2 = pg.plot()
			line_plot2.showGrid(x=True, y=True, alpha=1.)
			line_plot3 = pg.plot()
			line_plot3.showGrid(x=True, y=True, alpha=1.)
			line_plot4 = pg.plot()
			line_plot4.showGrid(x=True, y=True, alpha=1.)
			line_plot5 = pg.plot()
			line_plot5.showGrid(x=True, y=True, alpha=1.)
			curv0 = line_plot.plot(pen=None, symbolPen=None, symbolSize=3, symbolBrush=(255, 255, 255))#
			curv1 = line_plot.plot(pen=None, symbolPen=None, symbolSize=3, symbolBrush=(255, 0, 0))#
			curv2 = line_plot.plot(pen=None, symbolPen=None, symbolSize=3, symbolBrush=(0, 255, 0))#
			curv3 = line_plot2.plot(pen=None, symbolPen=None, symbolSize=3, symbolBrush=(255, 255, 0))#
			curv4 = line_plot3.plot(pen=None, symbolPen=None, symbolSize=3, symbolBrush=(255, 255, 100))#
			curv5 = line_plot4.plot(pen=None, symbolPen=None, symbolSize=3, symbolBrush=(255, 255, 200))#
			curv6 = line_plot4.plot(pen=None, symbolPen=None, symbolSize=3, symbolBrush=(0, 255, 200))#
			curv7 = line_plot4.plot(pen=None, symbolPen=None, symbolSize=3, symbolBrush=(0, 100, 250))#
			curv8 = line_plot5.plot(pen=None, symbolPen=None, symbolSize=2, symbolBrush=(255, 100, 250))#

		if self.save:
			filename = strftime("%d_%b_%H", localtime())
			if not os.path.exists(filename):
				os.mkdir(filename)
				os.mkdir(filename+"/norcharges")
				os.mkdir(filename+"/nortotalcharge")
				os.mkdir(filename+"/norelectricfield")
				os.mkdir(filename+"/nortemperature")
				os.mkdir(filename+"/nordt")
				os.mkdir(filename+"/norcentertemperature")
				os.mkdir(filename+"/norheatsource")
				os.mkdir(filename+"/norpotential")





	def Generation(self):
		"""
		Excitation intensity follows Gaussian profile in polar coordianate
		"""
		normc = 1./(self.sigma**2*(2.*np.pi)*self.T)
		self.G = self.gampvalues*normc*np.exp(-self.xv**2/(2*self.sigma**2))   # /(cm^3*s)


	def poissolv(self, psi, ro,  eps0):
		"""
		calculate spatial electric potential iteratively
		psi: previous spatical electric potential
		ro: spatical charge distribtuion
		eps0: tolerance factor for converge criteria
		"""
		eps = 1+eps0
		ind = 0
		while eps>eps0:
			newpsi = np.zeros(self.xn)
			newpsi[1:-1] =  (psi[2:]+psi[:-2])/2 + (1./(4*self.xv[1:-1]))*(psi[2:]-psi[:-2])*self.dx  + 0.5*self.dx2*ro[1:-1]
			newpsi[0] = newpsi[1] + 0.5*self.dx2*ro[0]
			newpsi[-1] = 0
			eps = max(abs(newpsi - psi))
			psi = newpsi.copy()
			ind+=1
		return psi, ind
		
	def do_timestep_charge(self, nl, pl, psil):
		"""
		calculate new spatial charge distribution after one time step 
		"""
		n = np.zeros(len(nl))
		p = np.zeros(len(pl))
		n00 = nl[1:-1]
		p00 = pl[1:-1]	
		cind = 0
		ronet = (pl-nl)*self.eoveps/self.T
		eps0 = self.dx2*max(ronet)/100.
		psi = psil
		Ef = np.zeros(self.xn)
		psi, niter_pois, = self.poissolv(psil, ronet, eps0)
		Ef[1:-1] = -(psi[2:]-psi[:-2])/(2*self.dx)
		Ef[-1] = sum(ronet*self.xv) *self.dx/self.xv[-1]
		Ef[0] = 0
		if max(abs(Ef))!= 0:
			dt =  min(1.*self.dx2/(4*self.Dce), 1.*self.dx/(10*self.mue*max(abs(Ef))))
		GR_common = self.G*self.T - self.Bscld*nl*pl - self.Cscld*(nl**2*pl + nl*pl**2)
		GRe = GR_common - self.Ascld*nl # Generation rate for electrons
		GRh = GR_common - self.Ascld*pl # Generation rate for holes
		Ef00 = Ef[1:-1]
		Ef_diff = (Ef[2:] - Ef[:-2])/(2*self.dx)
		n00_rate  = nl[2:]*(self.De_plus + self.mue_plus*Ef[2:]) + nl[:-2]*(self.De_minus - self.mue_minus*Ef[:-2]) + n00*(-2*self.Dce/self.dx2)
		p00_rate  = pl[2:]*(self.Dh_plus - self.muh_plus*Ef[2:]) + pl[:-2]*(self.Dh_minus + self.muh_minus*Ef[:-2]) + p00*(-2*self.Dch/self.dx2)		
		n[1:-1] = n00 + self.cdt*n00_rate + GRe[1:-1]*self.cdt
		p[1:-1] = p00 + self.cdt*p00_rate + GRh[1:-1]*self.cdt
		
		#---------------------------------------------------------------------
		n[0] = nl[0] + 2*(self.Dce/self.dx2)*(nl[1]-nl[0])*self.cdt + GRe[0]*self.cdt + nl[1]*self.mue/self.dx*Ef[1]*self.cdt
		p[0] = pl[0] + 2*(self.Dch/self.dx2)*(pl[1]-pl[0])*self.cdt + GRh[0]*self.cdt - pl[1]*self.muh/self.dx*Ef[1]*self.cdt
		#Boundary condition for charge density
		n[-1] = nl[-1] + GRe[-1]*self.cdt  - self.last_factor_mu*self.mue*nl[-2]*Ef[-2]*self.cdt/(2*self.dx) - self.Dce*self.cdt/self.dx2*(nl[-1] - nl[-2])*self.last_factor
		p[-1] = pl[-1] + GRh[-1]*self.cdt  + self.last_factor_mu*self.muh*pl[-2]*Ef[-2]*self.cdt/(2*self.dx) - self.Dch*self.cdt/self.dx2*(pl[-1] - pl[-2])*self.last_factor 
		return n, p, psi ,Ef, n00_rate, niter_pois


	def total_recombination(self, n, p):
		"""
		calculate total electron/hole recombination rate due to trapped state recombination, bimocular recombination and Auger recombination
		"""
		qe = n/self.T
		qh = p/self.T
		ple = self.kr*qe*qh #1/(cm^3*S)
		plh = ple
		nrad_e = self.kt*qe #1/(cm^3*S)
		nrad_h = self.kt*qh #1/(cm^3*S)
		augere = self.ka*qe**2*qh   #1/(cm^3*S)
		augerh = self.ka*qh**2*qe
		return (sum(ple*self.xv) + sum(plh*self.xv) + sum(augere*self.xv) + sum(nrad_e*self.xv)+sum(nrad_h*self.xv)+ sum(augerh*self.xv))*2*np.pi*self.dx

	def solvediffusion_charge(self):
		"""
		calculate spatial charge distrbution iteratively until it converges
		"""
		stime = time.time()
		Ef = np.zeros(self.xn)
		n = np.zeros(self.xn)
		p = np.zeros(self.xn)
		psi = np.zeros(self.xn)


		ntotal = np.zeros(self.tn/self.mpnt+1)
		tol2 = 1
		
		m = 0
		while abs(tol2) > self.tolerance_charge: 
			n, p, psi, Ef,n00_rate, niter_points = self.do_timestep_charge(n, p, psi)
			if m>self.tn:
				print "end by tn, tn is equal to" ,self.tn, "m is equal to",m
				break
				
			mshort = m/self.mpnt
			if m%self.mpnt==0:
				ntotal[mshort] = np.sum(n*self.xv)
				if  m>self.mpnt*2:
					recombination = self.total_recombination(n,p)*self.T
					tol2 = 1-recombination/(2*self.gampvalues)

			if m%1000==0 and m > 1:
				print 'Current charge tolerance = %1.11f'%tol2
				print "current time step", self.cdt
				print "Number of poissolv iterations = ", niter_points
				print "Current iteration = ", m
				print "----------------------------------------------------------"
				net_c = p-n
				if self.showprogress:
						curv0.setData(self.G*max(n)/(max(self.G)))
						curv1.setData(n)
						curv2.setData(p)
						curv3.setData(Ef)
						curv4.setData(ntotal[:mshort-1])
						curv8.setData(net_c*self.xv)
						curv5.setData(n00_rate)
						genv = self.G*self.T - self.Ascld*(n+p)*0.5 - self.Bscld*n*p - self.Cscld*(n**2*p + n*p**2)
						curv6.setData(genv[1:-1])
						curv7.setData(n00_rate+genv[1:-1])
						QtCore.QCoreApplication.processEvents()
			m+=1
		print "real time",m*self.cdt
		print "Charge equilibrium reached ... "  
		print "charge tolerance ", tol2
		print "***************************************\n"
		print "Time taken = %1.5f"%(time.time() - stime), ' s'
		print "Number of iterations = ", m+1, ' out of ', self.tn, '. Or ', (m+1.)/self.tn*100.,' %'
		ntotal = ntotal[:m/self.mpnt-1]
		print "Final time = ", self.cdt*m*10**9, ' ns'
		if self.save:
			plt.plot(self.xv,qe,'b', linewidth=2)
			plt.plot(self.xv,qh,'r', linewidth=2)
			plt.savefig(filename+"/norcharges/"+"_%03d_"%gind +"charge_distribution" , dpi=400)
			plt.close()
			plt.plot(self.xv, Ef, linewidth=2)
			plt.savefig(filename+"/norelectricfield/"+"_%03d_"%gind + "electric_field" , dpi=400)
			plt.close()
			charges = zip(n,p)
			np.savetxt(filename+"/norcharges/"+"_%03d_"%gind + '_grate=%1.2e'%power[gind]+'.txt', charges)
			np.savetxt(filename+"/norelectricfield/"+"_%03d_"%gind + '_grate=%1.2e'%power[gind]+'.txt', Ef)
			np.savetxt(filename+"/norpotential/"+"_%03d_"%gind + '_grate=%1.2e'%power[gind]+'.txt', psi)
		return n,p

	def heating_cooling_power(self, n, p):
		"""
		calculate heat/cooling power base on stable spatial charge distribution
		"""
		qe = n/self.T  
		qh = p/self.T
		#recombination rate
		pl = self.kr*qe*qh #1/(cm^3*S)
		nrad_e = 0.5*self.kt*qe #1/(cm^3*S)
		nrad_h = 0.5*self.kt*qh #1/(cm^3*S)
		auger = 0.5*(self.ka*qe**2*qh + self.ka*qh**2*qe) #1/(cm^3*S)
		plqyvals = np.sum(pl*self.xv)/self.total_recombination(n, p)
		Auheat = auger*self.Eg*e# W/cm^3
		nradheat = (nrad_e+nrad_h)*self.Eg*e# W/cm^3
		Heating = Auheat + nradheat # W/cm^3
		Cooling = self.G*self.dE*e  # W/cm^3
		self.hbalancedensity = Heating - Cooling
		hbminvals = self.hbalancedensity[0]
		if self.save:
			plt.plot(self.xv,Heating,'r', linewidth=2)
			plt.plot(self.xv,Cooling,'b', linewidth=2)
			plt.plot(self.xv,hbalancedensity,'g', linewidth=2)
			plt.savefig(filename+"/norheatsource/"+"_%03d_"%gind +"charge_distribution" , dpi=400)
			plt.close()
			heating_cooling = zip(Heating,Cooling)
			np.savetxt(filename+"/norheatsource/"+"_%03d_"%gind + '_grate=%1.2e'%power[gind]+'.txt', self.hbalancedensity)
			np.savetxt(filename+'/plqyvals'+"_%03d_"%gind+'.txt', plqyvals)


	def do_timestep_temperature(self,ul):
		"""
		calculate new temperature distribution after one time step 
		"""
		u = np.zeros(len(ul))
		u[0] = ul[0] + self.D0coef*(ul[1]-ul[0]) + self.hbalancedensity[0]*self.M*self.Tdt/(self.rho*self.Cp)
		u[1:-1] = ul[1:-1] + self.Dt1coef*(ul[2:] - 2*ul[1:-1]  + ul[:-2]) + self.Dt2coef*(ul[2:] - ul[:-2]) + self.hbalancedensity[1:-1]*self.M*self.Tdt/(self.rho*self.Cp)
		u[-1]  =  self.Tcool
		return  u

	def solvediffusion_temperature(self):
		"""
		calculate temperature distrbution iteratively until it converges
		"""
		stime = time.time()
		u = self.Tcool*np.ones(self.xn)
		utotal = np.zeros(self.tn/self.mpnt+1)
		tol = 1
		m = 0
		while tol > self.tolerance_T:
			u = self.do_timestep_temperature(u)
			
			if m>self.tn:
				print "end by tn, tn is equal to" ,self.tn, "m is equal to",m
				break
				
			mshort = m/self.mpnt
			if m%self.mpnt==0:
				utotal[mshort] = u[0]
				if  m>self.mpnt*2:
					tol = abs(utotal[mshort] - utotal[mshort-1])/abs(300-utotal[mshort])
			if m%10000==0 and m > 1:
				print 'Current temperature tolerance = %1.11f'%tol
			m+=1
		print "Temperature equilibrium reached ... "  
		print "*************************************** \n"
		print "Time taken = %1.5f"%(time.time() - stime), ' s'
		utotal = utotal[:m/self.mpnt]
		if self.save:
			plt.plot(self.xv, u, linewidth=2)
			plt.savefig(filename+"/nortemperature/"+"_%03d_"%gind + "temperature" , dpi=400)
			plt.close() 
			plt.plot(utotal, linewidth=2)
			plt.savefig(filename+"/norcentertemperature/"+"_%03d_"%gind + "center_temperature", dpi=400)
			plt.close()
			np.savetxt(filename+"/nortemperature/"+"_%03d_"%gind + '_grate=%1.2e'%power[gind]+'.txt', u)

   


if __name__ == "__main__":
	simulation = Small_GaAs()
	gen_values = simulation.Generation()
	n, p = simulation.solvediffusion_charge()
	simulation.heating_cooling_power(n, p)
	simulation.solvediffusion_temperature()
