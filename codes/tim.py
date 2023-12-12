# -*- coding: utf-8 -*-
"""
Created on Mon Jan  9 16:35:44 2023

@author: Tim Veenstra
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad
from scipy import signal
from parameters_channels import *

# Physical parameters

kT = 1.38e-23 * 300  # thermal energy [J (or equivalently kg m^2 s^-2)] 
e = 1.602e-19  # elementary charge [C]
N_A = 6.022e23  # Avogadro's number [1]
C = 5e-15  # capacitance [F]
D = 1.75e-9  # diffusion constant of both ions species [m^2 / s]
sigma = -0.0015*1e18  # [1 / m^2]
rho_b = 0.1*N_A  # [1 / m^3]
w = -9.5

# Simulation parameters

frequency = 40
dt = 7e-5  # time integration step [s]
dx = 1e-7  # spatial integration step [m]
Nsteps = 1000

# Initial configuration
t = 0

def wave_potential(t):
    return np.cos(2 * np.pi*t*frequency) 

def sawtooth_potential(t):
    return signal.sawtooth(2*np.pi*frequency*t-np.pi/2, width=0.5)

def square_potential(t):
    return signal.square(2*np.pi*t*frequency, duty=0.1)


class Memristor:
    def __init__(self, init_conductivity, L, Rb, dR):
        self.conductivity = init_conductivity
        self.L = L  # length
        self.tau = L**2/(12*D)

        self.old_conductivity = 0

        self.Rb = Rb  # base radius of cone
        self.dR = dR  # difference in radius between base and tip
        self.Rt = Rb- dR  # tip radius
        self.Du = sigma/(2*rho_b*self.Rt)  # Duhkin number
        self.dg = -2 * w * dR/Rb * self.Du
        # print(self.dg)
        self.g0 = np.pi * self.Rt*Rb/L * 2*rho_b*e**2*D / kT
        

    def radius(self, x):
        # radius of the cone at x
        return self.Rb - x * self.dR

    def conductivity_function(self, x, Pe):
        # local conductivity at x for Peclet number Pe
                
        return x/self.L *self.Rt/(self.Rb-x/self.L * self.dR) \
            - (np.exp(Pe*x/self.L*self.Rt**2/((self.Rb-x/self.L*self.dR)*self.Rb)) - 1)\
            / (np.exp(Pe * self.Rt/self.Rb) - 1)

    def conductivity_limit(self, V):
        # steady state limit of the conductivity at voltage V
        Pe = V * 16.5 / 4 * self.Rb/self.Rt  # see p2 of Kamsma et al. 2023
        integral = 0
        if np.abs(Pe) > 0.1:
            integral = quad(self.conductivity_function, 0, self.L,args=(Pe,), points=int(self.L/dx))[0]/self.L
            # print(integral)
        return (1 + self.dg * integral)

    def update_conductivity(self, voltage_t):
        # Simple forward Euler scheme 
        ginf = self.conductivity_limit(voltage_t)
        self.old_conductivity = self.conductivity
        self.conductivity += ((self.conductivity_limit(voltage_t)
                              - self.conductivity)/self.tau*dt)


def tim_solver():

    chonk = Memristor(1.4, 10e-6, 200e-9, 150e-9)

    t_list = np.arange(0,Nsteps)*dt
    V_list = np.array([])
    g_list = np.array([])
    ginfinity_list = np.array([])

    voltage_file = open(f"{DATA_PATH}voltage_file_tim.txt", "w")
    gininfinity_file = open(f"{DATA_PATH}average_density_tim.txt", "w")
    g_file = open(f"{DATA_PATH}solution_tim.txt", "w")
    file = open(f"{DATA_PATH}checking_file_tim.txt", "w")
    
    output = -1

    for t in t_list:
        V = sawtooth_potential(t)
        
        V_list = np.append(V_list, V)
        g_list = np.append(g_list, chonk.conductivity)
        # ginfinity_list = np.append(ginfinity_list, chonk.conductivity_limit)

        
        chonk.update_conductivity(V)
        if t>0.75/frequency and output<0:
            output = 1 if chonk.conductivity > 1 else 0

        voltage_file.write(f'{t} \t {V} \n')
        gininfinity_file.write(f'{t} \t {chonk.conductivity_limit(V)} \n')
        g_file.write(f'{t} \t {chonk.conductivity} \n')

        file.write(f'{t} \t {V} \t {chonk.conductivity_limit(V)} \t {chonk.old_conductivity} \t {chonk.tau} \t {dt} \n')


    
chonk = Memristor(1.4, 10e-6, 200e-9, 150e-9)  
print(chonk.g0)  



