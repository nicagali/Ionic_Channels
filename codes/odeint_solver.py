import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from parameters import *
from scipy import signal
import scipy.integrate as integrate

# Generate the desired potential shape

def triangular_potential(time, write=False):

    if write:

        t = np.linspace(initial_time,final_time,time_steps)
        
        triangle_waveform = signal.sawtooth(2*np.pi*frequency*t-np.pi/2, 0.5)

        voltage_file = open(f"{DATA_PATH}voltage_file.txt", "w")
        for time_index in range(len(t)):
            voltage_file.write(f'{t[time_index]} \t {triangle_waveform[time_index]} \n')


    return signal.sawtooth(2*np.pi*frequency*time-np.pi/2, 0.5)

# Integrand in the definition of average concentration

def integrand(x, potential):

    integrand=0

    if np.abs(peclet_number*potential)>=0.1:

        radius_x = radius_base - (x*delta_radius)/length_channel

        integrand1 = (x*radius_tip)/(radius_x*(length_channel))

        integrand2num = np.exp(peclet_number*potential*(x/length_channel)*((radius_tip)**2/(radius_base*radius_x))) - 1
        
        integrand2den = (np.exp(peclet_number*potential*radius_tip/radius_base) - 1)

        integrand = integrand1 - integrand2num/integrand2den

    # print(x, radius_tip, radius_x, length_channel, radius_base)

    # print(x, potential)

    # print(x, integrand1, integrand2num, integrand2den, integrand2num/integrand2den, integrand1 - integrand2num/integrand2den)

    # print(x, integrand2den, peclet_number, potential, radius_tip, radius_base)

    return integrand

# Compute integral and give the value of g/g_0 = \rho_s (=average concentration)

def g_infinity_func(potential): 

        
    integral_ginfty = integrate.quad(integrand, 0, length_channel, args=(potential,), points=length_channel/dx)[0]/length_channel
    
    # print(length_channel, integral_ginfty)

    g_infty = 1 + delta_g*integral_ginfty

    return g_infty

# Define the model in the differential equation dg/dt = model(g,t,potential)

def model(g, t, potential_shape):

    if potential_shape == 'triangular':

        potential = triangular_potential(t)

    g_infinity = g_infinity_func(potential)

    dgdt = (1/tau)*(g_infinity - g)

    return dgdt

# Copute the average salt concentration: perform integral in equation (5) 
# for a specific potential V(t)

def average_density_steady():

    time_interval = np.linspace(initial_time,final_time,time_steps)
    average_density = np.zeros((time_steps))

    avdensity_file = open(f"{DATA_PATH}average_density_steady.txt", "w")

    for time in range(time_steps):

        average_density[time] = g_infinity_func(
            triangular_potential(time_interval[time]))
        
        avdensity_file.write(f'{time_interval[time]} \t {average_density[time]} \n')

# Solve differential equation in (7) for a specific potential shape V(t)

def odeint_solver_function(potential_shape):

    time_interval = np.linspace(initial_time,final_time,time_steps)
    g=0
    g = odeint(model, inital_condition, time_interval, args=(potential_shape, ))
    # print(g)
    # g = odeint(trial_model, inital_condition, time_interval)

    solution_odeint = open(f"{DATA_PATH}solution_odeint.txt", "w")

    for time in range(len(time_interval)):
        solution_odeint.write(f'{time_interval[time]} \t {g[time][0]} \n')

def euler_forward_solver(potential_shape):
    
    file = open(f"{DATA_PATH}checking_file.txt", "w")
    
    time_interval = np.linspace(initial_time,final_time,time_steps)

    solution_euler = open(f"{DATA_PATH}solution_euler.txt", "w")

    g=inital_condition
    
    for time in time_interval:
        
        if potential_shape == 'triangular':
        
            potential = triangular_potential(time)
        
        # g += timestep_size*model(g, time, potential_shape)

        g += (g_infinity_func(potential) - g)/(tau)*timestep_size

        file.write(f'{time} \t {potential} \t {g_infinity_func(potential)} \t {g} \t {tau} \t {timestep_size} \n')
        
        solution_euler.write(f'{time} \t {g} \n')
        
