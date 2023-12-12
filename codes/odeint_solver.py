import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from parameters_channels import *
from scipy import signal
import scipy.integrate as integrate
import potential_shapes

# Integrand in the definition of average concentration

def integrand(x, potential):

    integrand=0

    if np.abs(peclet_number*potential)>=0.1:

        radius_x = radius_base - (x*delta_radius)/length_channel

        integrand1 = (x*radius_tip)/(radius_x*(length_channel))

        integrand2num = np.exp(peclet_number*potential*(x/length_channel)*((radius_tip)**2/(radius_base*radius_x))) - 1
        
        integrand2den = (np.exp(peclet_number*potential*radius_tip/radius_base) - 1)

        integrand = integrand1 - integrand2num/integrand2den

    return integrand

# Compute integral and give the value of g/g_0 = \rho_s (=average concentration)

def g_infinity_func(potential): 

    integral_ginfty = integrate.quad(integrand, 0, length_channel, args=(potential,), points=length_channel/dx)[0]/length_channel

    g_infty = 1 + delta_g*integral_ginfty

    return g_infty

# Define the model in the differential equation dg/dt = model(g,t,potential)

def model(g, t, potential_shape):

    potential = potential_shapes.potential(t, shape=f'{potential_shape}')

    g_infinity = g_infinity_func(potential)

    dgdt = (1/tau)*(g_infinity - g)

    return dgdt

# Copute the steady solution: perform integral in equation (5) 
# for a specific potential V(t)

def steady_solution(potential_shape):

    time_interval = np.linspace(initial_time,final_time,time_steps)
    steady_sol_vec = np.zeros((time_steps))

    steady_sol_file = open(f"{DATA_PATH}steady_solution.txt", "w")

    for time in range(time_steps):

        steady_sol_vec[time] = g_infinity_func(
            potential_shapes.potential(time_interval[time], shape=f'{potential_shape}'))
        
        steady_sol_file.write(f'{time_interval[time]} \t {steady_sol_vec[time]} \n')

# Solve differential equation in (7) for a specific potential shape V(t)

def odeint_solver_function(potential_shape):

    time_interval = np.linspace(initial_time,final_time,time_steps)
    g=0
    g = odeint(model, inital_condition, time_interval, args=(potential_shape, ))

    print(inital_condition)

    solution_odeint = open(f"{DATA_PATH}solution_odeint.txt", "w")

    for time in range(len(time_interval)):
        solution_odeint.write(f'{time_interval[time]} \t {g[time][0]} \n')

def euler_forward_solver(potential_shape):

    time_interval = np.linspace(initial_time,final_time,time_steps)

    solution_euler = open(f"{DATA_PATH}solution_euler.txt", "w")

    g=inital_condition
    
    for time in time_interval:
        
        potential = potential_shapes.potential(time, shape=f'{potential_shape}')

        g += (g_infinity_func(potential) - g)/(tau)*timestep_size
        
        solution_euler.write(f'{time} \t {g} \n')
        
