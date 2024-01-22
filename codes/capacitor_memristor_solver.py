import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from parameters_channels import *
from scipy import signal
import scipy.integrate as integrate
import potential_shapes
import odeint_solver

def model(variables, t):

    g_V_t, V_g = variables

    potential = potential_shapes.sine_potential(t, amplitude_sine, freq_sine)

    potential_der = freq_sine*potential_shapes.cosine_potential(t, amplitude_sine, freq_sine)

    dgdt = g_0*10**(12)*(odeint_solver.g_infinity_func(V_g) - g_V_t) / tau

    dVgdt = -(1/capacitance) * g_V_t * V_g + potential_der

    return [dVgdt, dgdt]

def odeint_solver_function():

    time_interval = np.linspace(initial_time,final_time,time_steps)
    Vg = 0
    inital_cond = [0, (inital_condition*10**(-12))/g_0]

    print(inital_cond)

    print(capacitance)

    solution = odeint(model, inital_cond, time_interval)

    Vg, g = solution.T

    solution_odeint_Vg = open(f"{DATA_PATH}solution_odeint_Vg.txt", "w")
    solution_odeint_g = open(f"{DATA_PATH}solution_odeint_g.txt", "w")

    for time in range(len(time_interval)):

        solution_odeint_Vg.write(f'{time_interval[time]} \t {solution[time][0]} \n')
        solution_odeint_g.write(f'{time_interval[time]} \t {g[time]} \n')


def euler_forward_solver():

    solution_euler_g = open(f"{DATA_PATH}solution_euler_g.txt", "w")
    solution_euler_Vg = open(f"{DATA_PATH}solution_euler_Vg.txt", "w")

    g=(inital_condition*10**(-12))/g_0

    V_g = 0
    
    solution_euler_g.write(f'{initial_time} \t {g} \n')
    solution_euler_Vg.write(f'{initial_time} \t {V_g} \n')

    time_interval = np.linspace(initial_time + (final_time-initial_time)/time_steps,final_time,time_steps)

    for time in time_interval:
        
        potential_der = freq_sine*potential_shapes.cosine_potential(time, amplitude_sine, freq_sine)

        g_old = g
        Vg_old = V_g

        g += (odeint_solver.g_infinity_func(potential_der) - g_old)/(tau)*timestep_size
        
        V_g += ((1/capacitance)*(g_old*Vg_old) + potential_der)*timestep_size
        
        solution_euler_g.write(f'{time} \t {g} \n')
        solution_euler_Vg.write(f'{time} \t {V_g} \n')

        
