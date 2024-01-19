import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from parameters_channels import *
from scipy import signal
import scipy.integrate as integrate
import potential_shapes
import odeint_solver

def model(g, t, potential_shape):

    potential = potential_shapes.potential(t, shape=f'{potential_shape}')

    g_V_t = odeint_solver.g_infinity_func(potential)

    dVgdt = -(1/capacitance)*

    return dgdt

def odeint_solver_function(potential_shape):

    time_interval = np.linspace(initial_time,final_time,time_steps)
    Vg = 0
    inital_cond = 0
    Vg = odeint(model, inital_cond, time_interval, args=(potential_shape, ))

    solution_odeint = open(f"{DATA_PATH}solution_odeint.txt", "w")

    for time in range(len(time_interval)):
        
        print(g[time][0]*g_0)
        solution_odeint.write(f'{time_interval[time]} \t {g[time][0]} \n')