import matplotlib.pyplot as plt
import numpy as np
from parameters import *

def plot_gsolution(ax, solver, which_code):
    
    data = 0

    if which_code=='my_code':
        if solver=='odeint':

            data = np.loadtxt(f"{DATA_PATH}solution_odeint.txt", unpack=True)
    
        if solver=='euler':
            
            data = np.loadtxt(f"{DATA_PATH}solution_euler.txt", unpack=True)
    if which_code=='tim':
            data = np.loadtxt(f"{DATA_PATH}solution_tim.txt", unpack=True)


    time_interval = data[0]
    solution = data[1]

    if which_code=='mycode':
        ax.plot(time_interval, solution, **gsolution_style)
    else:
        ax.plot(time_interval, solution, **gsolution_style_tim)


    ax.grid(ls=':')
    ax.set_ylabel(r'$g(t)/g_0$', fontsize = axis_fontsize)
    ax.set_xlabel(r't', fontsize = axis_fontsize)

    ax.tick_params('y', labelsize=size_ticks)
    ax.tick_params('x', labelsize=size_ticks)
    ax.set_xlim(np.min(time_interval), np.max(time_interval))

def plot_voltage(ax, which_code):

    if which_code=='my_code':
        data = np.loadtxt(f"{DATA_PATH}voltage_file.txt", unpack=True)
    if which_code=='tim':
        data = np.loadtxt(f"{DATA_PATH}voltage_file_tim.txt", unpack=True)
   
    time_interval = data[0]
    solution = data[1]

    if which_code=='my_code':
        ax.plot(time_interval, solution, **voltage_triangle_style)
    else:
        ax.plot(time_interval, solution, **voltage_triangle_style_tim)

    ax.grid(ls=':')
    ax.set_ylabel(r'$V(t)$', fontsize = axis_fontsize)
    ax.set_xlabel(r't', fontsize = axis_fontsize)

    ax.tick_params('y', labelsize=size_ticks)
    ax.tick_params('x', labelsize=size_ticks)
    ax.set_xlim(np.min(time_interval), np.max(time_interval))

def plot_avdensity(ax, which_code):

    if which_code=='my_code':
        data = np.loadtxt(f"{DATA_PATH}average_density_steady.txt", unpack=True)
    if which_code=='tim':
        data = np.loadtxt(f"{DATA_PATH}average_density_tim.txt", unpack=True)

    time_interval = data[0]
    solution = data[1]

    if which_code=='my_code':
        ax.plot(time_interval, solution, **averagedensity_style)
    else:
        ax.plot(time_interval, solution, **averagedensity_style_tim)


    ax.grid(ls=':')
    ax.set_ylabel(r'$V(t)$', fontsize = axis_fontsize)
    ax.set_xlabel(r't', fontsize = axis_fontsize)

    ax.tick_params('y', labelsize=size_ticks)
    ax.tick_params('x', labelsize=size_ticks)
    ax.set_xlim(np.min(time_interval), np.max(time_interval))

def plot_voltage_density_steady(ax):

    data_density = np.loadtxt(f"{DATA_PATH}average_density_steady.txt", unpack=True)
    data_voltage = np.loadtxt(f"{DATA_PATH}voltage_file.txt", unpack=True)

    time_interval = data_density[0]
    desnity = data_density[1]
    voltage = data_voltage[1]

    ax.plot(time_interval, voltage, **voltage_triangle_style)

    ax.grid(ls=':')
    ax.set_ylabel(r'$V(t)$', fontsize = axis_fontsize)
    ax.set_xlabel(r't', fontsize = axis_fontsize)

    ax.tick_params('y', labelsize=size_ticks)
    ax.tick_params('x', labelsize=size_ticks)
    ax.set_xlim(np.min(time_interval), np.max(time_interval))

    ax2 = ax.twinx()

    ax2.plot(time_interval, desnity, **averagedensity_style)

    ax2.set_ylabel(r'$\langle \rho_s \rangle / \left( 2\rho_b \right)$', fontsize = axis_fontsize)

def plot_conductance_voltage(ax, which_code):

    if which_code=='my_code':

        data_conductance = np.loadtxt(f"{DATA_PATH}solution_odeint.txt", unpack=True)
        data_voltage = np.loadtxt(f"{DATA_PATH}voltage_file.txt", unpack=True)

    if which_code=='tim':

        data_conductance = np.loadtxt(f"{DATA_PATH}solution_tim.txt", unpack=True)
        data_voltage = np.loadtxt(f"{DATA_PATH}voltage_file_tim.txt", unpack=True)

    conductance = data_conductance[1]
    voltage = data_voltage[1]

    if which_code=='my_code':
        ax.plot(voltage, conductance, **gsolution_style)
    else:
        ax.plot(voltage, conductance, **gsolution_style_tim)

    ax.grid(ls=':')
    ax.set_ylabel(r'$g(t)/g_0$', fontsize = axis_fontsize)
    ax.set_xlabel(r'V(t)', fontsize = axis_fontsize)

    ax.tick_params('y', labelsize=size_ticks)
    ax.tick_params('x', labelsize=size_ticks)
    # ax.set_xlim(np.min(time_interval), np.max(time_interval))

