import matplotlib.pyplot as plt
import numpy as np
from parameters_channels import *

def plot_gsolution(ax, solver, label="", time_in_steps=False, 
                   g0_ratio=True, capacitance=False):
    
    data = 0

    if solver=='odeint':

        data = np.loadtxt(f"{DATA_PATH}solution_odeint.txt", unpack=True)

    if solver=='euler':
        
        data = np.loadtxt(f"{DATA_PATH}solution_euler.txt", unpack=True)

    if capacitance and solver=='odeint':

        data = np.loadtxt(f"{DATA_PATH}solution_odeint_g.txt", unpack=True)

    if capacitance and solver=='euler':

        data = np.loadtxt(f"{DATA_PATH}solution_euler_g.txt", unpack=True)


    time_interval = data[0]
    solution = data[1]

    if time_in_steps:

        time_interval = list(range(0,len(time_interval)))

    if g0_ratio == False:

        solution *= g_0*10**(12)
        ax.set_ylabel(r'$g(t)[pS]$', fontsize = axis_fontsize)
        
    else:

        ax.set_ylabel(r'$g(t)/g_0$', fontsize = axis_fontsize)

    ax.plot(time_interval, solution, **gsolution_style)
 
    ax.grid(ls=':')
    ax.set_xlabel(r't[ms]', fontsize = axis_fontsize)

    ax.tick_params('y', labelsize=size_ticks)
    ax.tick_params('x', labelsize=size_ticks)
    ax.set_xlim(np.min(time_interval), np.max(time_interval))
    
    # label
    ax.text(-0.1, 1, label, transform=ax.transAxes, fontsize=size_labels, va='top', ha='right')

def plot_voltage(ax, label="", yaxis_label=False, capacitor=False, odeint=False, euler=False):

    data = np.loadtxt(f"{DATA_PATH}voltage_file.txt", unpack=True)
   
    if capacitor and odeint:

        data = np.loadtxt(f"{DATA_PATH}solution_odeint_Vg.txt", unpack=True)

    if capacitor and euler:

        data = np.loadtxt(f"{DATA_PATH}solution_euler_Vg.txt", unpack=True)

    if capacitor==False:

        data = np.loadtxt(f"{DATA_PATH}voltage_file.txt", unpack=True)

    time_interval = data[0]
    solution = data[1]

    if capacitor:
        ax.plot(time_interval, solution, **voltage_capacitor_style)
    else:
        ax.plot(time_interval, solution, **voltage_triangle_style)
    # ax.plot(list(range(0,len(time_interval))), solution, **voltage_triangle_style)

    ax.grid(ls=':')
    ax.set_xlabel(r't[ms]', fontsize = axis_fontsize)

    ax.tick_params('y', labelsize=size_ticks)
    ax.tick_params('x', labelsize=size_ticks)
    
    ax.set_xlim(np.min(time_interval), np.max(time_interval))

    if yaxis_label:
        ax.set_ylabel(r'$V[V]$', fontsize = axis_fontsize)
        
    # label
    ax.text(-0.15, 1, label, transform=ax.transAxes, fontsize=size_labels, va='top', ha='right')
    
    print("time at + ", delta_t_sq, ", time at -", period_sq-delta_t_sq)

def plot_steady_solution(ax, yaxis_label=False):

    data = np.loadtxt(f"{DATA_PATH}steady_solution.txt", unpack=True)

    time_interval = data[0]
    solution = data[1]

    ax.plot(time_interval, solution, **steady_solution_style)

    ax.grid(ls=':')
    ax.set_xlabel(r't[ms]', fontsize = axis_fontsize)

    ax.tick_params('y', labelsize=size_ticks)
    ax.tick_params('x', labelsize=size_ticks)
    ax.set_xlim(np.min(time_interval), np.max(time_interval))

    if yaxis_label:
        ax.set_ylabel(r'$g_{\infty}[V(t)]/g_0$', fontsize = axis_fontsize)

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

def plot_conductance_voltage(ax, label=""):

    data_conductance = np.loadtxt(f"{DATA_PATH}solution_odeint.txt", unpack=True)
    data_voltage = np.loadtxt(f"{DATA_PATH}voltage_file.txt", unpack=True)

    conductance = data_conductance[1]
    voltage = data_voltage[1]
    
    conductance = conductance[400:]
    voltage = voltage[400:]
    
    print(len(voltage))

    ax.plot(voltage, conductance, **gsolution_style)

    # ax.grid(ls=':')
    ax.set_ylabel(r'$g(V)/g_0$', fontsize = axis_fontsize)
    ax.set_xlabel(r'V[V]', fontsize = axis_fontsize)

    ax.tick_params('y', labelsize=size_ticks)
    ax.tick_params('x', labelsize=size_ticks)
    
    ax.text(-0.2, 1, label, transform=ax.transAxes, fontsize=size_labels, va='top', ha='right')

def plot_current_voltage(ax, solution_type, label="", capacitance=False):
    
    data = np.loadtxt(f"{DATA_PATH}voltage_file.txt", unpack=True)
    voltage = np.array(data[1])
    
    if solution_type=='static':
        data = np.loadtxt(f"{DATA_PATH}steady_solution.txt", unpack=True)
    if solution_type=='dynamic':
        data = np.loadtxt(f"{DATA_PATH}solution_odeint.txt", unpack=True)
    if capacitance and solution_type=='euler':
        data = np.loadtxt(f"{DATA_PATH}solution_euler_g.txt", unpack=True)

    conductance = np.array(data[1])
    
    current = voltage*conductance*(g_0*10**(12))
    
    ax.plot(voltage, current, **current_style)
    
    ax.axhline(0, color='dimgray', lw=1)
    ax.axvline(0, color='dimgray', lw=1)
    
    ax.grid(ls=':')
    ax.set_ylabel(r'I$[pA]$', fontsize = axis_fontsize)
    ax.set_xlabel(r'V[V]', fontsize = axis_fontsize)

    ax.tick_params('y', labelsize=size_ticks)
    ax.tick_params('x', labelsize=size_ticks)
    
    # label
    ax.text(-0.1, 1, label, transform=ax.transAxes, fontsize=size_labels, va='top', ha='right')
    
def plot_voltage_gsolution(ax, solution_type, label=""):
    
    data = np.loadtxt(f"{DATA_PATH}voltage_file.txt", unpack=True)
   
    time_interval = data[0]
    solution = data[1]

    ax.plot(time_interval, solution, **voltage_triangle_style)

    ax.grid(ls=':')
    ax.set_xlabel(r't[ms]', fontsize = axis_fontsize)

    ax.tick_params('y', labelsize=size_ticks)
    ax.tick_params('x', labelsize=size_ticks)
    
    ax.set_xlim(np.min(time_interval), np.max(time_interval))

    ax.set_ylabel(r'$V[V]$', fontsize = axis_fontsize)
        
    # label
    ax.text(-0.2, 1, label, transform=ax.transAxes, fontsize=size_labels, va='top', ha='right')
        
    ax2 = ax.twinx()
    
    if solution_type=='dynamic':
        
        data = np.loadtxt(f"{DATA_PATH}solution_odeint.txt", unpack=True)
        solution = data[1]
        ax2.plot(time_interval, solution, **gsolution_style)
        
        ax2.set_ylabel(r'$g(t)/g_0$', fontsize = axis_fontsize)
        ax2.spines['right'].set_color('navy')
        ax2.tick_params(axis='y', colors='navy', labelsize=size_ticks)
        
    if solution_type=='static':
        
        data = np.loadtxt(f"{DATA_PATH}steady_solution.txt", unpack=True)
        solution = data[1]
        ax2.plot(time_interval, solution, **steady_solution_style)
        
        ax2.set_ylabel(r'$g_{\infty}/g_0$', fontsize = axis_fontsize)
        ax2.spines['right'].set_color('firebrick')
        ax2.tick_params(axis='y', colors='firebrick', labelsize=size_ticks)
    
    
    

    