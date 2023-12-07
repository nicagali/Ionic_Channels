import numpy as np

# PATHS

DATA_PATH = '../data/'
PLOT_PATH = '../plots/'

# PARAMETERS ODE

inital_condition = 1.4
tau = 4.8e-3
initial_time = 0
final_time = 0.07
time_steps = 1000
timestep_size = final_time/1000
# initial_time = timestep_size
# final_time+=timestep_size
dx = 1e-7

# stepsize = final_time/time_steps
# initial_time = stepsize
# final_time = final_time + stepsize

# PARAMETERS G INFINITY

delta_g = -3.6
length_channel = 10e-6
radius_base = 200e-9
radius_tip = 50e-9
peclet_number = 16.5
rho_b = 0.1e-3
electron_charge = 1.60217663e-19
Dvalue = 1.75e-9
boltzman_const = 1.38e-23
temperature = 293.15

delta_radius = radius_base - radius_tip

g_1 = np.pi*radius_tip*radius_base/length_channel
g_2 = 2*rho_b*(electron_charge**2)*Dvalue/(boltzman_const*temperature)
g_0 = g_1*g_2

# PARAMETERS SAWTOOTH POTENTIAL

frequency = 40

# PLOT PREFERENCES

# size
axis_fontsize = 17
legend_size = 15
size_ticks = 13

# dictionaries

gsolution_style = dict(c = 'navy', lw=2)
voltage_triangle_style = dict(c='limegreen', lw=2, label=r'$V(t)$')
steady_solution_style = dict(c='firebrick', lw=2, label=r'$g_{\infty}/g_0(t)$')
