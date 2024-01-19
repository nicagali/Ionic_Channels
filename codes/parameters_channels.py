import numpy as np

# PATHS

DATA_PATH = '../data/'
PLOT_PATH = '../plots/'

PATH_TO_DIR_ION = "/Users/monicaconte/PhD/Projects/Ionic_Channels/"

# PARAMETERS ODE

inital_condition = 1
tau = 4.8
initial_time = 0
final_time = 70
# final_time = 5
time_steps = 1000
timestep_size = final_time/time_steps
dx = 1e-7

time_interval = np.linspace(initial_time,final_time,time_steps)

# PARAMETERS G INFINITY

delta_g = -3.6
length_channel = 10e-6
radius_base = 200e-9
radius_tip = 50e-9
peclet_number = 16.5
rho_b = 0.1
electron_charge = 1.60217663e-19
Dvalue = 1.75e-9
boltzman_const = 1.38e-23
temperature = 293.15

delta_radius = radius_base - radius_tip

g_1 = np.pi*radius_tip*radius_base/length_channel
g_2 = 2*rho_b*(electron_charge**2)*Dvalue/(boltzman_const*temperature)
g_0 = g_1*g_2

avogadro_number = 6.022e23

g_0 *= avogadro_number

inital_condition = 4

capacitance = 1e-14

# PARAMETERS SAWTOOTH POTENTIAL

# frequency = 40
period = 30
frequency_tr = 1/period

# PARAMETERS SQUARE POTENTIAL

# for increase bw -5 and 4
    # duty_sq = 0.2
    # frequency_sq = 0.1
    # period_sq=1/frequency_sq
    # amplitude_sq = -1/2
    # # phase_sq = np.pi*frequency_sq
    # phase_sq = 0
    # shift_sq = -4.5

    # delta_t_sq = period_sq*duty_sq

# for decrease bw 0 and 5

duty_sq = 0.1
frequency_sq = 0.1
period_sq=1/frequency_sq
amplitude_sq = -2
# phase_sq = np.pi*frequency_sq
phase_sq = 0
shift_sq = 2

delta_t_sq = period_sq*duty_sq

# PLOT PREFERENCES

# size
axis_fontsize = 17
legend_size = 15
size_ticks = 13
size_labels = 18 

# dictionaries

gsolution_style = dict(c = 'navy', lw=2)
voltage_triangle_style = dict(c='limegreen', lw=2, label=r'$V(t)$')
steady_solution_style = dict(c='firebrick', lw=2, label=r'$g_{\infty}(V)/g_0$')
current_style = dict(c='skyblue', lw=2)