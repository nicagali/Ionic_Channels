from scipy import signal
import numpy as np
from parameters_channels import *

def triangular_potential(time):
    return signal.sawtooth(2*np.pi*frequency*time-np.pi/2, 0.5)

def square_potential(time, parameters):
    return 4*((parameters[0])*(signal.square(2*np.pi*frequency*time, 0.1) +parameters[1])/parameters[2] +5)

def potential(time, shape, which_conductance, write=False):

    signal_returned = 0

    parameters = [0,0,0]

    if which_conductance=="g1":

        parameters = [-1, 9, 2]

    
    if which_conductance=="g2":

        parameters = [-1, 9, 2]

    if write:

        t = np.linspace(initial_time,final_time,time_steps)

        if shape=='triangular':
        
            waveform = triangular_potential(t)

        if shape=='square':

            waveform = square_potential(t, parameters)

        voltage_file = open(f"{DATA_PATH}voltage_file.txt", "w")
        for time_index in range(len(t)):
            voltage_file.write(f'{t[time_index]} \t {waveform[time_index]} \n')

    if shape=='triangular':

        signal_returned = triangular_potential(time)

    if shape=='square':

        signal_returned = square_potential(time, parameters)

    return signal_returned

