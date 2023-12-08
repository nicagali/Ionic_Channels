from scipy import signal
import numpy as np
from parameters import *

def triangular_potential(time):
    return signal.sawtooth(2*np.pi*frequency*time-np.pi/2, 0.5)

def square_potential(time):
    return (-1)*(signal.square(2*np.pi*frequency*time, 0.2) + 1)*5/2

def potential(time, shape, write=False):

    signal_returned = 0

    if write:

        t = np.linspace(initial_time,final_time,time_steps)

        if shape=='triangular':
        
            waveform = triangular_potential(t)

        if shape=='square':

            waveform = square_potential(t)

        voltage_file = open(f"{DATA_PATH}voltage_file.txt", "w")
        for time_index in range(len(t)):
            voltage_file.write(f'{t[time_index]} \t {waveform[time_index]} \n')

    if shape=='triangular':

        signal_returned = triangular_potential(time)

    if shape=='square':

        signal_returned = square_potential(time)

    return signal_returned

