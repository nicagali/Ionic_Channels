from scipy import signal
import numpy as np
from parameters import *

def potential(time, shape, write=False):

    signal_returned = 0

    if write:

        t = np.linspace(initial_time,final_time,time_steps)

        if shape=='triangular':
        
            waveform = signal.sawtooth(2*np.pi*frequency*t-np.pi/2, 0.5)

        # if shape=='step_free_desired':



        voltage_file = open(f"{DATA_PATH}voltage_file.txt", "w")
        for time_index in range(len(t)):
            voltage_file.write(f'{t[time_index]} \t {waveform[time_index]} \n')

        if shape=='triangular':

            signal_returned = signal.sawtooth(2*np.pi*frequency*time-np.pi/2, 0.5)

    return signal_returned

