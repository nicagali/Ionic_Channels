from scipy import signal
import numpy as np
from parameters_channels import *

def triangular_potential(time):
    return signal.sawtooth(2*np.pi*frequency_tr*time-np.pi/2, 0.5)

def square_potential(time, duty_sq, amplitude_sq, phase_sq, shift_sq):
    return amplitude_sq*signal.square(2*np.pi*frequency_sq*time + phase_sq, duty_sq) + shift_sq

def sine_potential(time, amplitude_sine, freq_sine):
    return amplitude_sine*np.sin(2 * np.pi * freq_sine * time)

def cosine_potential(time, amplitude_sine, freq_sine):
    return amplitude_sine*np.cos(2 * np.pi * freq_sine * time)

def potential(time, shape, which_conductance=0, write=False):

    signal_returned = 0

    parameters = [0,0,0,0]

    if which_conductance=="g1":

        duty_sq = 0.2
        frequency_sq = 0.1
        period_sq=1/frequency_sq
        amplitude_sq = -1/2
        phase_sq = 0
        shift_sq = -4.5
    
    if which_conductance=="g2":

        duty_sq = 0.1
        frequency_sq = 0.1
        period_sq=1/frequency_sq
        amplitude_sq = -2
        # phase_sq = np.pi*frequency_sq
        phase_sq = 0
        shift_sq = 2

    if write:

        t = np.linspace(initial_time,final_time,time_steps)

        if shape=='triangular':
        
            waveform = triangular_potential(t)

        if shape=='square':

            waveform = square_potential(t, parameters)

        if shape=='sine':

            waveform = sine_potential(t, amplitude_sine, freq_sine)

        voltage_file = open(f"{DATA_PATH}voltage_file.txt", "w")
        for time_index in range(len(t)):
            voltage_file.write(f'{t[time_index]} \t {waveform[time_index]} \n')

    if shape=='triangular':

        signal_returned = triangular_potential(time)

    if shape=='square':

        signal_returned = square_potential(time, duty_sq, amplitude_sq, phase_sq, shift_sq)

    return signal_returned

