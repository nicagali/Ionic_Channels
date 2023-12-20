import sys
import matplotlib.pyplot as plt
sys.path.append('../codes/')
import odeint_solver
import plotting_functions
from parameters_channels import *
import potential_shapes

print(g_0)

print(initial_time, inital_condition)

odeint_solver.euler_forward_solver('square', conductance_numb='g2')
potential_shapes.potential(0, 'square', which_conductance='g2', write=True)
fig, ax = plt.subplots(1,2, figsize=(10,5))
plotting_functions.plot_voltage(ax[0])
plotting_functions.plot_gsolution(ax[1], solver = 'euler', time_in_steps=True, g0_ratio=False)
fig.tight_layout()
plt.savefig(f'{PLOT_PATH}voltage_gt.pdf')