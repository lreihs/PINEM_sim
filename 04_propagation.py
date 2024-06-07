import numpy as np
import matplotlib.pyplot as plt
import importlib
import scipy.constants as const
import matplotlib as mpl
from matplotlib import cm

func = importlib.__import__('02_functions')

z = func.z
N_steps = int(func.N_steps)
dt = func.dt * 10**-15
psi_p_evol = []

t_L = 5*10**-10

fig, ax = plt.subplots()
ax.plot(func.z_vt(z, t_L), func.electron_density(z, t_L))
plt.show()
