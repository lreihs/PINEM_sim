import numpy as np
import matplotlib.pyplot as plt
import scipy.constants as const
import matplotlib as mpl
from matplotlib import cm
import warnings
#warnings.filterwarnings('ignore')
import importlib
func = importlib.__import__('02_functions')

z = func.z
k = func.k(z)

### 3.1 Single spectrum ###

spec = func.spec(z)
fig, ax = plt.subplots()
ax.plot(func.photon_order(z), spec)
ax.set(xlabel = 'photon order', ylabel = 'normlaized spectral density')
plt.show()


### 3.2 Field dependent spectra ###

fields = np.linspace(func.min_field, func.max_field, func.N_spec)
specs = []
X, Y = np.meshgrid(fields, k)

for field_strength in fields:
    func.drawProgressBar(list(fields).index(field_strength)/len(fields))
    func.field = field_strength
    specs.append(np.abs(func.psi_p(z))**2 / max(np.abs(func.psi_p(z))**2))

specs = np.array(specs)
specs = specs.T

fig, ax = plt.subplots()
cs = ax.contourf(X, Y, specs, cmap = cm.Reds)
cbar = fig.colorbar(cs)
plt.show()


