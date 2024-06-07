import numpy as np
import matplotlib.pyplot as plt
from read_input import * 

### 2.3 Functions modelling the system ###

def psi_in(t):
    return np.exp(-t**2/(2*pulse_duration**2)) 

def spec_in_t(t):
    psi_w_squared =np.abs(np.fft.fftshift(np.fft.fft(psi_in(t))))**2
    return psi_w_squared / max(psi_w_squared)

def A(t):
    return np.exp(-2j * np.abs(g) * np.sin(w_opt*t + np.angle(g)))

def psi_t_0(t):
    return A(t) * psi_in(t)

def broadening_fct():
    return np.exp(-w**2/((np.amax(w)-np.amin(w))/1000)**2)

def psi_w_0(t):
    return np.fft.fftshift(np.fft.fft(psi_t_0(t))) 

def spec(t):
    psi_squared = np.abs(psi_w_0(t))**2
    return psi_squared / np.amax(psi_squared)

def propagator(delta_t):
    return np.exp(1j*w*delta_t)

def psi_w_t(t, delta_t):
    # psi_w_0(t) * propagator(delta_t)
    return np.sqrt(spec(t)) * propagator(delta_t)

def spec_t(t, delta_t):
    psi_squared = np.abs(psi_w_t(t, delta_t)**2)
    return psi_squared / np.amax(psi_squared)

def psi_z_t(t, delta_t):
    return np.fft.ifft(psi_w_t(t, delta_t))

def z_vt(delta_t):
    return np.fft.fftshift(np.fft.fftfreq(N_samp, d = np.abs(w[1]-w[0]))) + v_el*delta_t

def electron_density(t, delta_t):
    psi_squared = np.abs(psi_z_t(t, delta_t)**2)
    return psi_squared / np.amax(psi_squared)

### 2.5 Troubleshooting ###

# fig, ax = plt.subplots()
# ax.plot(t, psi_in(t))
# ax.plot(t, psi_in_corr(t))

fig, ax = plt.subplots()
ax.plt(photon_order, spec(t))

delta_t = 1*10**-10

# fig, ax = plt.subplots()
# ax.plot(z_vt(delta_t), electron_density(t, delta_t))
# # ax.set(ylim=(-1e-20,1e-19))
# plt.show()

# fig, ax = plt.subplots()
# ax.plot(photon_order, spec_t(t, delta_t))
# plt.show()

# fig, ax = plt.subplots()
# ax.plot(w, np.sin(w*10**-16))

plt.show()



