import numpy as np
from read_input import *

def psi_in(t):
    return np.exp(-t**2/(2*pulse_duration**2)) 

def A(t):
    return np.exp(-2j * np.abs(g) * np.sin(w_opt*t + np.angle(g)))

def psi_t_0(t):
    return A(t) * psi_in(t)

def psi_w_0(t):
    return np.fft.fftshift(np.fft.fft(psi_t_0(t))) 

def spec_0(t):
    psi_squared = np.abs(psi_w_0(t))**2
    return psi_squared / np.trapz(psi_squared, w)

def propagator(delta_t):
    return np.exp(1j*np.pi*(hbar*w)**2/((hbar*w_opt)**2)*delta_t)

def psi_w_t(t, delta_t):
    return psi_w_0(t) * propagator(delta_t)

def psi_z_t(t, delta_t):
    return np.fft.ifft(psi_w_t(t, delta_t))

def electron_density(t, delta_t):
    return np.abs(psi_z_t(t, delta_t))**2

def total_prob(t, delta_t):
    return np.trapz(np.abs(psi_z_t(t, delta_t))**2, t)

