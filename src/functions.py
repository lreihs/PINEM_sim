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
    return psi_squared / np.amax(psi_squared)

def propagator(delta_t):
    # p_sqared = hbar**2*w**2/(c**2)-m_e**2*c**2
    # np.exp(1j*p_sqared*delta_t/(2*m_e*hbar))

    # p = np.sign(w)*np.sqrt(2*m_e*hbar*np.abs(w))
    # np.exp(-1j*v_el*p*delta_t/hbar)
    return np.exp(1j*w*delta_t)

def psi_w_t(t, delta_t):
    return psi_w_0(t) * propagator(delta_t)

def psi_z_t(t, delta_t):
    return np.fft.ifft(psi_w_t(t, delta_t))

# def z_vt(delta_t):
#     return np.fft.fftshift(np.fft.fftfreq(N_samp, d = np.abs(w[1]-w[0])))+v_el*delta_t

def electron_density(t, delta_t):
    psi_squared = np.abs(psi_z_t(t, delta_t)**2)
    return psi_squared / np.amax(psi_squared)

