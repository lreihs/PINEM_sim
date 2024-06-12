import numpy as np
from read_input import *

current_path = os.getcwd()
end_index = current_path.find('PINEM_sim') + len('PINEM_sim')
anchor = current_path[:end_index]

def psi_in():
    results = np.exp(-t**2/(2*pulse_duration**2))
    with open(f'{anchor}/docs/psi_in.log', 'w') as f:
        for result in results:
            f.write(f'{result}\n')

def A():
    return np.exp(-2j * np.abs(g) * np.sin(w_opt*t + np.angle(g)))

def psi_t_0():
    psi_in()
    psi_in_ = np.loadtxt(f'{anchor}/docs/psi_in.log', dtype=complex)
    results = A() * psi_in_
    with open(f'{anchor}/docs/psi_t_0.log', 'w') as f:
        for result in results:
            f.write(f'{result}\n')

def psi_w_0():
    psi_t_0()
    psi_t_0_ = np.loadtxt(f'{anchor}/docs/psi_t_0.log', dtype=complex)
    results = np.fft.fftshift(np.fft.fft(psi_t_0_))
    with open(f'{anchor}/docs/psi_w_0.log', 'w') as f:
        for result in results:
            f.write(f'{result}\n')

def spec_0():
    psi_w_0()
    psi_w_0_ = np.loadtxt(f'{anchor}/docs/psi_w_0.log', dtype=complex)
    psi_squared = np.abs(psi_w_0_)**2
    return psi_squared / np.amax(psi_squared)

def propagator(delta_t):
    return np.exp(1j*w*delta_t)

def psi_w_t(delta_t):
    psi_w_0()
    psi_w_0_ = np.loadtxt(f'{anchor}/docs/psi_w_0.log', dtype=complex)
    results = psi_w_0_ * propagator(delta_t)
    with open(f'{anchor}/docs/psi_w_t.log', 'w') as f:
        for result in results:
            f.write(f'{result}\n')

def psi_z_t(delta_t):
    psi_w_t(delta_t)
    psi_w_t_ = np.loadtxt(f'{anchor}/docs/psi_w_t.log', dtype=complex)
    results = np.fft.ifft(psi_w_t_)
    with open(f'{anchor}/docs/psi_z_t.log', 'w') as f:
        for result in results:
            f.write(f'{result}\n')

def electron_density(delta_t):
    psi_z_t(delta_t) 
    psi_z_t_ = np.loadtxt(f'{anchor}/docs/psi_z_t.log', dtype=complex)
    return np.abs(psi_z_t_)**2

def total_prob(t, delta_t):
    return np.trapz(np.abs(psi_z_t(t, delta_t))**2, t)