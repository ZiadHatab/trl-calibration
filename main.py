"""
@author: Ziad Hatab (zi.hatab@gmail.com)

Example of using TRL calibration. 
"""

import os

# pip install numpy scikit-rf matplotlib -U
import numpy as np
import skrf as rf   # to load S-parameters from .s2p files
import matplotlib.pyplot as plt # for plotting

import TRL # this should be located in same folder as this script.

# useful functions 
c0 = 299792458   # speed of light in vacuum (m/s)
mag2db = lambda x: 20*np.log10(abs(x))  
db2mag = lambda x: 10**(x/20)           
gamma2ereff = lambda x,f: -(c0/2/np.pi/f*x)**2  # convert gamma to effective relative permittivity
ereff2gamma = lambda x,f: 2*np.pi*f/c0*np.sqrt(-(x-1j*np.finfo(float).eps)) # the eps is to ensure positive square-root
gamma2dbmm  = lambda x: mag2db(np.exp(x.real*1e-3))  # convert losses (real part of gamma) in db/mm

# retrieve the path of this file and modify it to point to the s2p folder 
s2p_path = os.path.dirname(os.path.realpath(__file__)) + '\\s2p\\'

# load uncalibrated S-parameters
thru = rf.Network(s2p_path + 'thru.s2p')       # microstrip with impedance approx 50 ohm
line = rf.Network(s2p_path + 'line_15mm.s2p')  # microstrip with impedance approx 50 ohm
open_A = rf.Network(s2p_path + 'open_A.s1p')
open_B = rf.Network(s2p_path + 'open_B.s1p')

dut = line  # I will update the measurements later with a proper dut.

# switch terms
switch_f = rf.Network(s2p_path + 'sw_forward.s1p')
switch_r = rf.Network(s2p_path + 'sw_reverse.s1p')

f = thru.frequency.f
line_length = 15e-3  # in meters
gamma_est   = ereff2gamma(2.6, f[0])  # use effective relative permittivity as initial guess
reflect_est = 1 # open  standard. 
# Hint: if the reflect is offseted, then you should include an offset by multiplying the estimated reflect by np.exp(-2*gamma_est*offset)

# iterate through frequencies
As = []
Bs = []
ks = []
gammas = []
reflects = []
for thru_S,line_S,reflect_A,reflect_B,GF,GR in zip(thru.s.squeeze(), line.s.squeeze(), open_A.s.squeeze(), open_B.s.squeeze(), switch_f.s.squeeze(), switch_r.s.squeeze()):
    thru_S = TRL.correct_switch(thru_S, GF,GR)
    line_S = TRL.correct_switch(line_S, GF,GR)
    A,B,k,gamma_est,reflect_est = TRL.trl(thru_S, line_S, line_length, gamma_est, reflect_A, reflect_B, reflect_est)
    As.append(A)
    Bs.append(B)
    ks.append(k)
    gammas.append(gamma_est)
    reflects.append(reflect_est)
As = np.array(As)
Bs = np.array(Bs)
ks = np.array(ks)
gammas = np.array(gammas)
reflects = np.array(reflects)

# apply cal (also shift plane and change impedance)
dut_cal_S = []
for dut_S,A,B,k,GF,GR,gamma in zip(dut.s.squeeze(),As,Bs,ks,switch_f.s.squeeze(),switch_r.s.squeeze(),gammas):
    dut_S = TRL.correct_switch(dut_S, GF, GR)
    # A,B,k = TRL.shift_plane(A,B,k,line_length/2,gamma) # shift the reference plane to the middle of the Line standard
    # A,B,k = TRL.change_impedance(A,B,k,50,75)  # change reference impedance from 50 to 75 ohm
    dut_cal_S.append(TRL.apply_cal(dut_S,A,B,k))
dut_cal_S = np.array(dut_cal_S)

# plots
plt.figure()
plt.plot(f*1e-9, mag2db(dut_cal_S[:,0,0]))
plt.xlabel('Frequency (GHz)')
plt.ylabel('S11 (dB)')
plt.title('Calibrated DUT')

plt.figure()
plt.plot(f*1e-9, mag2db(dut_cal_S[:,1,0]))
plt.xlabel('Frequency (GHz)')
plt.ylabel('S21 (dB)')
plt.title('Calibrated DUT')

plt.figure()
plt.plot(f*1e-9, gamma2ereff(gammas, f).real)
plt.xlabel('Frequency (GHz)')
plt.ylabel('Effect relative permittivity')

plt.figure()
plt.plot(f*1e-9, gamma2dbmm(gammas))
plt.xlabel('Frequency (GHz)')
plt.ylabel('Losses (dB/mm)')

plt.figure()
plt.plot(f*1e-9, mag2db(reflects))
plt.xlabel('Frequency (GHz)')
plt.ylabel('S11 (dB)')
plt.title('Calibrated reflect standard')

plt.show()

# EOF