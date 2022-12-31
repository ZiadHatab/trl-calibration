"""
@author: Ziad Hatab (zi.hatab@gmail.com)

This is an implementation of the TRL (Thru-Reflect-Line) calibration. 
To learn more about mathematical derivation of the calibration, see my post on the topic:
https://ziadhatab.github.io/posts/trl-calibration  

##########-NOTE-##########
This script is written as collection of functions to process only one frequency point. 
Therefore, you need to import this script in your main script and iterate through all frequency points.
##########-END-##########
"""

import numpy as np   # pip install numpy -U

def s2t(S, norm=False):
    '''
    Convert 2x2 S-parameters to T-parameters.

    If norm is set to True, the output is normalized such that T22 = 1. 
    This is not a correct conversion, but it is useful when applying error-box calibration.
    '''
    T = S.copy()
    T[0,0] = -(S[0,0]*S[1,1]-S[0,1]*S[1,0])
    T[0,1] = S[0,0]
    T[1,0] = -S[1,1]
    T[1,1] = 1
    return T if norm else T/S[1,0]

def t2s(T, norm=False):
    '''
    Convert 2x2 T-parameters to S-parameters.
    
    If norm is set to True, the output is normalized such that S21 = 1. 
    This is not a correct conversion, but it is useful when applying error-box calibration.
    '''
    S = T.copy()
    S[0,0] = T[0,1]
    S[0,1] = T[0,0]*T[1,1]-T[0,1]*T[1,0]
    S[1,0] = 1
    S[1,1] = -T[1,0]
    return S if norm else S/T[1,1]

def correct_switch(S, GF, GR):
    '''
    S  : 2x2 complex array; S-parameters
    GF : complex scalar; forward switch-term
    GR : complex scalar; reverse switch-term

    ####-NOTE-####
    Equations are taken from Eqs. (18)-(21) in below reference:
    R. B. Marks, "Formulations of the Basic Vector Network Analyzer Error Model including Switch-Terms,"
    50th ARFTG Conference Digest, 1997, pp. 115-126, doi: 10.1109/ARFTG.1997.327265.
    https://ieeexplore.ieee.org/document/4119948 
    ####-END-####
    '''
    S_new = S.copy()
    S_new[0,0] = (S[0,0]-S[0,1]*S[1,0]*GF)/(1-S[0,1]*S[1,0]*GF*GR)
    S_new[0,1] = (S[0,1]-S[0,0]*S[0,1]*GR)/(1-S[0,1]*S[1,0]*GF*GR)
    S_new[1,0] = (S[1,0]-S[1,1]*S[1,0]*GF)/(1-S[0,1]*S[1,0]*GF*GR)
    S_new[1,1] = (S[1,1]-S[0,1]*S[1,0]*GR)/(1-S[0,1]*S[1,0]*GF*GR)
    return S_new

def trl(thru_S, line_S, line_length, gamma_est, reflect_A, reflect_B, reflect_est=-1):
    '''
    thru_S      : 2x2 complex array; the measured S-parameters of the Thru standard.
    line_S      : 2x2 complex array; the measured S-parameters of the Line standard.
    line_length : float scalar; the length of the line standard in meters (in reference to the Thru).
    gamma_est   : complex scalar; estimated propagation constant
    reflect_A   : complex scalar, measured reflection coefficient at the left port (aka port A)
    reflect_B   : complex scalar, measured reflection coefficient at the right port (aka port B)
    reflect_est : complex scalar, estimated reflection coefficient of the reflect standard (e.g., short: -1, open: 1)

    ####-NOTE-####
    It is assumed that the S-parameters are already corrected for the switch terms.
    ####-EMD-####
    '''
    ## Convert S- to T-parameters
    line_T = s2t(line_S)
    thru_T = s2t(thru_S)
    thru_T_inv = np.linalg.inv(thru_T)
    
    ## Start solving for the normalized coefficients
    # left port, i.e., A
    eval, evec = np.linalg.eig(line_T@thru_T_inv)
    mininx = np.argmin( abs(eval-np.exp(-gamma_est*line_length)) ) # get index of exp(-gamma*l)
    v1 = evec[:,mininx]    # eigenvector corresponds to exp(-gamma*l)
    v2 = evec[:,~mininx]   # eigenvector corresponds to exp(+gamma*l)
    a21_a11 = v1[1]/v1[0]
    a12 = v2[0]/v2[1]
    # right port, i.e., B
    eval, evec = np.linalg.eig((thru_T_inv@line_T).T)
    mininx = np.argmin( abs(eval-np.exp(-gamma_est*line_length)) ) # get index of exp(-gamma*l)
    u1 = evec[:,mininx]    # eigenvector corresponds to exp(-gamma*l)
    u2 = evec[:,~mininx]   # eigenvector corresponds to exp(+gamma*l)
    b12_b11 = u1[1]/u1[0]
    b21 = u2[0]/u2[1]
    # build normalized error-boxes 
    A_ = np.array([ [1, a12], [a21_a11, 1] ])
    B_ = np.array([ [1, b12_b11], [b21, 1] ])

    ## New estimated propagation constant
    E = (eval[mininx]  + 1/eval[~mininx])/2 # average for exp(-gamma*l)
    P = np.round( (gamma_est*line_length + np.log(E)).imag/2/np.pi ) # phase unwrap factor
    gamma = (-np.log(E) + 1j*2*np.pi*P)/line_length
    
    ## Compute k and a11b11 from the thru measurements
    ka11b11,_,_,k = ( np.linalg.pinv(A_)@thru_T@np.linalg.pinv(B_) ).flatten()
    a11b11 = ka11b11/k   # a11*b11
    
    ## solve for a11/b11, and hence a11 and b11
    a11_b11 = (reflect_A - a12)/(1 - reflect_A*a21_a11)*(1 + reflect_B*b12_b11)/(reflect_B + b21)  # a11/b11
    a11 = np.sqrt(a11_b11*a11b11)
    b11 = a11b11/a11
    reflect_cal_A = (reflect_A - a12)/(1 - reflect_A*a21_a11)/a11  # calibrated reflect measurement from port A
    reflect_cal_B = (reflect_B + b21)/(1 + reflect_B*b12_b11)/b11  # calibrated reflect measurement from port B
    reflect_cal = (reflect_cal_A + reflect_cal_B)/2  # average from both port measurements
    # choose the answer closest to the the estimate
    if abs(reflect_cal - reflect_est) > abs(-reflect_cal - reflect_est):
        a11 = -a11
        b11 = -b11
        reflect_cal = -reflect_cal
    
    ## De-normalize the error-boxes
    A = A_@np.array([[a11, 0], [0, 1]])
    B = np.array([[b11, 0], [0, 1]])@B_
    
    return A, B, k, gamma, reflect_cal  
    # Hint: gamma and reflect_cal can be passed back to the function for the next frequency point.

def apply_cal(S_dut, A, B, k, left=True):
    '''
    apply calibration to a DUT.
    S_dut : 2x2 complex array; DUT S-parameters
    A     : 2x2 complex array; left error-box (T-parameters)
    B     : 2x2 complex array; right error-box (T-parameters)
    k     : complex scalar;   7-term of the calibration
    left  : boolean; True means using left port (A) if single port DUT was provided, else right port (B).

    ####-NOTE-####
    It is assumed that the S-parameters are already corrected for the switch terms.
    ####-EMD-####
    '''
    nports = np.sqrt(S_dut.size).astype('int') # number of ports

    # if 1-port, convert to 2-port (later convert back to 1-port)
    if nports < 2:
        S_dut = np.eye(2)*S_dut  # i.e., S11 = S22
    
    S_new = S_dut.copy()
    
    # apply cal
    T_dut_norm = s2t(S_dut, norm=True)
    T_dut_cal_norm = np.linalg.inv(A)@T_dut_norm@np.linalg.inv(B)
    s21_cal = k*S_dut[1,0]/T_dut_cal_norm[1,1]
    T_dut_cal_norm = T_dut_cal_norm/T_dut_cal_norm[1,1]  # set last element to 1
    
    S_new[0,0] =  T_dut_cal_norm[0,1]
    S_new[1,0] =  s21_cal
    S_new[0,1] =  (T_dut_cal_norm[0,0] - T_dut_cal_norm[1,0]*T_dut_cal_norm[0,1])/s21_cal
    S_new[1,1] = -T_dut_cal_norm[1,0]

    # revert to 1-port device if the input was a 1-port device
    if nports < 2:
        S_new = S_new[0,0] if left else S_new[1,1]

    return S_new


## Below are auxiliary functions that could be used to update the error-boxes if desired.

def shift_plane(A,B,k,offset,gamma):
    '''
    Shift the reference plane of the calibration.
    A      : 2x2 complex array; left error-box (T-parameters)
    B      : 2x2 complex array; right error-box (T-parameters)
    k      : complex scalar; 7-term of the calibration
    offset : float scalar; positive offset away from ports; negative offset towards the ports
    gamma  : complex scalar; propagation constant of the line standards.
    '''
    exp  = np.exp(-gamma*offset)
    L = np.array([[exp, 0], [0, 1/exp]])

    A_new = A@L # temp
    B_new = L@B # temp

    k_new = k*A_new[1,1]*B_new[1,1]
    A_new = A_new/A_new[1,1]  # set last element to 1
    B_new = B_new/B_new[1,1]  # set last element to 1

    return A_new, B_new, k_new

def Qnm(Zn, Zm):
    '''
    Impedance transformer in T-parameters from on Eqs. (86) and (87) in
    R. Marks and D. Williams, "A general waveguide circuit theory," 
    Journal of Research (NIST JRES), National Institute of Standards and Technology,
    Gaithersburg, MD, no. 97, 1992.
    https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4914227/

    Zn : complex scalar; impedance from which we are transforming
    Zm : complex scalar; impedance to which we are transforming to.
    '''
    Gnm = (Zm-Zn)/(Zm+Zn)
    return np.sqrt(Zn.real/Zm.real*(Zm/Zn).conjugate())/np.sqrt(1-Gnm**2)*np.array([[1, Gnm],[Gnm, 1]])

def change_impedance(A,B,k,Z_line,Z_new):
    '''
    Change the reference impedance of the calibration from Z_line to Z_new.
    A      : 2x2 complex array; left error-box (T-parameters).
    B      : 2x2 complex array; right error-box (T-parameters).
    k      : complex scalar; 7-term of the calibration.
    Z_line : complex scalar; the characteristic impedance of the Line standard.
    Z_new  : complex scalar; the new reference impedance.
    '''
    A_new = A@Qnm(Z_line,Z_new) # temp
    B_new = Qnm(Z_new,Z_line)@B # temp

    k_new = k*A_new[1,1]*B_new[1,1]
    A_new = A_new/A_new[1,1]  # set last element to 1
    B_new = B_new/B_new[1,1]  # set last element to 1

    return A_new, B_new, k_new
