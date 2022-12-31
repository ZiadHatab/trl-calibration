"""
@author: Ziad Hatab (zi.hatab@gmail.com)

A script to compute the length and frequency limits of TRL kit.
The length computed here is of the Line standard is in reference to the Thru standard.
That is, the Thru standard by definition has zero-length.
"""

import numpy as np  # pip install numpy -U


def freqs2length(fmin, fmax, ereff, phi=20):
    """
    Compute the required line length given the frequency limits.
    The phase margin (phi) is not enforced. The computed length comes with its exact phase margin.

    Parameters
    -------
    fmin  : float scalar; lower frequency limit in Hz.
    fmax  : float scalar; upper frequency limit in Hz.
    ereff : float scalar; effective relative permittivity (real part only).
    phi   : float scalar; desired phase margin in degrees (between 0 and 90 degrees).

    Returns
    -------
    list_of_lengths : list of of tuples of length and corresponding phase margin (l, phi).

    """
    c0 = 299792458   # speed of light in vacuum (m/s)
    q = fmin/fmax
    list_of_lengths = []
    n = int( (q-(q+1)*phi/180)/(1-q) )
    for m in range(n+1):
        phi = 180*(m*q-m+q)/(q+1)
        l = c0/2/np.sqrt(ereff.real)/fmin*(m+phi/180)
        list_of_lengths.append((l, phi))

    return list_of_lengths

def length2freqs(l, ereff, num_bands=2, phi=20):
    """
    Compute the frequency limits of a line standard for a given length.
    The output is given for multiple frequency bands (phase wrapping). 
    The phase margin (phi) is enforced and all frequency bands exactly satisfy the given phase margin.

    Parameters
    -------
    l         : float scalar; length in meters.
    ereff     : float scalar; effective relative permittivity (real part only).
    phi       : float scalar; specified phase margin in degrees (between 0 and 90 degrees).
    num_bands : int scalar; number of frequency bands to return.

    Returns
    -------
    list_of_freqs : list of of tuples of frequencies and the band index (fmin,fmax,n).

    """
    c0 = 299792458   # speed of light in vacuum (m/s)
    list_of_freqs = []
    for n in range(num_bands):
        fmin = c0/2/np.sqrt(ereff.real)/l*(n+phi/180)
        fmax = c0/2/np.sqrt(ereff.real)/l*(n+1-phi/180)
        list_of_freqs.append((fmin,fmax,n))

    return list_of_freqs

if __name__=="__main__":

    # example of WR-42 waveguide
    ereff = 0.5
    fmin  = 18e9    # in Hz
    fmax  = 26.5e9  # in Hz
    phi   = 20      # in deg (phase margin)
    print(f'Line lengths and phase margins for fmin = {fmin*1e-9:.2f} GHz, fmax = {fmax*1e-9:.2f} GHz, ereff = {ereff:.2f}:')
    for x in freqs2length(fmin, fmax, ereff, phi=phi):
        print(f'l = {x[0]*1e3:.2f} mm, phi = {x[1]:.2f} deg')

    print('-----------------')
    # example of microstrip line
    ereff = 2.6
    phi   = 20      # in deg (phase margin)
    l     = 15e-3   # in meters
    num_bands = 4
    print(f'Frequency bands for l = {l*1e3:.2f} mm, ereff = {ereff:.2f}, phi = {phi:.2f} deg:')
    for x in length2freqs(l, ereff, num_bands=num_bands, phi=phi):
        print(f'fmin = {x[0]*1e-9:.2f} GHz, fmax = {x[1]*1e-9:.2f} GHz')
    
    # EOF
