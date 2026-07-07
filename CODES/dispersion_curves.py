# Functions 

## Dispersion curve estimative 

import numpy as np
from disba import PhaseDispersion,EigenFunction,PhaseSensitivity
import pandas as pd
from scipy.ndimage import gaussian_filter
import matplotlib.pyplot as plt

from CODES.modeling import create_velocity_model_from_profile_vs

# -----------------------------------------------------------

def estimate_disp_from_velocity_model(vel_mol, min_freq, max_freq, number_samples):
    """
    Calculates the fundamental mode Rayleigh wave phase velocity dispersion curve 
    from a 1-D velocity model using the disba library.

    Parameters
    ----------
    vel_mol : numpy.ndarray
        A 2-D array of shape (N, 4) representing the velocity model.
        Each row corresponds to a layer with the following columns:
        [thickness, velocity_p, velocity_s, density]
            - Layer thickness (in km).
            - Layer P-wave velocity (in km/s).
            - Layer S-wave velocity (in km/s).
            - Layer density (in g/cm³).
    min_freq : float
        The minimum frequency for the dispersion curve calculation, in Hertz.
    max_freq : float
        The maximum frequency for the dispersion curve calculation, in Hertz.
    number_samples : int
        The number of frequency samples to generate evenly between min_freq 
        and max_freq.
    algorithm_str : str, optional
        The root-finding algorithm to be used by disba (e.g., 'dunkin' or 'thomson'). 
        Default is 'dunkin'.

    Returns
    -------
    cpr : disba.DispersionCurve
        An object containing the computed dispersion curve data for the fundamental 
        Rayleigh wave mode. It typically contains `.period` and `.velocity` attributes.

    Examples
    --------
    >>> import numpy as np
    >>> from disba import PhaseDispersion
    >>> velocity_model = np.array([
    ...    [10.0, 7.00, 3.50, 2.00],
    ...    [10.0, 6.80, 3.40, 2.00],
    ...    [10.0, 7.00, 3.50, 2.00],
    ...    [10.0, 7.60, 3.80, 2.00]
    ... ])
    >>> # Estimate dispersion from 10 Hz to 100 Hz
    >>> cpr = estimate_disp_from_velocity_model(velocity_model, 10.0, 100.0)
    """

    # Generate a linear array of frequencies from min_freq to max_freq
    hz = np.linspace(min_freq, max_freq, number_samples) 
    
    # Convert frequency (Hertz) to period (seconds).
    # Disba requires periods to be sorted in ascending order (low to high).
    # Reversing the frequency array [::-1]
    t = 1 / hz[::-1] 
    
    # Initialize the PhaseDispersion object from disba.
    # *vel_mol.T transposes the 2D array and unpacks it into four separate 1D arrays 
    # (thickness, Vp, Vs, density) as required by disba.
    # dc=0.001 sets the phase velocity phase search step size.
    pdisp = PhaseDispersion(*vel_mol.T, dc=0.001)
    
    # Compute the phase dispersion curve for the given periods.
    # mode=0 corresponds to the fundamental mode.
    # wave="rayleigh" specifies that we are modeling Rayleigh waves (instead of Love waves).
    cpr = pdisp(t, mode=0, wave="rayleigh")

    return cpr

# -----------------------------------------------------------

def compute_dispersion(row,vs_col,depth_col,min_freq,max_freq,number_samples):
    """
    Computes a simulated Rayleigh wave phase dispersion curve for a single 
    geological profile (row of a DataFrame).

    This function extracts the necessary physical profile parameters, builds a 
    1-D velocity model, estimates the dispersion curve, and formats the output 
    with standard units (m/s and Hz).

    Parameters
    ----------
    row : pandas.Series or dict-like
        A single row from a pandas DataFrame containing profile data.
    vs_col : str
        The column or key name containing the shear wave velocity (Vs) array. 
    depth_col : str
        The column or key name containing the layer thickness array. 

    Returns
    -------
    pandas.Series
        A Series containing two arrays:
        - 'simulated_dispersion': The calculated phase velocities in m/s.
        - 'simulated_frequency': The corresponding frequencies in Hz.
    """
    
    # Extract the physical property arrays from the provided row/dictionary
    vs = row[vs_col]
    thick = row[depth_col]
    
    layer_model = [[t, v] for t, v in zip(thick, vs)]
    
    # Construct the 1-D velocity model. 
    # Note: Assumes create_velocity_model_from_profile_vs can unpack a list 
    # containing [thickness, vs].
    simulated_velocity_model = create_velocity_model_from_profile_vs(layer_model)
    
    # Estimate the fundamental mode dispersion curve using the generated model
    simulated_cpr = estimate_disp_from_velocity_model(simulated_velocity_model,min_freq,max_freq,number_samples)
    
    # Extract the simulated velocities and convert from standard modeling units 
    # (km/s) back to meters per second (m/s)
    simulated_dispersion = simulated_cpr.velocity * 1000  
    
    # Extract the periods (in seconds) and convert to frequency (in Hertz)
    simulated_frequency = 1 / simulated_cpr.period         
    
    # Return the results neatly packed into a pandas Series for easy DataFrame integration
    return pd.Series({
        'simulated_dispersion': simulated_dispersion[::-1],
        'simulated_frequency': simulated_frequency[::-1]
    })

# -----------------------------------------------------------

def calculate_dispersion_image(data, nt, dt, offsets, c_min, c_max, dc, f_min, f_max):
    
    """
    Calculates the dispersion image using the Phase-Shift Method.

    Based on MASWavesPy package (https://github.com/Mazvel/maswavespy/tree/main):

    - Ólafsdóttir, Elín Ásta & Bessason, Bjarni & Erlingsson, Sigurdur & Kaynia, Amir. (2024). 
      A Tool for Processing and Inversion of MASW Data and a Study of Inter-session Variability of MASW. 
      Geotechnical Testing Journal. 47. 1006-1025. 10.1520/GTJ20230380. 
    
    Parameters
    ----------
    data : numpy.ndarray
        2D numpy array containing the seismic data. Expected shape: (Time Samples, Number of Receivers).
    nt : float
        Total number of discrete time steps.
    dt : float
        Time sampling interval (seconds).
    offsets : numpy.ndarray
        1D array with the distance of each receiver from the source (meters).
    c_min, c_max, dc : float
        Range and step size of the testing phase velocities (m/s).
    f_min, f_max : float
        Frequency range to analyze (Hz).
        
    Returns
    -------
    freqs_interest : numpy.ndarray
        Frequency array analyzed [Hz].
    c_test : numpy.ndarray
        Phase velocity testing array [m/s].
    dispersion_image : numpy.ndarray
        Normalized dispersion image matrix (Energy).
    """
    
    # FFT along the time axis (axis 0)
    data_f = np.fft.rfft(data, axis=0)
    freqs = np.fft.rfftfreq(n=nt, d=dt)
    
    # Filter only the frequencies of interest to optimize computation
    f_idx = np.where((freqs >= f_min) & (freqs <= f_max))[0]
    freqs_interest = freqs[f_idx]
    data_f_interest = data_f[f_idx, :]
    
    # Amplitude normalization (eliminates geometric attenuation effects, focusing solely on the phase)
    # A small value (epsilon) is added to avoid division by zero
    data_f_norm = data_f_interest / (np.abs(data_f_interest) + 1e-12)
    
    # Set up the phase velocity array for testing
    c_test = np.arange(c_min, c_max + dc, dc)
    
    # Matrix to store the dispersion panel (Energy)
    dispersion_image = np.zeros((len(freqs_interest), len(c_test)))
    
    # Phase-Shift and Stacking
    for i, f in enumerate(freqs_interest):
        omega = 2.0 * np.pi * f
        for j, c in enumerate(c_test):
            # Wavenumber for the current testing velocity
            k_test = omega / c
            
            # Apply the phase shift: e^(i * k * x)
            phase_shift = np.exp(1j * k_test * offsets)
            
            # Multiply frequency domain data by the shift and sum along the offsets
            # Equivalent to integrating over dx
            stack = np.sum(data_f_norm[i, :] * phase_shift)
            
            # Energy is the magnitude of the stacked signal
            dispersion_image[i, j] = np.abs(stack)
            
    # Normalize the final image for easier visualization (scale from 0 to 1 per frequency row)
    dispersion_image = dispersion_image / np.max(dispersion_image, axis=1, keepdims=True)
    
    return freqs_interest, c_test, dispersion_image

def pick_dispersion_curve_hessian(dispersion_image, freqs, v_phase, sigma, min_freq,plot_curve=True):

    """
    Extracts the fundamental surface wave dispersion curve from a dispersion image 
    using Hessian matrix eigenvalue analysis. 
    
    This method identifies local energy ridges by evaluating the second-order 
    derivatives (local curvature) of the dispersion power spectrum based on:
        - Xiaoping, H., Jiashun, Y., Jianlong, Y. et al. 
          An automatic algorithm for surface wave dispersion curve picking based on
          Hessian matrix attributes. Sci Rep 15, 21595 (2025).
          https://doi.org/10.1038/s41598-025-04954-w

    Parameters
    ----------
    dispersion_image : numpy.ndarray
        2D array representing the dispersion power spectrum (energy).
    freqs : numpy.ndarray
        1D array of frequency values [Hz].
    v_phase : numpy.ndarray
        1D array of phase velocity values [m/s].
    sigma : float
        Standard deviation for the Gaussian filter used to smooth the image.
    min_freq : float
        Minimum frequency threshold [Hz]. Ridge points below this are discarded to 
        prevent picking low-frequency noise artifacts.
    plot_curve : Boolean    
        To plot the results of the extracted dispersion curve.

    Returns
    -------
    picked_freqs : numpy.ndarray
        1D array of the frequencies for the extracted dispersion curve.
    picked_vels : numpy.ndarray
        1D array of the corresponding phase velocities.
    """
    
    # -------------------
    # Get data parameters
    # -------------------
    shape = dispersion_image.shape

    # Define the coordinates and axes indices based on the presumed input orientation
    grad_coords = [freqs, v_phase]
    f_axis, v_axis  = 0, 1

    # Calculate the average grid spacing for frequency (df) and velocity (dv).
    # This is later used to check if the exact ridge peak falls within the current pixel.
    df = np.mean(np.diff(freqs))
    dv = np.mean(np.diff(v_phase))
    
    # -------------------------------------
    # Data preparation to extract the curve
    # -------------------------------------

    # Apply a 2D Gaussian filter to smooth the energy spectrum.
    # This suppresses small local noise spikes and creates continuous ridge structures.
    z = gaussian_filter(dispersion_image, sigma=sigma)

    # Estimate the first-order derivatives (gradient) using central differences.
    grads = np.gradient(z, *grad_coords)
    zf = grads[f_axis]  # First derivative with respect to frequency
    zv = grads[v_axis]  # First derivative with respect to velocity

    # Estimate the second-order partial derivatives starting from zv
    grads_v = np.gradient(zv, *grad_coords)
    zvf = grads_v[f_axis]
    zvv = grads_v[v_axis]
    
    # Estimate the second-order partial derivatives starting from zf
    grads_f = np.gradient(zf, *grad_coords)
    zff = grads_f[f_axis]
    zfv = grads_f[v_axis]
    
    # ------------------
    # Ridge point search 
    # ------------------
    ridge_points = []
    
    # Iterate over every pixel (grid point) in the dispersion image
    for i in range(shape[0]):
        for j in range(shape[1]):

            idx_f, idx_v = i, j

            # Construct the 2x2 Hessian Matrix (H) for the current grid point.
            # The Hessian describes the local curvature of the energy surface.
            H = np.array([[zff[i, j], zfv[i, j]],
                          [zvf[i, j], zvv[i, j]]])
            
            # Compute eigenvalues and eigenvectors of the Hessian matrix.
            eigvals, eigvecs = np.linalg.eig(H)
            
            # Sort the eigenvalues by absolute magnitude in descending order.
            # lambda_1 (largest absolute value) points in the direction of steepest curvature.
            idx = np.argsort(np.abs(eigvals))[::-1]
            eigvals = eigvals[idx]
            eigvecs = eigvecs[:, idx]
            
            lambda_1, lambda_2 = eigvals[0], eigvals[1]
            n1 = eigvecs[:, 0]     # Principal eigenvector corresponding to lambda_1
            n1f, n1v = n1[0], n1[1]
            
            # ------------
            # Condition 1: 
            # Check if it's a local maximum (ridge) rather than a valley.
            # The primary curvature (lambda_1) must be negative and stronger than the secondary.
            if lambda_1 <= -np.abs(lambda_2):
                
                # Calculate sub-pixel distance 't' to the exact extreme point 
                # along the direction of the principal eigenvector (n1).
                numerator = -(zf[i, j] * n1f + zv[i, j] * n1v)
                denominator = (zff[i, j] * n1f**2 + 2 * zfv[i, j] * n1f * n1v + zvv[i, j] * n1v**2)
                
                if denominator != 0:
                    t = numerator / denominator
                    
                    # ------------
                    # Condition 2: 
                    # Check if this theoretical exact maximum falls 
                    # within the bounds of the current pixel 
                    # (half a grid spacing).
                    if np.abs(n1f * t) <= 0.5 * df and np.abs(n1v * t) <= 0.5 * dv:
                        
                        # If both conditions are met, record the point as a valid ridge candidate.
                        ridge_points.append({
                            'f_idx': idx_f, 
                            'v_idx': idx_v,
                            'f_val': freqs[idx_f],
                            'v_val': v_phase[idx_v],
                            'energy': z[i, j]
                        })
    # -------------------
    # Return empty arrays 
    # if no ridge points were found in the entire image
    if not ridge_points:
        return np.array([]), np.array([])

    # --------------------------------------
    # Curve extraction with frequency cutoff
    # --------------------------------------
    picked_freqs = []
    picked_vels = []
    
    # Get a unique, sorted list of frequency indices that contain valid ridge points
    f_indices = sorted(list(set(p['f_idx'] for p in ridge_points)))
    
    for f_idx in f_indices:

        current_freq = freqs[f_idx]
        
        # Skip this frequency bin entirely if it falls below the user-defined cutoff
        if current_freq < min_freq:
            continue
            
        # Isolate all valid ridge candidates found at this specific frequency
        points_at_f = [p for p in ridge_points if p['f_idx'] == f_idx]
        
        # Select the single ridge point with the highest energy.
        # (This implicitly assumes we are tracking the strongest/fundamental mode).
        best_point = max(points_at_f, key=lambda p: p['energy'])
        
        picked_freqs.append(best_point['f_val'])
        picked_vels.append(best_point['v_val'])


    if plot_curve:

        plt.figure(figsize=(10, 6))

        # Notice the .T added to dispersion_image to transpose the matrix for matplotlib
        plt.pcolormesh(freqs, v_phase, dispersion_image.T, shading='auto', cmap='Spectral_r')

        # Overlay the extracted ridge
        plt.plot(picked_freqs, picked_vels, color='k', linewidth=1.5, ls='--',label='Picked Curve (Hessian)')

        plt.xlabel('Frequency (Hz)')
        plt.ylabel('Phase Velocity (m/s)')
        plt.title('Dispersion Image')
        plt.colorbar(label='Normalized Energy')
        plt.legend()
        plt.gca()
        plt.show()

    return np.array(picked_freqs), np.array(picked_vels)
        
