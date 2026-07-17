import numpy as np
from scipy.signal import butter, filtfilt
from scipy.signal import windows

# ---------------------------------------------------------------------
# 1. FORWARD MODELING HELPER FUNCTIONS
# ---------------------------------------------------------------------

def add_white_noise(trace, noise_percent):
    """
    Add Gaussian white noise to a seismic trace based on a percentage of the
    signal's RMS amplitude.

    Parameters
    ----------
    trace : np.ndarray
        Input seismic trace (1D or 2D array).
    noise_percent : float
        Percentage of noise relative to the signal's RMS amplitude
        (e.g., 1, 2, 5, 10 for 1%, 2%, 5%, 10% noise).

    Returns
    -------
    np.ndarray
        Noisy trace (same shape as input).
    """
    # 1. Compute the RMS amplitude of the signal
    signal_rms = np.sqrt(np.mean(trace**2))

    # 2. Convert percentage to noise standard deviation
    noise_std = (noise_percent / 100.0) * signal_rms

    # 3. Generate Gaussian white noise with that standard deviation
    noise = np.random.normal(0, noise_std, trace.shape)

    return trace + noise

# ---------------------------------------------------------------------


def geometric_spreading_correction(trace, dt, mode="cylindrical", t_ref=None):
    """
    Compensates for geometric-spreading attenuation in the observed trace,
    bringing its amplitude envelope to the same "no-spreading" regime
    assumed by the 1D convolutional model (which has no geometric
    decay built in).

    mode="cylindrical" (2D, gain ~ sqrt(t)): use if Deepwave's scalar()
        was run with v as a 2D grid (most common case).
    mode="spherical" (3D, gain ~ t): use if v was modeled in 3D.

    Parameters
    ----------
    trace : ndarray
        Observed trace (already muted, with direct arrival removed).
    dt : float
        Time sampling interval (s).
    mode : str
        "cylindrical" (2D) or "spherical" (3D).
    t_ref : float, optional
        Reference time (s) used to avoid an exploding gain near t=0
        (where t->0 would make the gain diverge). Default: 1 sample.

    Returns
    -------
    corrected_trace : ndarray (float32)
    """

    nt = len(trace)
    t = np.arange(nt) * dt
    t_ref = t_ref if t_ref is not None else dt
    t_safe = np.maximum(t, t_ref)
    exponent = 0.5 if mode == "cylindrical" else 1.0
    gain = (t_safe / t_ref) ** exponent

    return (trace * gain).astype(np.float32)

# ---------------------------------------------------------------------

def bandpass_filter(trace, dt, f_peak, low_frac=0.3, high_frac=1.0, order=4):
    """
    Applies a zero-phase Butterworth band-pass filter to remove
    low-frequency drift and high-frequency numerical noise that the
    1D convolutional model (with a single Ricker wavelet) cannot
    physically explain.

    The passband is defined relative to the wavelet's peak frequency,
    since that is the only frequency content the forward model can
    reproduce -- energy well outside that band is necessarily noise
    or trend, not usable reflectivity information.

    Parameters
    ----------
    trace : ndarray
        Input trace (1D).
    dt : float
        Time sampling interval (s).
    f_peak : float
        Wavelet peak frequency (Hz), used as the reference for the
        passband limits.
    low_frac, high_frac : float
        Passband limits as a fraction of f_peak (default: 0.3x-3x f_peak).
    order : int
        Butterworth filter order.

    Returns
    -------
    filtered_trace : ndarray (float32)
    """
    fs = 1.0 / dt
    nyq = fs / 2.0

    low = (f_peak * low_frac) / nyq
    high = min((f_peak * high_frac) / nyq, 0.99)

    b, a = butter(order, [low, high], btype="band")
    filtered = filtfilt(b, a, trace)

    return filtered.astype(np.float32)

# ---------------------------------------------------------------------

def apply_edge_taper(trace, taper_samples=200):
    """
    Applies a smooth decay (half of a Hanning window) to the last samples
    of the trace to avoid numerical edge artifacts.
    """
    # np.squeeze ensures the trace is (N,) instead of (N, 1)
    trace_tapered = np.copy(np.squeeze(trace))
    
    # Create the descending half of a Hanning window
    taper = windows.hann(taper_samples * 2)[taper_samples:]
    
    # Multiply the end of the trace by the taper (forcing it to zero)
    trace_tapered[-taper_samples:] *= taper
    
    return trace_tapered
# ---------------------------------------------------------------------

def estimate_wavelet_duration(f_peak, n_cycles=3.0):
    """
    Estimate a physically motivated total duration for a Ricker wavelet,
    based on its peak frequency, instead of an arbitrary fixed value.

    A Ricker wavelet's energy is concentrated within a few periods of its
    peak frequency: lower frequencies require a proportionally longer
    time window to avoid truncation (which introduces spectral leakage
    when the wavelet is convolved with the reflectivity series), while
    higher frequencies need a shorter one. Tying the duration to 1/f_peak
    keeps the wavelet support consistent regardless of what frequency is
    chosen, without requiring a separate user-defined length parameter.

    Parameters
    ----------
    f_peak : float
        Peak (dominant) frequency of the wavelet (Hz).
    n_cycles : float
        Number of periods (1/f_peak) to include in the total duration.
        3.0 is a common practical choice: it captures the main lobe and
        the first side lobes, where the amplitude has already decayed to
        a small fraction of the peak.

    Returns
    -------
    length : float
        Total wavelet duration (s).
    """
    return n_cycles / f_peak

# ---------------------------------------------------------------------

def ricker_wavelet(f_peak, dt):
    """
    Generate a Ricker wavelet (second derivative of a Gaussian).

    Parametrization follows the convention used by Deepwave
    (deepwave.wavelets.ricker): the wavelet is defined by an explicit
    number of samples and an explicit peak time, rather than by a
    symmetric time window built with a floating-point step. This avoids
    off-by-one sample counts caused by floating-point rounding in
    np.arange, and it gives explicit control over where the peak sits,
    which matters for phase alignment when the wavelet is later
    convolved with the reflectivity series.

    Parameters
    ----------
    f_peak : float
        Peak (dominant) frequency of the wavelet (Hz).
    dt : float
        Sampling interval (s).
    
    Returns
    -------
    wavelet : ndarray
        Ricker wavelet sampled at dt, with `n_samples` samples.
    """

    length = estimate_wavelet_duration(f_peak)
    n_samples = int(round(length / dt)) + 1
    peak_time = (n_samples - 1) * dt / 2.0

    t = np.arange(n_samples) * dt - peak_time
    pi2_f2_t2 = (np.pi * f_peak * t) ** 2
    wavelet = (1.0 - 2.0 * pi2_f2_t2) * np.exp(-pi2_f2_t2)

    return wavelet

# ---------------------------------------------------------------------

def mute_direct_arrival(
    trace,
    offset,
    dt,
    freq,
    peak_time=0.0,
    velocity=1500.0,
    flat_len=300,
    taper_len=100,
):
    """
    Applies a cosine-tapered mute to the direct arrival of a 1D seismic trace.

    The direct-arrival time is estimated with a linear moveout:
        t_arrival = peak_time + |offset| / velocity

    A flat (fully zeroed) window is applied around this arrival time,
    surrounded by a smooth cosine taper that transitions the mask
    from 0 back to 1.

    Parameters
    ----------
    trace : np.ndarray
        1D seismic trace, shape (nt,)
    offset : float
        Source-receiver distance (m)
    dt : float
        Sampling interval (s)
    freq : float
        Peak/dominant frequency of the source wavelet (Hz), used to
        estimate a default mute window width
    peak_time : float, optional
        Source wavelet delay (time-to-peak), in seconds
        (e.g. 1.5/freq for a typical Ricker wavelet). Default is 0.0.
    velocity : float, optional
        Direct-wave velocity (m/s). Default is 1500.0.
    flat_len : float, optional
        Number of samples in the fully-muted (flat) region.
        If None, defaults to one wavelet period (1/freq/dt).
    taper_len : float, optional
        Number of samples in the taper (transition) region.
        If None, defaults to one wavelet period (1/freq/dt).

    Returns
    -------
    muted_trace : np.ndarray (float32)
        Trace with the direct arrival muted.
    """
    trace = np.asarray(trace, dtype=np.float32)
    nt = trace.shape[0]

    # dominant wavelet period, used as a default scale reference
    period = 1.0 / freq
    if flat_len is None:
        flat_len = period / dt
    if taper_len is None:
        taper_len = period / dt

    # direct-arrival time (linear moveout, "whole offset")
    t_arrival = peak_time + abs(offset) / velocity
    arrival_sample = t_arrival / dt

    # time grid in samples
    t_grid = np.arange(nt)
    abs_dist = np.abs(t_grid - arrival_sample)

    half_flat = flat_len / 2.0
    taper_start = half_flat
    taper_end = half_flat + taper_len

    mask = np.ones(nt, dtype=np.float32)

    # fully muted (flat) region
    mask[abs_dist <= taper_start] = 0.0

    # cosine taper region (0 -> 1)
    in_taper = (abs_dist > taper_start) & (abs_dist <= taper_end)
    normalized_dist = (abs_dist[in_taper] - taper_start) / taper_len
    mask[in_taper] = 0.5 * (1.0 - np.cos(normalized_dist * np.pi))

    return trace * mask

# ---------------------------------------------------------------------

def normalize_data(
    data,
    method="trace_max",
    axis=-1,
    eps=1e-10,
):
    """
    Normalizes seismic data using different strategies.

    Parameters
    ----------
    data : np.ndarray
        Seismic data. Can be:
        - 1D (nt,)                -> single trace
        - 2D (n_traces, nt)       -> shot gather
        - 3D (n_shots, n_rec, nt) -> full dataset
    method : str, optional
        Normalization method:
        - "trace_max" : each trace divided by its own max absolute
          amplitude (default). Preserves relative amplitude between
          samples within a trace, but not between traces.
        - "global_max" : entire array divided by a single global max
          absolute amplitude. Preserves relative amplitude between
          traces/shots.
        - "rms" : each trace divided by its own RMS value.
        - "global_rms" : entire array divided by the global RMS value.
    axis : int, optional
        Axis along which "trace_max"/"rms" are computed (default -1,
        i.e. the time axis, assuming trace shape (..., nt)).
    eps : float, optional
        Small value to avoid division by zero.

    Returns
    -------
    normalized_data : np.ndarray (float32)
        Normalized data, same shape as input.
    """
    data = np.asarray(data, dtype=np.float32)

    if method == "trace_max":
        norm_factor = np.max(np.abs(data), axis=axis, keepdims=True)
        norm_factor = np.maximum(norm_factor, eps)
        return data / norm_factor

    elif method == "global_max":
        norm_factor = np.max(np.abs(data))
        norm_factor = max(norm_factor, eps)
        return data / norm_factor

    elif method == "rms":
        norm_factor = np.sqrt(np.mean(data**2, axis=axis, keepdims=True))
        norm_factor = np.maximum(norm_factor, eps)
        return data / norm_factor

    elif method == "global_rms":
        norm_factor = np.sqrt(np.mean(data**2))
        norm_factor = max(norm_factor, eps)
        return data / norm_factor

    else:
        raise ValueError(
            f"Unknown method '{method}'. "
            "Choose from: 'trace_max', 'global_max', 'rms', 'global_rms'."
        )
    
# ---------------------------------------------------------------------

def crop_seismogram_by_depth(seismogram, base_model, dz, dt, z_max):
    """
    Crops a 1D seismogram based on a maximum physical depth (z_max),
    converting the limit to time using the reference velocity profile.

    Parameters
    ----------
    seismogram : ndarray
        The seismic trace (observed or synthetic) to be cropped.
    base_model : ndarray
        Reference velocity profile (Vp) in m/s.
    dz : float
        Depth sampling interval (m).
    dt : float
        Time sampling interval (s).
    z_max : float
        Maximum physical depth (m) used to define the cutoff.

    Returns
    -------
    cropped_seismogram : ndarray
        Seismic trace cropped at the TWT corresponding to z_max.
    """
    # 1. Limit the base model to the maximum depth for the calculation
    max_depth_idx = int(z_max / dz)
    vp_cut = base_model[:max_depth_idx]
    
    # 2. Compute the Two-Way Travel Time (TWT) up to z_max
    # Integration of slowness (1/v) along depth
    twt_max = 2.0 * np.sum(dz / vp_cut)
    
    # 3. Convert TWT to the corresponding index in the trace (sample)
    sample_max = int(twt_max / dt)
    
    # 4. Ensure the index does not exceed the length of the provided seismogram
    sample_max = min(sample_max, len(seismogram))
    
    # 5. Perform the crop
    return seismogram[:sample_max]

# ---------------------------------------------------------------------

def calculate_synthetic_trace(vp_profile, rho_profile, wavelet, dz, dt, nt):
    """
    1D convolutional seismic modeling with strict phase alignment.

    Computes acoustic impedance and reflection coefficients from the
    velocity/density profiles, converts the reflectivity series from
    depth to two-way travel time (TWT) using interval slowness, maps
    it onto a regular time grid, and convolves it with the source
    wavelet to produce a synthetic seismic trace.

    Parameters
    ----------
    vp_profile : ndarray
        P-wave velocity profile (m/s), sampled at depth interval dz.
    rho_profile : ndarray
        Density profile (kg/m^3 or g/cm^3), same length as vp_profile.
    wavelet : ndarray
        Source wavelet used for convolution.
    dz : float
        Depth sampling interval (m).
    dt : float
        Time sampling interval (s) of the output trace.
    nt : int
        Number of time samples in the output trace.

    Returns
    -------
    synthetic_trace : ndarray of shape (nt,)
        Synthetic seismic trace, as float32.

    Notes
    -----
    - Reflection coefficients (RC) are computed at each interface
      (n-1 values for n samples) using the normal-incidence formula.
    - TWT for each interface is obtained by cumulatively summing
      2*dz/v_avg (average velocity between consecutive samples),
      rather than by linear interpolation, to avoid time jitter.
    - Reflectivity is placed at the nearest time sample (via rounding
      and np.add.at) instead of being linearly interpolated, which
      preserves reflector amplitude and avoids smearing.
    - Convolution uses mode='same'; this assumes the wavelet is
      zero-phase (centered on zero). If the wavelet is not centered,
      the output may be shifted in time — consider using mode='full'
      and manually cropping around the center instead.
    """
    # 1. Impedance and Reflection Coefficient (in the depth domain)
    z = vp_profile * rho_profile
    rc = np.diff(z) / (z[:-1] + z[1:])  # RC computed at each interface (n-1)
    
    # 2. Conversion to Time (TWT) using slowness-based interpolation
    # dt_sample = (2 * dz) / v_i
    dt_samples = (2.0 * dz) / ((vp_profile[:-1] + vp_profile[1:]) / 2.0)
    twt_interface = np.cumsum(dt_samples)
    
    # 3. Direct mapping onto the time grid (avoids linear-interpolation jitter)
    # Create a zero vector and insert reflection energy at the nearest sample
    reflectivity_time = np.zeros(nt)
    indices = np.round(twt_interface / dt).astype(int)
    
    # Filter out-of-bounds indices for safety
    valid_mask = indices < nt
    np.add.at(reflectivity_time, indices[valid_mask], rc[valid_mask])
    
    # 4. Convolution with the wavelet (use 'same' only if the wavelet is zero-centered)
    # To avoid shifting, it's often better to use 'full' and crop the center
    synthetic_trace = np.convolve(reflectivity_time, wavelet, mode='same')
    
    return synthetic_trace.astype(np.float32)

# ---------------------------------------------------------------------
# 2. DEAP OBJECTIVE (FITNESS) FUNCTION
# ---------------------------------------------------------------------

def phase_misfit_objective(individual, base_model, eof_basis, s_obs, f_peak, dt, dz, nt, z_max, rho_model):
    """
    Fitness function for the genetic algorithm, based on maximizing the
    normalized cross-correlation (i.e., minimizing phase misfit) between
    the observed and synthetic seismic traces.

    Parameters
    ----------
    individual : list
        DEAP chromosome (candidate EOF expansion coefficients).
    base_model : ndarray
        Background (reference) P-wave velocity profile (vp0).
    eof_basis : ndarray
        Matrix containing the retained empirical orthogonal functions (EOFs).
    s_obs : ndarray
        Observed (or target synthetic) seismic trace.
    f_peak : float
        Central frequency of the source wavelet / low-pass filter in Hz
    dz : float
        Depth sampling interval of vp_profile/rho_profile (m).
    dt : float
        Time sampling interval of the desired output trace (s).
    nt : int
        Number of time samples of the desired output trace.
    z_max : int
        Maximum depth for the inversion grid in meters (defines the TWT window).
    rho_model : ndarray
        Background density profile (rho0), kept fixed during inversion.

    Returns
    -------
    tuple
        (misfit,) — DEAP requires fitness values to be returned as a tuple.
    """

    # 1. Reconstruct the Vp profile from the EOF basis and the individual's genes
    coefs = np.array(individual)
    vp_reconstructed = base_model + np.dot(eof_basis, coefs)

    # Optional: apply a severe penalty if the algorithm suggests physically
    # unrealistic water velocities (e.g., < 1400 m/s or > 1600 m/s), in
    # order to avoid numerical instability.
    if np.any(vp_reconstructed < 1400) or np.any(vp_reconstructed > 1600):
        return 10.0,  # Return a very poor fitness value

    wavelet = ricker_wavelet(f_peak, dt)

    # 2. Forward Modeling
    s_syn = calculate_synthetic_trace(vp_reconstructed, rho_model, wavelet, dz, dt, nt)

    # 3. Mute
    s_syn = mute_direct_arrival(trace=s_syn,offset=25,dt=dt,freq=f_peak,peak_time=1.5/f_peak)
    s_syn = bandpass_filter(s_syn, dt=dt, f_peak=f_peak) 

    # 4. Crop
    s_syn = crop_seismogram_by_depth(seismogram=s_syn, base_model=base_model, dz=dz, dt=dt, z_max=z_max)
    s_syn = apply_edge_taper(s_syn)

    # 5. Normalização consistente: 
    s_syn = normalize_data(s_syn, method="trace_max")

    # 6. NCC
    norm_obs = np.linalg.norm(s_obs)
    norm_syn = np.linalg.norm(s_syn)

    # Evita divisão por zero
    if norm_syn < 1e-6 or norm_obs < 1e-6:
        return 2.0,

    ncc = np.dot(s_obs, s_syn) / (norm_obs * norm_syn)
    mse = np.mean((s_obs - s_syn)**2)

    # Normalização Científica:
    # O termo de fase (1 - NCC) está naturalmente limitado entre 0 e 2.
    # O termo de MSE, após a normalização do traço, costuma ser bem menor.
    # Para igualar, escalamos o MSE pelo seu desvio padrão ou pelo valor médio esperado.
    # Uma forma simples de equilibrar sem parâmetros adicionais é:

    misfit_phase = (1.0 - ncc)
    misfit_amplitude = mse / np.mean(s_obs**2) # MSE relativo à energia do sinal observado

    # Agora ambos estão normalizados em termos de "energia relativa"
    misfit = 0.25 * misfit_phase + 0.75 * misfit_amplitude


    return misfit,