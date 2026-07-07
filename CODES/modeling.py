# Functions 

## modeling 

import numpy as np
import xarray as xr
from scipy.ndimage import gaussian_filter1d
import torch

import deepwave
from deepwave import elastic
import deepwave.common as dwc

import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np


from parameters_py.config import (
					FL_VPVS,LL_VPVS)

def calculate_parameters_from_vs(vel_s,vp_vs):
    
    '''
    Function to estimate the P velocities and Densities in function of initial Shear wave velocity and Vp/Vs ratio
    
    == Empirical Relations between Elastic Wavespeeds and Density == 
    
    Thomas M. Brocher; Empirical Relations between Elastic Wavespeeds and Density in the Earth's Crust. 
    Bulletin of the Seismological Society of America 2005;; 95 (6): 2081–2092. doi: https://doi.org/10.1785/0120050077
    
    The purpose of this paper is present empirical relations among Vp, Vs, and rho, that can be used to infer Vs for
    the entire Earth’s crust from either Vp or rho. Vp/Vs is dependent on many factors, including fluid content and 
    bulk chemical composition, and there is no direct relation between it and Vp.
    
    Brocher (2005) presented a series of empirical relationships between velocity and density:

    - rho = 1.6612*Vp - 0.4721*Vp**2 + 0.0671*Vp**3 - 0.0043*Vp**4 + 0.000106*Vp**5

    - here rho is in g/cm³, and Vp is in km/s (Nafe–Drake equation).    

    == Vp/Vs ratio == 

    Vp/Vs velocity ratio is a strong function of:
    - water saturation, 
    - porosity, 
    - crack intensity, and
    - clay content. 
    
    Recently, it became important factor to study underground properties. 
    Vp/Vs velocity ratio has the variation interval as 1.45 to 8 and 
    it have been used as a lithological indicators in studies of soil amplification 
    and soil classification, acquifers and hydrocarbon reservoirs.

    
    Parameters:
    ----------
    - vs – Velocities in m/s.    
    - vp/vs – Ratio between Shear and Compressional wave velocity.    
    
    Returns
    -------
    Array corresponding to P-wave velocity and density:
        - vp velocities (in m/s) 
        - density (in g/cm³)  
    ''' 
    
    vp = (vel_s*vp_vs)/1000

    rho = 1.6612*vp - 0.4721*vp**2 + 0.0671*vp**3 - 0.0043*vp**4 + 0.000106*vp**5

    return vp*1000,rho

# -----------------------------------------------------------------------------------------------------

def create_velocity_model_from_profile_vs(model_profile,FL_boundvpvs=FL_VPVS,LL_boundvpvs=LL_VPVS):
    """
    Constructs a 1D seismic velocity model from a profile of layer thicknesses, 
    shear-wave velocities (vs), and ratio between Shear and Compressional wave 
    velocity (vpvs) .

    This function calculates P-wave velocity (vp) and density from the provided 
    vs and vpvs values. It also converts the input units (m, m/s) 
    to standard modeling units (km, km/s). The final layer is treated as 
    an infinite half-space with a thickness of 0.

    Parameters
    ----------
    model_profile : tuple or list of array-like
        A sequence containing three arrays/lists: 
        - [0]: Layer thicknesses (in m).
        - [1]: Layer S-wave velocities (vs) (in m/s).

    Returns
    -------
    velocity_model : numpy.ndarray
        A 2D array of shape (N, 4) representing the 1-D velocity model. 
        Each row corresponds to a layer with the following columns:
        [thickness, velocity_p, velocity_s, density]
            - Layer thickness (in km).
            - Layer P-wave velocity (in km/s).
            - Layer S-wave velocity (in km/s).
            - Layer density (in g/cm³).
    """
    
    # Initialize an empty list to store the parameters for each layer
    vmodel = []

    # Initialize Vp/Vs ratio for each layer
    vpvs = np.linspace(start=FL_boundvpvs, stop=LL_boundvpvs, num=len(model_profile))

    # Unpack and iterate through the thicknesses, and vs values simultaneously.
    for i, (thickness, vs) in enumerate(model_profile):
 
        # Calculate P-wave velocity and density using an external helper function
        vp, dens = calculate_parameters_from_vs(vs,vpvs[i])
        
        # Check if the current layer is NOT the last layer in the profile
        if not i == len(model_profile) - 1:
            # Append standard layer: convert units from (m, m/s, m/s, g/m³) to (km, km/s, km/s, g/cm³)
            # by dividing all values by 1000
            vmodel.append([thickness / 1000, vp / 1000, vs / 1000, dens])
        else: 
            # Append final layer (half-space): set thickness to 0.
            # Convert vp, vs, and density units as done for the previous layers.
            vmodel.append([0.0, vp / 1000, vs / 1000, dens])

    # Convert the resulting list of layer parameters into a NumPy array
    velocity_model = np.array(vmodel)    

    return velocity_model

# -----------------------------------------------------------------------------------------------------

def create_seismic_model(depth_ranges,vs_ranges,vpvs_ranges,num_layers,max_total,DX,DZ,NX,NZ,basement_vs,basement_vpvs,seed):
    """
    Generates a 2D synthetic seismic velocity model with discrete, blocky lateral 
    variations in layer depth. 

    The model consists of horizontal segments where the depth of interfaces remains 
    constant within a segment but jumps abruptly at segment boundaries. The physical 
    properties (Vs, Vp, density) vary vertically by layer but remain constant 
    laterally within each layer.

    Parameters
    ----------
    depth_ranges : list of tuple of float
        Min and max depth boundaries for each layer interface [(min, max), ...].
    vs_ranges : list of tuple of float
        Min and max S-wave velocities (m/s) for each layer.
    vpvs_ranges : list of tuple of float
        Min and max Vp/Vs ratios for each layer.
    num_layers : int
        Number of geological layers above the basement.
    max_total : float
        Maximum allowable depth (m) for any interface.
    DX : float
        Grid spacing in the horizontal (X) direction (m).
    DZ : float
        Grid spacing in the vertical (Z) direction (m).
    NX : int
        Number of grid points in the horizontal direction.
    NZ : int
        Number of grid points in the vertical direction.
    basement_vs : float
        S-wave velocity (m/s) of the basement rock.
    basement_vpvs : float
        Vp/Vs ratio of the basement rock.
    seed : int
        Seed for the random number generator to ensure reproducibility.

    Returns
    -------
    ds : xarray.Dataset
        A dataset containing the 2D grids ('z', 'x') for Vp, Vs, Vp/Vs ratio, 
        density, and formation indices, along with spatial coordinates and attributes.
    """

    # ---------------------------------------------------------
    # 1. Initialization and Grid Setup
    # ---------------------------------------------------------
    
    # Initialize the random number generator for reproducible results
    rng = np.random.default_rng(seed)

    # Create spatial coordinate vectors based on grid spacing and size
    z_vector = np.linspace(0, NZ * DZ, NZ)
    x_vector = np.linspace(0, NX * DX, NX)

    # ---------------------------------------------------------
    # 2. Interface (Topography) Generation
    # ---------------------------------------------------------
    
    # Initialize the interfaces array: shape (num_layers, NX).
    # This stores the depth of the bottom of each layer across all columns.
    interfaces = np.zeros((num_layers, NX))

    # Define the number of discrete lateral segments (e.g., 4 for quarters)
    num_segments = 10

    # Split the column indices into the specified number of segments.
    # np.array_split ensures all columns are covered even if NX isn't perfectly divisible.
    segment_indices = np.array_split(np.arange(NX), num_segments)

    for l in range(num_layers):
        for segment_cols in segment_indices:
            # Generate ONE random depth for this specific horizontal segment
            segment_depth = rng.uniform(*depth_ranges[l])
            
            # Apply that uniform depth to all columns within the current segment
            interfaces[l, segment_cols] = segment_depth
            
        # Prevent non-physical layer crossings: 
        # A layer's bottom must be strictly deeper than the layer above it.
        if l > 0:
            interfaces[l] = np.maximum(interfaces[l], interfaces[l - 1] + 0.01)
            
        # Ensure no interface exceeds the predefined maximum model depth
        interfaces[l] = np.clip(interfaces[l], 0, max_total)
    
    # ---------------------------------------------------------
    # 3. Layer Property Generation
    # ---------------------------------------------------------
    
    # Generate constant physical properties for each layer (no lateral variation)
    vs_layer   = np.array([round(rng.uniform(*vs_ranges[l])) for l in range(num_layers)])
    vpvs_layer = np.array([round(rng.uniform(*vpvs_ranges[l]), 2) for l in range(num_layers)])
    
    # Derive Vp and Density from Vs and Vp/Vs using the external function
    vp_layer, density_layer = calculate_parameters_from_vs(vs_layer, vpvs_layer)

    # Calculate basement parameters
    basement_vp, basement_density = calculate_parameters_from_vs(
        np.array([basement_vs]), np.array([basement_vpvs])
    )
    
    # ---------------------------------------------------------
    # 4. 2D Grid Population
    # ---------------------------------------------------------
    
    # Initialize 2D arrays with basement values as the default background
    vs_grid        = np.full((NZ, NX), basement_vs)
    vp_grid        = np.full((NZ, NX), basement_vp[0])
    vpvs_grid      = np.full((NZ, NX), basement_vpvs)
    density_grid   = np.full((NZ, NX), basement_density[0])
    formation_grid = np.zeros((NZ, NX), dtype=int) # 0 = basement, 1..N = layers

    # Iterate through every column (X) and every depth point (Z) to assign values
    for ix in range(NX):
        for iz, z in enumerate(z_vector):
            
            # Optimization: If current Z is deeper than the last interface, 
            # we are in the basement. Since arrays were initialized with 
            # basement values, we can break early to save computation time.
            if z >= interfaces[-1, ix]:      
                break                        
            
            # Determine which layer the current (Z, X) point falls into
            for l in range(num_layers):
                top    = interfaces[l - 1, ix] if l > 0 else 0.0
                bottom = interfaces[l, ix]
                
                if top <= z < bottom:
                    # Assign the corresponding layer properties
                    vs_grid[iz, ix]        = vs_layer[l]
                    vp_grid[iz, ix]        = vp_layer[l]
                    vpvs_grid[iz, ix]      = vpvs_layer[l]
                    density_grid[iz, ix]   = density_layer[l]
                    formation_grid[iz, ix] = l + 1
                    break # Move to the next depth point once assigned

    # ---------------------------------------------------------
    # 5. Data Packaging
    # ---------------------------------------------------------
    
    # Define dimensions and coordinates for the xarray dataset
    dims   = ("z", "x")
    coords = {
        "distance": ("x", x_vector, {"units": "m", "long_name": "distance"}),
        "depth": ("z", z_vector, {"units": "m", "long_name": "depth"})
    }

    # Pack the numpy grids into an xarray.Dataset with rich metadata
    ds = xr.Dataset(
        {
            "vp":        xr.DataArray(vp_grid,        dims=dims, coords=coords,
                                      attrs={"units": "m/s",   "long_name": "P-wave velocity"}),
            "vs":        xr.DataArray(vs_grid,        dims=dims, coords=coords,
                                      attrs={"units": "m/s",   "long_name": "S-wave velocity"}),
            "vpvs":      xr.DataArray(vpvs_grid,      dims=dims, coords=coords,
                                      attrs={"units": "-",     "long_name": "Vp/Vs ratio"}),
            "density":   xr.DataArray(density_grid,   dims=dims, coords=coords,
                                      attrs={"units": "g/cm³", "long_name": "Bulk density"}),
            "formation": xr.DataArray(formation_grid, dims=dims, coords=coords,
                                      attrs={"long_name": "Formation index (0=basement, 1..N=soil)"}),
        },
        attrs={
            "DX": DX, 
            "DZ": DZ,
            "NX": NX, 
            "NZ": NZ,
            "num_layers":             num_layers,
            "vs_per_layer_m_s":       vs_layer.tolist(),
            "vp_per_layer_m_s":       vp_layer.tolist(),
            "density_per_layer_gcm3": density_layer.tolist(),
            "basement_vp_m_s":        basement_vp,
            "basement_vs_m_s":        basement_vs,
            "basement_vpvs":          basement_vpvs,
            "basement_density_g_cm3": basement_density,
        },
    )
    
    return ds

# -----------------------------------------------------------------------------------------------------

def create_acquisition_geometry(n_shots,n_receivers_per_shot,d_source,d_receivers,first_source,first_receiver,source_depth,receiver_depth,n_shots_per_shot,device,plot_profile,data_to_plot,figures_path):
        
    """
    Creates the source and receiver location tensors for Deepwave forward modeling,
    simulating a roll-along acquisition geometry where the receiver array moves 
    with each shot.

    Parameters
    ----------
    n_shots : int
        Total number of independent seismic shots.
    n_receivers_per_shot : int
        Number of active receivers recording each shot.
    d_source : int
        Horizontal spacing between consecutive shot locations (in grid points).
    d_receivers : int
        Horizontal spacing between adjacent receivers (in grid points).
    first_source : int
        Horizontal grid index for the first source location.
    first_receiver : int
        Horizontal grid index for the first receiver (relative to the first shot).
    source_depth : int
        Vertical grid index for the sources (default is 0 for surface).
    receiver_depth : int
        Vertical grid index for the receivers (default is 0 for surface).
    n_shots_per_shot : int
        Number of sources activated simultaneously per shot (default is 1).
    device : torch.device
        The PyTorch device to allocate the tensors on.
    plot_profile : Boolean    
        To plot the seismic profile.
    data_to_plot : xarray
        XARRAY with the 2D soil model.
    figures_path : str
        Directory for saving figures.

    Returns
    -------
    source_locations : torch.Tensor
        Tensor of shape (n_shots, n_shots_per_shot, 2) containing source [z, x] indices.
    receiver_locations : torch.Tensor
        Tensor of shape (n_shots, n_receivers_per_shot, 2) containing receiver [z, x] indices.
    """

    # =========================
    # Source locations
    # =========================

    source_locations = torch.zeros(
        n_shots,
        n_shots_per_shot,
        2,
        dtype=torch.long,
        device=device
    )

    # Set z-coordinates (depth)
    source_locations[..., 0] = source_depth

    # Set x-coordinates
    source_locations[:, 0, 1] = (
        torch.arange(n_shots, device=device) * d_source
        + first_source
    )

    # =========================
    # Receiver locations
    # =========================
    receiver_locations = torch.zeros(
        n_shots,
        n_receivers_per_shot,
        2,
        dtype=torch.long,
        device=device
    )

    # Set z-coordinates (depth)
    receiver_locations[..., 0] = receiver_depth

    # 1. Base line of receivers for the FIRST shot
    base_receivers = (
        torch.arange(n_receivers_per_shot, device=device) 
        * d_receivers 
        + first_receiver
    )

    # 2. Offset applied to the receiver array for subsequent shots
    # Shape transformed to (n_shots, 1) for correct broadcasting
    shot_offsets = (torch.arange(n_shots, device=device) * d_source).unsqueeze(1)

    # 3. Add the base geometry to the offsets
    receiver_locations[:, :, 1] = base_receivers + shot_offsets

    if plot_profile:

        fig, ax = plt.subplots(1, 1, figsize=(20, 5))

        # -----
        # Model
        # -----

        im = ax.imshow(
            data_to_plot.density.data,
            cmap='viridis',
            extent=[0, data_to_plot.attrs['NX'] * data_to_plot.attrs['DX'], -data_to_plot.attrs['NZ'] * data_to_plot.attrs['DZ'], 0],
            vmin=data_to_plot.density.data.min(),
            vmax=data_to_plot.density.data.max(),
            aspect="auto",
            origin='upper'
        )

        # -------
        # Sources
        # -------

        src_x = source_locations[:, 0, 1].cpu().numpy() * data_to_plot.attrs['DX']
        src_z = source_locations[:, 0, 0].cpu().numpy() * data_to_plot.attrs['DZ']

        # ---------
        # Receivers
        # ---------

        rec_x = receiver_locations[:, :, 1].cpu().numpy() * data_to_plot.attrs['DX']
        rec_z = receiver_locations[:, :, 0].cpu().numpy() * data_to_plot.attrs['DZ']

        # ====
        # Plot
        # ====

        # Sources
        ax.scatter(
            src_x,
            -src_z,
            c='red',
            marker='*',
            s=200,
            edgecolors='k',
            label='Sources',
            zorder=1
        )

        # Receivers
        ax.scatter(
            rec_x,
            -rec_z,
            c='k',
            marker='v',
            s=50,
            label='Receivers',
            zorder=-10
        )

        ax.hlines(y=0,xmin=0,xmax=data_to_plot.attrs['NX'] * data_to_plot.attrs['DX'],colors='k',lw=2,ls='-',alpha=1)

        # labels
        ax.set_xlabel('Distance (m)')
        ax.set_ylabel('Depth (m)')

        # grid
        ax.grid(which='major', color='gray', linestyle='-', linewidth=0.5)
        ax.grid(which='minor', color='gray', linestyle=':', linewidth=0.5, alpha=0.5)

        ax.xaxis.set_ticks_position('both')
        ax.yaxis.set_ticks_position('both')

        ax.tick_params(
            which='both',
            direction='in',
            top=True,
            right=True
        )

        ax.set_xlim(0, data_to_plot.attrs['NX'] * data_to_plot.attrs['DX'])
        ax.set_ylim(-(data_to_plot.attrs['NZ'] * data_to_plot.attrs['DZ']), 0.15)

        # legend
        ax.legend()

        # colorbar
        plt.colorbar(
            im,
            ax=ax,
            fraction=0.15,
            shrink=0.5,
            label='Density (kg/cm³)'
        )

        fig.savefig(figures_path + 'model_slice_and_receptors.png')

    return source_locations, receiver_locations

# -----------------------------------------------------------------------------------------------------

def create_source_amplitudes(freq,nt,dt,n_shots,n_shots_per_shot,device,plot_wavelet,figures_path):
    
    """
    Generates a Ricker wavelet and formats the source amplitudes tensor 
    required for Deepwave forward modeling.

    Parameters
    ----------
    freq : float
        Central frequency of the Ricker wavelet in Hertz.
    nt : int
        Total number of discrete time steps.
    dt : float
        Time step duration in seconds.
    n_shots : int
        Total number of independent seismic shots.
    n_shots_per_shot : int
        Number of sources activated simultaneously per shot.
    device : torch.device
        The PyTorch device to allocate the tensors on.
    plot_wavelet : bool
        Flag to plot the generated wavelet.
    figures_path : str
        Directory path for saving the figure (e.g., 'outputs/figures/').

    Returns
    -------
    source_amplitudes : torch.Tensor
        Tensor of shape (n_shots, n_shots_per_shot, nt) containing the 
        source wavelet amplitudes over time.
    """

    # Calculate the peak time based on the central frequency
    peak_time = 1 / freq

    # Generate the 1D Ricker wavelet using Deepwave
    wavelet = deepwave.wavelets.ricker(freq, nt, dt, peak_time)

    # Reshape and repeat the wavelet for all shots and sources per shot
    source_amplitudes = (wavelet.repeat(n_shots, n_shots_per_shot, 1).to(device=device, dtype=torch.float32))

    # =========================
    # Optional Plotting
    # =========================
    if plot_wavelet:
        fig, ax = plt.subplots(1, 1, figsize=(17, 3))
        
        t_axis = np.arange(nt) * dt
        
        ax.plot(t_axis, wavelet.numpy(), color='black', lw=1.5)
        
        ax.set_xlabel('Time [s]')
        ax.set_ylabel('Amplitude')
        ax.set_title(f'Ricker wavelet  f={freq} Hz,  peak at t0={peak_time:.4f} s')
        
        # Grid and ticks formatting for consistency
        ax.grid(which='major', color='gray', linestyle='-', linewidth=0.5)
        ax.grid(which='minor', color='gray', linestyle=':', linewidth=0.5, alpha=0.5)
        ax.tick_params(which='both', direction='in', top=True, right=True)
        
        plt.tight_layout()
        
        # Save the figure using the provided path
        fig.savefig(figures_path + 'ricker_wavelet.png', dpi=300, bbox_inches='tight')
        plt.show()

    return source_amplitudes

# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------

def simulate_elastic_shots(vp,vs,rho,source_amplitudes,source_locations,receiver_locations,dx,dz,dt,pml_width,pml_freq,accuracy,device,plot_shot_gathers,figures_path):

    """
    Runs Deepwave's elastic forward modeling to generate synthetic shot gathers.
    
    Converts P-wave velocity, S-wave velocity, and density into Lamé parameters 
    and buoyancy before running the wave propagator. Extracts the vertical 
    particle velocity component simulating vertical geophones.
    
    Parameters
    ----------
    vp : torch.Tensor
        P-wave velocity model.
    vs : torch.Tensor
        S-wave velocity model.
    rho : torch.Tensor
        Density model.
    source_amplitudes : torch.Tensor
        The source wavelet amplitudes over time.
    source_locations : torch.Tensor
        Tensor containing [z, x] grid indices for the sources.
    receiver_locations : torch.Tensor
        Tensor containing [z, x] grid indices for the receivers.
    dx : float
        Horizontal grid cell size in meters.
    dz : float
        Vertical grid cell size in meters.
    dt : float
        Time step duration in seconds.
    pml_width : list
        Number of absorbing boundary layers [top, bottom, left, right].
    pml_freq : float
        Dominant frequency of the source for PML tuning.
    accuracy : int
        Finite difference spatial accuracy order.
    device : torch.device
        PyTorch device to run the computation on.
    plot_shot_gathers : bool
        If True, plots a Common Shot Gather (CSG) for the first shot.
    figures_path : str
        Directory path for saving the generated figure.
        
    Returns
    -------
    numpy.ndarray
        The recorded vertical particle velocities (synthetic shot gathers) 
        as a numpy array.
    """
    
    dtype = torch.float32
    
    # Ensure all inputs are on the correct device and use float32
    vp = vp.to(device=device, dtype=dtype)
    vs = vs.to(device=device, dtype=dtype)
    rho = rho.to(device=device, dtype=dtype)
    source_amplitudes = source_amplitudes.to(device=device, dtype=dtype)
    source_locations = source_locations.to(device)
    receiver_locations = receiver_locations.to(device)

    # Convert to Lamé parameters and buoyancy
    lam, mu, buoyancy = dwc.vpvsrho_to_lambmubuoyancy(vp, vs, rho)

    # Run Elastic Propagator
    out = elastic(
        lam, 
        mu, 
        buoyancy,
        grid_spacing=[dz, dx],
        dt=dt,
        source_amplitudes_y=source_amplitudes,
        source_locations_y=source_locations,
        receiver_locations_y=receiver_locations,
        pml_width=pml_width,
        pml_freq=pml_freq,
        accuracy=accuracy
    )
    
    # out[-2] corresponds to receiver_amplitudes_y (vertical component)
    shot_gathers = out[-2].cpu().numpy()
    
    # ====================
    # Plotting shot gather
    # ====================
    
    if plot_shot_gathers:
        # We plot the first shot (index 0) as a representative CSG
        shot_idx = 0
        nt = shot_gathers.shape[-1]
        t = np.arange(nt) * dt
        
        # Convert grid indices to physical distances (meters)
        src_x = source_locations[shot_idx, 0, 1].cpu().item() * dx
        rec_x = receiver_locations[shot_idx, :, 1].cpu().numpy() * dx
        
        # Determine receiver spacing to scale wiggles so they don't overlap too much
        rec_spacing = np.abs(rec_x[1] - rec_x[0]) if len(rec_x) > 1 else dx
        wiggle_scaling = rec_spacing * 0.8  # Max amplitude spans 80% of the spacing
        
        fig, ax = plt.subplots(figsize=(20, 5))
        
        # Plot individual traces
        for i in range(len(rec_x)):
            trace = shot_gathers[shot_idx, i, :]
            max_amp = np.max(np.abs(trace))
            
            # Normalize trace and apply dynamic scaling
            if max_amp > 1e-12:
                trace_norm = (trace / max_amp) * wiggle_scaling
            else:
                trace_norm = trace
            
            # Plot trace centered at its physical receiver X-coordinate
            ax.plot(trace_norm + rec_x[i], t, c='k', lw=0.5)

        # Plot Source and Receivers at the top (t=0 line)
        ax.scatter(src_x, 0, c='red', marker='*', s=250, edgecolors='k', label='Source', zorder=10)
        ax.scatter(rec_x, np.zeros_like(rec_x), c='blue', marker='v', s=40, edgecolors='k', label='Receivers', zorder=9)

        # ---------------------------------------------------------
        # NEW: Draw Offset Dimension Line
        # ---------------------------------------------------------
        # Place the line down slightly in time (e.g., 8% of the total time window)
        # so it sits between the surface markers and the wave arrivals.
        offset_y_pos = dt * nt * 0.08 
        offset_val = rec_x[0] - src_x
        
        # The |-| arrow style creates a dimension line
        ax.annotate(
            '', 
            xy=(src_x, offset_y_pos), 
            xytext=(rec_x[0], offset_y_pos),
            arrowprops=dict(arrowstyle='|-|', color='black', lw=1.5, shrinkA=0, shrinkB=0)
        )
        
        # Add the text label exactly in the middle of the line
        ax.text(
            (src_x + rec_x[0]) / 2, 
            offset_y_pos - (dt * nt * 0.015), # Slightly above the line (remember Y is inverted)
            f'Offset = {offset_val:.1f} m',
            ha='center', 
            va='bottom', 
            fontsize=12, 
            fontweight='bold', 
            color='black',
            bbox=dict(facecolor='white', edgecolor='none', pad=2.0, alpha=0.9) # White background so traces don't clash with text
        )
        # ---------------------------------------------------------

        ax.set_ylabel("Time (s)", fontsize=12)
        ax.set_xlabel("Distance (m)", fontsize=12)
        
        # Invert Y-axis so time goes downward
        ax.set_ylim(dt * nt, -dt * (nt * 0.05))  # Slight negative buffer to show the markers
        
        # Set X-axis limits with a buffer
        ax.set_xlim(min(src_x, rec_x.min()) - 2*rec_spacing, max(src_x, rec_x.max()) + 2*rec_spacing)
        
        ax.set_title(f"Common Shot Gather (CSG) - Shot {shot_idx + 1}", fontsize=14, fontweight='bold')
        ax.legend(loc='lower left') # Moved legend so it doesn't overlap the new offset text
        
        ax.grid(which='major', color='gray', linestyle='-', linewidth=0.5, alpha=0.5)
        
        plt.tight_layout()
        
        if figures_path:
            fig.savefig(figures_path+'common_shot_gather.png', dpi=300, bbox_inches='tight')
            
        plt.show()
    
    return shot_gathers

# -----------------------------------------------------------------------------------------------------------------------

def plot_seismic_model(ds, variables=("density", "vs", "vp", "vpvs"),receiver_x=None, figsize=(14, 10),path_output=None):
    """
    Plot seismic model variables from an xarray.Dataset.

    ds          : xarray.Dataset with dims (z, x) and coords 'distance', 'depth'
    variables   : tuple of variable names to plot (one subplot each)
    receiver_x  : array of receiver x-positions (m) for triangle markers; 
                  if None, auto-spaced every 10 m
    figsize     : figure size
    output      : path to save the figure
    """

    # Coordinates
    x = ds["distance"].values
    z = ds["depth"].values

    if receiver_x is None:
        receiver_x = np.arange(0, x.max() + 1, 10)

    n = len(variables)
    fig, axes = plt.subplots(n, 1, figsize=figsize, sharex=True)
    if n == 1:
        axes = [axes]

    # Colormaps per variable
    cmaps = {"density": "YlGnBu_r", "vs": "RdYlBu", "vp": "RdYlBu", "vpvs": "PuOr"}
    labels = {
        "density": "Density (g/cm³)",
        "vs":      "Vs (m/s)",
        "vp":      "Vp (m/s)",
        "vpvs":    "Vp/Vs",
    }

    for ax, var in zip(axes, variables):
        data = ds[var].values
        cmap = cmaps.get(var, "viridis")

        pcm = ax.pcolormesh(
            x, -z, data,
            cmap=cmap, shading="auto",
        )
        cbar = fig.colorbar(pcm, ax=ax, pad=0.01, fraction=0.02)
        cbar.set_label(labels.get(var, var), fontsize=10)

        # dashed vertical lines at receivers
        for rx in receiver_x:
            ax.axvline(rx, color="black", linestyle="--", linewidth=0.6, alpha=0.5)

        ax.set_ylabel("Depth (m)", fontsize=10)
        ax.set_ylim(-z.max(), 0)
        ax.set_xlim(x.min(), x.max())
        ax.set_title(labels.get(var, var), fontsize=11, loc="left")

    axes[-1].set_xlabel("Distance (m)", fontsize=11)
    fig.suptitle("Seismic Model", fontsize=13)
    fig.savefig(path_output)
    plt.tight_layout()
    plt.show()

# -----------------------------------------------------------------------------------------------------
