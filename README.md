# 🧬 Inversion of NFH data Using Evolutionary Algorithms 🧬

This repository contains reproducible material for the study *"...."* by **Diogo Luiz de Oliveira Coelho, and ....**, submitted to *Applied Computing and Geosciences*.

The provided scripts and notebooks demonstrate the generation and inversion of ... data, enabling the retrieval of the Sound Speed Profile as a function of depth using Evolutionary Algorithms.

## 📦 Required Libraries 📦

The following libraries are used in this project:

- [NumPy (2.2.6)](https://numpy.org/): Fundamental package for numerical computing in Python.
- [Pandas (2.3.1)](https://pandas.pydata.org/): Data analysis and manipulation tool.
- [Matplotlib (3.10.0)](https://matplotlib.org/): Visualization library for creating static, animated, and interactive plots.
- [tqdm (4.67.1)](https://github.com/tqdm/tqdm): Library for displaying progress bars in loops and scripts.
- [Verde (1.8.1)](https://www.fatiando.org/verde/latest/): Processing and interpolation of spatial data.
- [GemPy (2025.1.0)](https://www.gempy.org/): 3D structural geological modeling library.
- [GemGIS (1.1.8)](https://github.com/cgre-aachen/gemgis): Geospatial processing library.
- [PyVista (0.45.2)](https://pyvista.org/): 3D visualization and mesh analysis library.
- [DEAP (1.4.3)](https://deap.readthedocs.io/): Evolutionary algorithm framework for optimization tasks.
- [Disba (0.7.0)](https://github.com/keurfonluu/disba): Surface wave dispersion analysis.

## 📀 Installation 📀

This project provides an `environment.yml` file to ensure reproducibility of results.  
To set up the environment, install [Conda](https://docs.conda.io/) and run:

```bash
conda env create -f environment.yml
conda activate Evolutionary_inversion_soil


## 🏗️ Project structure 🏗️
This repository is organized as follows:

* 🗃️ **CODES**: These scripts collectively enable the generation, simulation, inversion, and visualization of surface wave dispersion data. 🚀
    * 🗒️ **dispersion_curves.py**: Computes **surface wave dispersion curves** based on velocity profiles, essential for inversion analysis.
    * 🗒️ **evolutionary_algorithm.py**: Implements an **evolutionary algorithm** for inverting dispersion curves and estimating S-wave velocity profiles.  
    * 🗒️ **modeling.py**: Defines **geological models** and manages layer properties used in the inversion process. 
    * 🗒️ **pyvista_create_gif.py**: Uses **PyVista** to generate **3D visualizations** and create animated GIFs of the geological models. 

* 🗃️ **NOTEBOOKS**: set of jupyter notebooks reproducing the experiments in the paper (see below for more details);

## 📑 Notebooks 📑
The following notebooks are provided:

- 📔 ``creating_model_with_gempy.ipynb``: notebook used to generate the model
- 📔 ``generating_observed_data_via_disba.ipynb``: notebook performing the simulation on synthetic dataset
- 📔 ``inversion_evolutionary_algorithm.ipynb``: notebook to perform inversion using evolutionary algorithms
    * 🗒️ **pyvista_create_gif.py**: Uses **PyVista** to generate **3D visualizations** and create animated GIFs of the geological models. 

## 🖱️ Usage 🖱️

1. Clone this repository:
   ```bash
   git clone <repository_url>
   cd Evolutionary_inversion_soil
   ```
2. Open the Jupyter Notebook environment:
   ```bash
   jupyter-lab
   ```
3. Run the following notebooks to reproduce the results:
   - `creating_model_with_gempy.ipynb`: Generates a 3D geological model.
   - `generating_observed_data_via_disba.ipynb`: Simulates observed data (2D profile).
   - `inversion_evolutionary_algorithm.ipynb`: Performs the inversion using an evolutionary algorithm to retrieve the 2D velocity profile.

## 📺 Models utilized in Inversion 📺

- #### **MWC Model:** *No compaction*

Represents the reference model with normal density variation across layers.

<img src="Models/MWC.gif" width="400" align="center">

- #### **MCP1 Model:** *Compaction in the first layer*

Characterized by a higher density in the first layer.

<img src="Models/MCP1.gif" width="400" align="center">

- #### **MCP2_1 Model:** *Compaction in the second layer (with lateral variation)*

Features lateral density variations in the second layer.

<img src="Models/MCP2_1.gif" width="400" align="center">

- #### **MCP2_2 Model:** *Compaction in the second layer (without lateral variation)*

Uniform density in the second layer with no lateral variation.

<img src="Models/MCP2_2.gif" width="400" align="center">

#### These animations illustrate the 3D soil models used in the inversions.

## 📝 License 📝 

This project is licensed under the BSD-3 License. See the `LICENSE` file for details.

## 📚 References 📚  

The implementation of the algorithms and methods in this repository is based on the following key references:  

- Gallagher, K., & Sambridge, M. (1994). **Genetic algorithms: a powerful tool for large-scale nonlinear optimization problems**. *Computers & Geosciences*, 20(7–8), 1229–1236.  
- Fortin, F. A., Rainville, F. M., Gardner, M., Parizeau, M., & Gagné, C. (2012). **DEAP: Evolutionary Algorithms Made Easy**. *Journal of Machine Learning Research*, 13, 2171-2175.  
- Haskell, N. A. (1953). **The dispersion of surface waves on multi-layered media**. *Bulletin of the Seismological Society of America*, 43, 17-34.  
- Xia, J., Miller, R. D., & Park, C. B. (1999). **Estimation of near-surface shear-wave velocity by inversion of Rayleigh waves**. *Geophysics*, 64(3), 691-700.  
- Yamanaka, H., & Ishida, H. (1996). **Application of genetic algorithms to an inversion of surface-wave dispersion data**. *Bulletin of the Seismological Society of America*, 86, 436–444.


## 🔖 Disclaimer 🔖  

All experiments were conducted on two different setups running **Debian GNU/Linux 12 (Bookworm)**:  

- 💻 **AMD Ryzen 7 5700U** with **10 GB RAM**  
- 💻 **Intel® Core™ Ultra 9** with **64 GB RAM**  
- 💻 **13th Gen Intel® Core™ i9-13900KF** with **128 GB RAM**  

📣 **Multiprocessing is implemented.**  

---
For further details, refer to the paper associated with this repository.

