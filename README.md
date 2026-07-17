# 🌊 Inversion of NFH Data Using Evolutionary Algorithms 🌊

This repository contains reproducible material for the study *"[PAPER TITLE — TO CONFIRM]"* by **Diogo Luiz de Oliveira Coelho, and [...]**, submitted to *[journal — TO CONFIRM, e.g. Applied Computing and Geosciences]*.

The provided scripts and notebooks demonstrate the generation and inversion of seismic oceanography (NFH) data, combining synthetic seismic modeling via **Deepwave**, EOF-based analysis of **GLORYS** ocean reanalysis data, and inversion via **Evolutionary Algorithms**, enabling the retrieval of ocean sound speed / thermohaline structure as a function of depth.


## 📦 Required Libraries 📦

The libraries below were inferred from the project structure (use of Deepwave, GLORYS NetCDF data, EOF analysis, and evolutionary algorithms). **Please check and adjust against the actual `environment.yml`:**

- [NumPy](https://numpy.org/): Fundamental package for numerical computing in Python.
- [Pandas](https://pandas.pydata.org/): Data analysis and manipulation tool.
- [Matplotlib](https://matplotlib.org/): Visualization library for creating static, animated, and interactive plots.
- [tqdm](https://github.com/tqdm/tqdm): Library for displaying progress bars in loops and scripts.
- [xarray](https://docs.xarray.dev/): Handling of multi-dimensional labeled arrays, used for GLORYS NetCDF data.
- [Deepwave](https://github.com/ar4/deepwave): PyTorch-based library for seismic wave propagation and full-waveform modeling, used to generate observed data.
- [PyTorch](https://pytorch.org/): Deep learning framework, backend required by Deepwave.
- [eofs](https://ajdawson.github.io/eofs/) *(or equivalent)*: Empirical Orthogonal Function (EOF) analysis of ocean reanalysis fields.
- [DEAP](https://deap.readthedocs.io/): Evolutionary algorithm framework for optimization tasks.


## 📀 Installation 📀

This project provides an `environment.yml` file to ensure reproducibility of results.
To set up the environment, install [Conda](https://docs.conda.io/) and run:

```bash
conda env create -f environment.yml
conda activate NFH_inversion_seismic_oceanography
```

## 🏗️ Project structure 🏗️

This repository is organized as follows:

```
.
├── config_file.cnf
├── environment.yml
├── LICENSE.txt
├── notebooks
│   ├── 01_Generating_observed_data_via_deepwave.ipynb
│   ├── 02_EOF_glorys_analysis.ipynb
│   ├── 03_Inversion_evolutionary_algorithm.ipynb
│   └── 04_Results_analysis_inversion.ipynb
├── ocean_inversion
│   ├── config.py
│   ├── evolutionary_algorithm.py
│   ├── modeling.py
│   └── __init__.py
└── README.md
```

* 🗃️ **ocean_inversion**: Python package containing the source code that supports the entire generation, modeling, and inversion workflow. 🚀
    * 🗒️ **config.py**: Centralizes the experiment's configuration parameters (read from `config_file.cnf`).
    * 🗒️ **evolutionary_algorithm.py**: Implements the **evolutionary algorithm** used to invert the observed data and estimate the velocity/medium property profile.
    * 🗒️ **modeling.py**: Defines the **models** and manages the layer/structure properties used in the inversion process.

* 🗃️ **notebooks**: set of Jupyter notebooks reproducing the experiments in the paper (see details below).

* 📄 **config_file.cnf**: external configuration file read by `ocean_inversion/config.py`, defining the experiment's parameters (grid, layers, evolutionary algorithm hyperparameters, etc.). **[TO CONFIRM description]**

## 📑 Notebooks 📑

The following notebooks are provided, numbered in the order they should be executed:

- 📔 **`01_Generating_observed_data_via_deepwave.ipynb`**: generates the synthetic observed data through wave propagation modeling with **Deepwave**.
- 📔 **`02_EOF_glorys_analysis.ipynb`**: performs **Empirical Orthogonal Function (EOF)** analysis on **GLORYS** ocean reanalysis data.
- 📔 **`03_Inversion_evolutionary_algorithm.ipynb`**: performs the **inversion** of the observed data using the evolutionary algorithm implemented in `ocean_inversion/evolutionary_algorithm.py`.
- 📔 **`04_Results_analysis_inversion.ipynb`**: performs the **analysis and visualization of the inversion results**.

## 🖱️ Usage 🖱️

1. Clone this repository:
   ```bash
   git clone <repository_url>
   cd NFH_inversion_seismic_oceanography
   ```
2. Open the Jupyter Notebook environment:
   ```bash
   jupyter-lab
   ```
3. Run the following notebooks, in order, to reproduce the results:
   - `01_Generating_observed_data_via_deepwave.ipynb`: Generates synthetic observed data via seismic wave modeling with Deepwave.
   - `02_EOF_glorys_analysis.ipynb`: Performs EOF analysis on GLORYS ocean reanalysis data.
   - `03_Inversion_evolutionary_algorithm.ipynb`: Performs the inversion using an evolutionary algorithm to retrieve the target profile.
   - `04_Results_analysis_inversion.ipynb`: Analyzes and visualizes the inversion results.

## 📝 License 📝

This project is licensed under the BSD-3 License. See the `LICENSE.txt` file for details.

## 📚 References 📚

The implementation of the algorithms and methods in this repository is based on the following key references:

- Gallagher, K., & Sambridge, M. (1994). **Genetic algorithms: a powerful tool for large-scale nonlinear optimization problems**. *Computers & Geosciences*, 20(7–8), 1229–1236.
- Fortin, F. A., Rainville, F. M., Gardner, M., Parizeau, M., & Gagné, C. (2012). **DEAP: Evolutionary Algorithms Made Easy**. *Journal of Machine Learning Research*, 13, 2171-2175.
- Richardson, A. (2018). **Seismic Full-Waveform Inversion Using Deep Learning Tools and Techniques**. *arXiv preprint*. 

## 🔖 Disclaimer 🔖

All experiments were conducted on running **Debian GNU/Linux 12 (Bookworm)**.

📣 **Multiprocessing is implemented.**

---
For further details, refer to the paper associated with this repository.