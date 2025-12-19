# SeisTracer

SeisTracer is a Python and C package for computing travel times and ray paths of seismic waves in 1D spherical Earth models. It uses dense precomputed tables, the Newton-Raphson method, and analytical expressions for efficient calculations without Earth-flattening transformations. Supports teleseismic, core, diffracted, and higher-order phases with amplitude computations and visualizations.

* **Author**: Caio Ciardelli, Northwestern University
* **Supervisor**: Prof. Suzan van der Lee
* **License**: GNU General Public License (GPL) version 3 or later

If you use SeisTracer, please cite the following AGU poster:

Ciardelli, C. & Van der Lee, S. (2023). PyTracer: Ray Tracing and High-Frequency Synthetic Seismograms for 1D Earth Models. American Geophysical Union, Fall Meeting 2023, Board 0330, San Francisco - USA.

---

## Features

* Load and visualize 1D velocity models (Vp, Vs, density vs. depth).
* Precompute tables of ray parameters and distances for various models and depths.
* Compute travel times, ray paths, take-off angles, and amplitudes for specified or all phases.
* Support for minor/major arcs, multiple loops, and diffracted phases (e.g., Pdiff, Sdiff).
* Plot ray paths with colored branches and discontinuity circles.
* Generate travel-time curves for phases.
* Compute and plot amplitude coefficients at discontinuities (surface, Moho, etc.).

---

## Installation

SeisTracer requires no traditional installation beyond compiling the C components and setting up the Python path. Follow these steps:

**1. Obtain the Package**

Extract the package directory to your desired location (e.g., `~/seistracer`).

**2. Install Dependencies**

Ensure the following are installed:

* **Bash**: For running scripts (e.g., `create_tables.bash`).
* **GCC**: For compiling C binaries and libraries (`tables.c`, `tracer.c`, etc.).
* **Python 3.8+**: For the main module (`seistracer.py`).
* **Python Libraries**: Install via pip:
    ```bash
    pip install numpy scipy matplotlib
    ```
* **Make**: For compilation via Makefile.

**3. Compile the Package**

Navigate to the package directory and compile the C components using the provided `Makefile`:

```bash
cd seistracer
make all
```

This creates a `bin/` directory with `tables` and a `lib/` directory with shared libraries (`libcoefficients.so`, `libio.so`, `libtracer.so`). To clean previous builds:

```bash
make clean
```

**4. Set Python Path**

Add the package directory to PYTHONPATH to import the module:

```bash
export PYTHONPATH="$HOME/seistracer"
```
Add this to your `~/.bashrc` or `~/.zshrc` for persistence.

**5. Build Precomputed Tables**

After compilation, generate tables for the provided models:

```bash
./create_tables.bash ak135f iasp91 prem rem1d stw105
```

This creates tables in `tables/MODEL/` for depths 0 to 700 km (in 1 km increments) for each model. Tables are required for computations.

-----

## Directory Structure

```
seistracer/
├── create_tables.bash # Script to build precomputed tables
├── docs/ # Documentation and examples
│ └── example.py # Example Python script
├── __init__.py # Makes directory a Python package
├── LICENSE # GPL v3 license
├── Makefile # Compilation instructions
├── models/ # 1D Earth models
│ ├── ak135f.model
│ ├── iasp91.model
│ ├── prem.model
│ ├── rem1d.model
│ └── stw105.model
├── seistracer.py # Main Python module
├── README.md
├── setup/ # Configuration for C components
│ └── constants.h # Constants header
├── src/ # C source code
│ ├── headers/ # Header files (amplitudes.h, enums.h, etc.)
│ ├── shared/ # Shared utilities (amplitudes.c, io.c, etc.)
│ ├── tables.c # Table generation program
│ └── tracer.c # Ray tracing backend
└── tables/ # Precomputed tables (generated)
    ├── ak135f/
    ├── iasp91/
    ├── prem/
    ├── rem1d/
    └── stw105/
```

-----

## Usage

After setup, use the Python module for computations. Import `seistracer` and create a `Tracer` object with a model (e.g., 'iasp91'). The module handles loading precomputed tables, phase calculations, and visualizations.

### Running the Example Script

A complete demonstration of SeisTracer's capabilities is provided in `docs/example.py`. This script showcases model loading, phase computations, visualizations (model, ray paths, travel-time curves, amplitude coefficients), and data retrieval methods.

To run the example:

```bash
python docs/example.py
```
The script performs the following:
- Loads the IASP91 model and plots Vp, Vs, and density vs. depth.
- Computes and visualizes simple phases (P, S), diffracted phases (Pdiff, Sdiff), and complex phases (e.g., reflected, converted, core phases).
- Demonstrates minor/major arcs, multi-loop phases, amplitude computations, and custom phase queries.
- Prints phase information and retrieves path/data dictionaries.

Detailed explanations of the example script and its output are provided in Section 7 of the Manual (Manual.pdf).

-----

## Velocity Models

The `models/` directory includes standard 1D Earth models:

  * `ak135f.model`
  * `iasp91.model` (default)
  * `prem.model`
  * `rem1d.model`
  * `stw105.model`

Users can add custom models in the same format (depth, density, Vp, Vs, with header lines).

-----

## Development

The package is under active development with planned features:

  * Support for 3D velocity models.
  * Inversion for velocity structure.
  * Enhanced amplitude computations.
  * Parallel processing for large datasets.
  * GUI for visualizations.

Development is hosted on GitHub in the [caio.ciardelli/seistracer repository](https://github.com/caiociardelli/seistracer).

-----

## Contact

For questions, suggestions, or bug reports, contact:

**Caio Ciardelli**
*Email*: `caio.ciardelli@gmail.com`
Northwestern University, Department of Earth and Planetary Sciences

-----

## Acknowledgments

  * Prof. Suzan van der Lee for supervision.
  * Northwestern University for support.
  * Open-source communities for NumPy, SciPy, Matplotlib, and ctypes.
