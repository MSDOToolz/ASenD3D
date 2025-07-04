# ASenD3D

ASenD3D (Adjoint-based Sensitivity and Design in 3D) is a tool for performing finite element analysis of physical systems and computing sensitivities of objectives to design parameters/variables using the adjoint method.  High-efficiency, high-fidelity gradient-based optimization and uncertainty quantification in multidisciplinary design are the primary intended applications, but the package functions as a general open-source finite element analysis tool.  Capabilities include analysis of elastic, diffusive and thermal response of 3D solid, shell and beam structures under loading, as well as modal analysis for natural frequency and buckling problems.

Originally created at the University of Wyoming under government-funded research projects and formerly known as AStrO (Adjoint-based Structural Optimizer), ASenD3D continues to expand, with recent efforts in improving efficiency and usability of interface.  A python-based interface provides easy-to-use tools for model/input generation, job submission to the core solver, post-processing and visualization.  It also facilitates integration with other packages such as optimizers.

Basic example scripts are provided in the examples directory.  Users are encouraged to try these out and familiarize themselves, and seek to make contributions as they feel inclined.

# Setup

## Prerequisities

Depending on their scope of use and current setup, there are a few prerequisites that users may or may not need before setting up ASenD3D.

### 1. Install Anaconda/Python Environment

It is recommended to run ASenD3D within Anaconda, which provides all needed tools conveniently packaged and accessible, as well as support for creating specialized environments for different projects.  If not already installed, it can be downloaded for free from https://anaconda.org.

### 2. Install rustup to locally build the core analysis engine

The core analysis engine/solver of ASenD3D is written in the Rust programming language for speed, resource efficiency and robustness. To build locally, the offical package for build and compilation of rust projects, or rustup (www.rust-lang.org/tools/install), must be installed. On MacOS or Linux, it can be installed from the command line with

        curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh

Windows users can visit the above url for a download, or use the included build, in which case installing rustup is not required. Building locally is generally recommended, but not necessarily required since most users would only need to modify the core solver's source code for user-defined, progressive analysis.

### 3. Install Git

Especially if running on Windows, users may want to install Git, an application for conveniently organizing and tracking changes to complex code projects.  It is not strictly necessary, but it will be helpful if one wants to make changes to the interface modules or the core solver itself.  Downloads can be found at https://git-scm.com.

## Main Install Steps

### 1. Clone the Repository

To clone the ASenD3D repository onto your local machine, navigate to a directory of your choice in the git bash, or command line terminal in Linux and run

        git clone https://github.com/MSDOToolz/ASenD3D.git
		cd ASenD3D
        
Currently there is only one branch called main.  If preferred, users can just download and extract the package from the above url without utilizing git commands.  They just won't have the functionality for updating and tracking changes that the git interface offers.  

### 2. Installing the Python Interface

If desired, create an anaconda environment for ASenD3D, by running

        conda create -n asend-env (asend-env can be changed to any environment name of your choice)
        conda activate asend-env

from the anaconda powershell prompt, or just the command line terminal in Linux.  To install the interface for model and input file generation, job submission and post-processing, run

        pip install .
		
### 3. Initializing the Core Solver

To initialize the core solver to be invoked through the interface tools, run

        python solverSetup.py
		
from the root directory in the anaconda powershell prompt or the command line in Linux.
        
# Running ASenD3D

Upon completing the steps for setup, ASenD3D can be run by creating Python scripts such as those given in the examples directory, and running them in your environment from the command prompt or an IDE like Spyder.

# Documentation

Full documentation is planned to be made available at https://readthedocs.org, but it is presently under construction.  Please stand by for further updates.
