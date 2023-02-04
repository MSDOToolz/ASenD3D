# ASenD3D

ASenD3D (Adjoint-based Sensitivity and Design in 3D) is a tool for performing finite element analysis of structural systems and computing sensitivities of structural objectives to design parameters/variables using the adjoint method.  High-efficiency, high-fidelity gradient-based optimization and uncertainty quantification in structural design are the primary intended applications, but the package functions as a general open-source finite element analysis tool.  Capabilities include analysis of elastic and thermal response of 3D solid, shell and beam structures under loading, as well as modal analysis for natural frequency and buckling problems.

Originally created at the University of Wyoming under government-funded research projects and formerly known as AStrO (Adjoint-based Structural Optimizer), ASenD3D continues to expand, with recent efforts in improving efficiency and usability of interface.  Several python-based utilities are under development to promote ease and covenience of use and integration with other packages such as optimizers and meshing tools.  These remain a work-in-progress and will come online in the months ahead.

Basic examples of usage are provided for starters as the public repository undergoes construction.  Users are encouraged to try these out and familiarize themselves, and seek to make contributions as they feel inclined.

# Setup

Presently, ASenD3D is set up to be run in a Linux environment, though support for Windows and MacOS is intended to be added in the near future.  That said, it's use of external libraries/packages is relatively minimal so setup is rather straight-forward, and should need only minor modifications for other operating systems.

The core solver and sensitivity tool of ASenD3D is written in Fortran, and users should be able to run from the pre-compiled binary included in the 'bin' directory.  But if the user has a preferred compiler, or wishes to modify the main source code, they will need to compile themselves from source.  Either way, the user must first clone the repository to their local machine, my navigating from the command line to a directory of their choice and running

        $ git clone https://github.com/MSDOToolz/ASenD3D

## Building From Source

        $ sudo apt update
        $ sudo apt install gfortran

        $ pip install plotly==5.13.0
        $ pip install ruamel.yaml

## Running With Included Binary

If the user is content to simply run ASenD3D from the included binary and there are no issues with doing so, they need only make sure they have the needed python packages for the python-based utilites that are meant to assist with pre and post-processing, and interfacing with the core solver by running

        $ pip install plotly==5.13.0
        $ pip install ruamel.yaml

# Documentation
