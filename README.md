# ASenD3D

ASenD3D (Adjoint-based Sensitivity and Design in 3D) is a tool for performing finite element analysis of structural systems and computing sensitivities of structural objectives to design parameters/variables using the adjoint method.  High-efficiency, high-fidelity gradient-based optimization and uncertainty quantification in structural design are the primary intended applications, but the package functions as a general open-source finite element analysis tool.  Capabilities include analysis of elastic and thermal response of 3D solid, shell and beam structures under loading, as well as modal analysis for natural frequency and buckling problems.

Originally created at the University of Wyoming under government-funded research projects and formerly known as AStrO (Adjoint-based Structural Optimizer), ASenD3D continues to expand, with recent efforts in improving efficiency and usability of interface.  Several python-based utilities are under development to promote ease and covenience of use and integration with other packages such as optimizers and meshing tools.  These remain a work-in-progress and will come online in the months ahead.

Basic examples of usage are provided for starters as the public repository undergoes construction.  Users are encouraged to try these out and familiarize themselves, and seek to make contributions as they feel inclined.

# Setup

Presently, ASenD3D is set up to be run in a Linux environment, though support for Windows and MacOS is intended to be added in the near future.  That said, it's use of external libraries/packages is relatively minimal so setup is rather straight-forward, and should need only minor modifications for other operating systems.

The core solver and sensitivity tool of ASenD3D is written in Fortran, and users should be able to run from the pre-compiled binaries included in the 'bin' directory.  But if the user has a preferred compiler, or wishes to modify the main source code, they will need to compile it themselves from source.  Either way, the user must first clone the repository to their local machine as a first step.

## 1. Clone the repository

To clone the ASenD3D repository onto your local machine, navigate to a directory of your choice in the command line terminal and run

        $ git clone https://github.com/MSDOToolz/ASenD3D.git
        
Currently there is only one branch called main.  To build the core solver from source, proceed to step 2.  To run using the included binaries, skip to step 3.

## 2. Building From Source

If a user wishes to build ASenD3D themselves from source code, a Fortran compiler is needed.  Seasoned coders may have their own preference of compiler, and any choice should do.  If you don't have one installed already, a simple and convenient choice on Linux is the GNU Fortran compiler, or gfortan, which can be installed with

        $ sudo apt update
        $ sudo apt install gfortran

Once the compiler is installed, run the shell script buildFromSource.sh in the main directory.  You will likely need to give it permission first, with

        $ chmod +x buildFromSource.sh    (this command only needs called on first setup, and should not required for subsequent compiles)
        
And then execute, with

        $ ./buildFromSource.sh
        
The included script is set up for gfortran as the compiler, if another is being used the script will have to be be modified accordingly.

## 3.  Installing Python-Based Packages

A library of utilities are provided to assist with pre and post-processing, and other useful functions periferral to the main solver.  These are Python-based, so the user must have Python installed in some form as a pre-requisite.  A popular avenue is the Anaconda environment (https://www.anaconda.com), but standard installation as guided from https://www.python.org is more than adequate.  The needed packages can be installed with pip:

        $ pip install plotly==5.13.0
        $ pip install ruamel.yaml
        
or conda:

        $ conda install -c plotly plotly
        $ conda install -c conda-forge ruamel.yaml
        
# Running ASenD3D

The basic usage of ASenD3D is to create the necessary input files along with a job script detailing the analysis steps desired, and call the main solver executable in the form:

        $ ./bin/ASenD3D.exe path/to/job/script/jobScriptName.txt
        
There are three basic types of input files: 1) model input, 2) design variable input, and 3) objective function input.  Examples of all of these are provided for various basic test cases in the examples directory, along with several sample job scripts.  The main solver can be invoked directly from the command line, though many tools are provided in the pythonUtilities directory to help streamline the workflow process and make the most of the results.  Users are encouraged to build such tools of their own as their needs call for, and consider sharing them for others who may benefit.

# Documentation

Full documentation is planned to be made available at https://readthedocs.org, but it is presently under construction.  Please stand by for further updates.
