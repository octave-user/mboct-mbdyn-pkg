# mboct-mbdyn-pkg<sup>&copy;</sup>
**mboct-mbdyn-pkg** belongs to a suite of packages which can be used for pre- and postprocessing of MBDyn models (https://www.mbdyn.org) with GNU-Octave (http://www.gnu.org/software/octave/). This package provides interfaces between the multibody dynamics software MBDyn and the GNU Octave programming language.

# List of features
  - Generate MBDyn input files for arbitrary curved beam structures represented by Non-Uniform Rational B-Splines (NURBS).
  - Generate MBDyn input files for modal elements based on Finite Element models.
  - Run the multibody dynamics solver MBDyn.  
  - Load output files from MBDyn.
  - Compute inertia properties of groups of several rigid and flexible bodies.
  - Compute frequency response functions of linearized equations of motion.
  - Scale elastic deformations of flexible bodies for post-processing.
  - Do post-processing of elastohydrodynamic bearing data from MBDyn.

Copyright<sup>&copy;</sup> 2019-2020

[Reinhard](mailto:octave-user@a1.net)

# Installation

## GNU Octave installation
  - Follow the instructions on (http://www.gnu.org/software/octave/) to install GNU Octave.  
  - Make sure, that `mkoctfile` is installed.  
    `mkoctfile --version` 

### MBDyn installation:
  - Clone the source tree of MBDyn.  
    `git clone https://public.gitlab.polimi.it/DAER/mbdyn.git -b develop`
  - Compile and install MBDyn.  
    `cd mbdyn`  
    `./bootstrap.sh`  
    `./configure CXXFLAGS=-O3 --enable-octave --enable-autodiff --with-static-modules --with-umfpack`  
    `make`  
    `make install`

### GNU Octave package installation:
  - Make sure that the GNU Octave nurbs package is installed.  
    `octave --eval 'pkg install -forge nurbs'`
  - Install the following packages from github.  
    `for pkg in octave mbdyn; do`    
        `git clone https://github.com/octave-user/mboct-${pkg}-pkg.git && make -C mboct-${pkg}-pkg install_local`	  
    `done`

### Usage
  - Run Octave.  
    `octave`
  - At the Octave prompt load the package.   
    `pkg load mboct-mbdyn-pkg`
  - At the Octave prompt execute a demo.  
    `demo mbdyn_post_ehd_load_output`
	
