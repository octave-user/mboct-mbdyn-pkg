# mboct-mbdyn-pkg<sup>&copy;</sup>
**mboct-mbdyn-pkg** belongs to a suite of packages which can be used for pre- and postprocessing of MBDyn models (https://www.mbdyn.org) with GNU-Octave (http://www.gnu.org/software/octave/). This package provides interfaces between the multibody dynamics software MBDyn and the GNU Octave programming language.

# List of features
  - Generate MBDyn input files for arbitrary curved beam structures represented by Non-Uniform Rational B-Splines (NURBS).
  - Generate MBDyn input files for modal elements based on Finite Element models.
  - Generate MBDyn input files for solid elements based on Finite Element meshes.
  - Run the multibody dynamics solver MBDyn.  
  - Load output files from MBDyn.
  - Compute inertia properties of groups of several rigid and flexible bodies.
  - Compute frequency response functions of linearized equations of motion.
  - Scale elastic deformations of flexible bodies for post-processing.
  - Perform post-processing of elastohydrodynamic bearing data from MBDyn.
  - Perform post-processing of MBDyn models with modal elements.
  - Perform post-processing of MBDyn models with solid elements.

Copyright<sup>&copy;</sup> 2019-2025

[Reinhard](mailto:octave-user@a1.net)

# Installation
  - See [simple.yml](https://github.com/octave-user/mboct-mbdyn-pkg/blob/master/.github/workflows/simple.yml) as an example on how to install mboct-mbdyn-pkg.

### Usage
  - Run Octave.  
    `octave`
  - At the Octave prompt load the package.   
    `pkg load mboct-mbdyn-pkg`
  - At the Octave prompt execute a demo.  
    `demo mbdyn_post_ehd_load_output`
	
