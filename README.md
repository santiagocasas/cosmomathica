What is cosmomathica?
=====================

Wolfram Mathematica is a powerful and convenient software package used by
many cosmologists everywhere. However, since it is not always the most
efficient for low-level computing, most popular algorithms for numerical
computations in cosmology have been written in C or Fortran, since these
languages are typically much better suited for the task at hand. This
package bundles the functionality of some of these algorithms in
a Mathematica package, which makes them easier to use and avoids the need to
learn C or Fortran.

It makes use of the MathLink technology, which means that the original
source code remains untouched. The algorithms are not reimplemented either,
because that would negate the advantage of having them implemented in
a lower level language in terms of efficiency. On top of that, it would
inevitably introduce new bugs given the complexity of these algorithms.

For further details on how MathLink works, see
http://arxiv.org/abs/1107.4379.



Requirements
============

  * Wolfram Mathematica version >=10

  * The GNU Compiler Compilation (gcc) and (gfortran) version >=4.6
  
  * GSL (GNU Scientific Library) 
    You need the developer libraries for GSL (gsl-bin, libgsl-dev)
    
  * Universally Unique ID library 
    (uuid-dev) and (libuuid1)
    
  * Note: Library names could depend on your operating-system. Check compiler errors for the needed libraries.



External programs
=================

Due to copyright issues, the external packages cannot be delivered as part
of this package. You need to download the source code of the following
programs yourself. To make sure you got the right version, the MD5 sums of
the files in question are given here.

For CAMB and CLASS, the best option is to clone them from their public git repositories.

  * Transfer function: 
    tf_fit.c
    38ef737fe3bab405bac17db78815559f
    power.c
    32c6a5912acd36f5522d0a81c8f8c5d7
    http://background.uchicago.edu/~whu/transfer/transferpage.html
    Daniel Eisenstein and Wayne Hu 

  * Halofit+:
    halofit+.tar
    4ce61ec6504a2f7cf750ab007142dea8
    http://www.roe.ac.uk/~jap/haloes/
    Robert Smith et al, Martin Kilbinger

  * CosmicEmulator "Coyote" (v1.1):
    CosmicEmu_v1.1.tar.gz   
    a3c7da2b41152b7d30ba458e56f7e4ab  
    http://www.lanl.gov/projects/cosmology/CosmicEmu/emu.html
    Earl Lawrence

  * FrankenEmu
    CosmicEmu_v2.tar.gz
    1018867d84d9e6820a9ae8ee2f213b06
    http://www.hep.anl.gov/cosmology/CosmicEmu/CosmicEmu_v2.tar.gz
    Earl Lawrence
  
  * CAMB
    https://github.com/cmbant/CAMB.git
    Anthony Lewis et al

  * CLASS
    https://github.com/lesgourg/class_public.git
    Julien Lesgourges et al 

The (extracted) files need to placed in directories with names `tf`,
`halofit`, `camb`, `CosmicEmulator`, `FrankenEmu`, and `class` 
respectively, inside the `ext` directory.

If you use a modified
version of CAMB, cosmomathica might stil work, provided you have the same set of parameters. In general, if you modify things like the
type `CAMBparams`, you almost certainly will have to modify cosmomathica as
well.


How to compile the MathLink
===========================

The paths to the required Mathematica libraries are often different on every
system, so first you should adjust the Mathematica-related lines in the
`Makefile_mathlink` file, which is located in the `ext` directory. 
For example:
`MLINKDIR = /usr/local/Wolfram/Mathematica/11.2/SystemFiles/Links/MathLink/DeveloperKit`
Next, change into the `ext` directory and type `make`. 

CAMB uses the Intel Fortran compiler by default. Due to compatibility issues
with linking the different object files, CAMB needs to be compiled with GNU
gfortran. Please make sure that CAMB inside the folder camb is compiled with gfortran, by changing the respective Makefile. 


In `ext/halofit/Makefile` do this if you get an error message of libraries:
put the library flags beyond the .o files:

`$(CC) -o smith2demo smith2demo.o smith2.o $(lflags)`

Note that all warnings and errors that you see may be cause by the external
programs. Make sure the warning is not from compiling the MathLink and
contact the respective authors if you have concerns. Otherwise, feel free to
file a bug report on GitHub.


How to use Cosmomathica
=======================

In general, each one of the functions `CAMB`, `Transfer`, `Halofit`, `Class`
and `FrankenEmu` returns a list of replacement rules containing the raw data
as computed by the respective program. See the notebook `demo.nb` for
a demonstration. 


Version
=======

This is version 1.0. Since this is an alpha status, backwards compatibility
may be broken in future releases.

This is a forked version by Santiago Casas, with some major and minor modifications.


Copyright and licensing
=======================

Cosmomathica is released under the GPL2. Contributions are welcome. Note
that the copyright of the external software packages belong to their
respective owners. Read their copyright remarks before using them.
