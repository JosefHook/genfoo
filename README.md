# GenFoo : a generalized Fokker-Planck solver 

Contents

    1. Installation
    2. Building an operator
    3. License

-------------------------------------------------------------------------------
Installation

1.

From a fresh ubuntu 12.04

Install the following packages:
sudo apt-get install bzr cmake cmake-curses-gui fenics libxslt1.1 libxslt1-dev
 

Footprint: around 600MB


If you intend to use the Monte Carlo factories you will need to execute

 ./bootstrap.sh 
 
This script will fetch the additional packages, 
ANN, h5Part and boost-numeric-bindings.


If you encounter problems with the h5Part package,prompting 
something about MPI then you probably have a serial version 
of HDF5 installed. To fix this type: sudo apt-get install libhdf5-openmpi-dev

1.1

Select solvers with:  ccmake .
By default only the FEM factory will be built.


1.2
Compile the source:
1. cmake .
2. make install


IF everything went well you should have a working copy of GenFoo in

${GENFOO_INSTALL_DIR}/Build


Test GenFoo by 

cd int Build/ and 

./GenFoo ../Params/Orn2D.xml for a FEM factory example
./GenFoo ../Params/Duffing2D.xml for a MonteCarlo factory example

Note that the DeltaFMonteCarlo and Adaptive FEM factories are
currently under development and will most likely not work 
as expected. However the FEM and MonteCarlo factories
are "stable"







------------------------------------------------------------------------------
2. Building an operator

GenFoo is shipped with a couple of operator demos in the Operator 
directory that can be glanced upon when building a new operator. 
However many of them are currently under active development.

When building a new operator one can either write in c++ or directly 
in XML.

Below is an example XML operator of the Ornstein-Uhlenbeck process:

<?xml version="1.0"?>
<Operator>
    <Private>
       	<variable type="double" static="true" const="true">
          <documentation> Our private parameter </documentation>
          <name> kappa </name>
          <value> 9.0 </value>
          <unit> 1 </unit>
       </variable>
    </Private>
    <!-- See RabbitEarsModel.xml for implementation of private functions -->
    <!-- Number of dimension, here a 2D problem --> 
    <Dimensions>
            2
    </Dimensions>

    <!-- The FP drift vector -->
    <!-- The variable x is here two-dimensional x[0],x[1] -->
    <Drift>
        <component index="0">
                <value>   -2.0*(x[0]-1.5) </value>
        </component>
        <component index="1">
                <value>   -2.0*(x[1]-1.5) </value>
        </component>
    </Drift>
    <!-- The FP diffusion matrix -->
    <!-- with elements | 2kappa  0     |
                       | 0      2kappa |
    -->
    <Diffusion>
      <code>
        int foo = 1.0;
        </code>
        <component indexColumn="0" indexRow="0">
                <!-- the option diffusion  -->
                <value> 2.0*kappa </value>
        </component>
        <component indexColumn="1" indexRow="1"> <!-- the volatility diffusion -->
                <value> 2.0*kappa </value>
        </component>
    </Diffusion>

    <!-- If a source term is present add it here -->
    <Source>
      <value>0.0</value>
    </Source>

    <!-- The initial condition here a Gaussian centered in (0,0)
    <InitialCondition>
      <value>
            20.0*1/(2.0*M_PI)*sqrt(1/(2.0*M_PI))*exp(-( std::pow(x[0]-0.0, 2)
      	     + std::pow(x[1]-0.0,2) )/0.2 )
      </value>
    </InitialCondition>

</Operator>




------------------------------------------------------------------------------
3. License

GenFoo is licensed under the GNU LGPL Version 2.1, see http://www.gnu.org 
for non-commmercial and academic use. 



-------------------------------------------------------------------------------
