# Chromosome_modelling
C++ software tool for 3D chromosome modelling using HiC data

--------------------------------
## Table of Contents

### •	 About ‘Chromosome-condensation’

### •	 Getting Started

### • Tutorial

### •	Output & visualization

### •	Documentation

### •	License

### •	Contact

### •	Acknowledgements

------------------------------------

### About ‘Chromosome-reconstruction’

‘Chromosome-reconstruction’ code is a C++ biophysical software designed for 3D chromosome modelling from HiC data. 


### Getting Started

This is an example of how you may set up the code running for your projects locally.

#### Prerequisites
-	C++ built-in
-	Visualization software such as Chimera, VMD or PyMOL

#### Installation
1.	You can download the code (`initConfig.hpp, initConfig.cpp, initDynamics.hpp, initDynamics.cpp, chromoCell.cpp, Makefile, chromosome_real.numberBead, chromosome_real.normMatrix`)
2.	Or you can clone the repo: `git clone` https://github.com/FrancisCrickInstitute/Chromosome-Condensation


### Tutorial

First, create a folder where you want to run simulations:

`mkdir test_reconstruct`


`cd test_reconstruct`

Make sure you copy all the code files from chromo_3D_code to the folder where you want to run simulations:

`cp path_to_code_folder/chromo_3D_code/* ./ `

Next, set up parameters in the initConfig.hpp (for parameter description see **Documentation**) and then compile the code:

`make`

Prepare normalized matrix of a single chromosome or whole genome to be modelled and name it *_chromosome_real.normMatrix_*. You can find a template of single chromosome normalized Matrix together with other files.
Also, do include number of reconstructed chromosomes and number of bins for each in *_chromosome_real.numberBead_*. 

Run the simulation:

`./chromo _chromosome_real`


### Output & visualization

The main output of code are PDB files from specific timepoints with coordinates of reconstructed chromosomes. These can be visalized using Chimera, VMD or PyMOL.  


### Documentation – parameters 

The ‘Chromosome-condensation’ code consists of several files:

- **initConfig.hpp** = includes parameter settings and declaration of functions that initialise the chromatin chain
- **initConfig.cpp** = includes definition of functions that initialise the chromosome(s)
- **initDynamics.hpp** = includes declaration of functions that describe the dynamical evolution of chromosome chain(s)
- **initDynamics.cpp** = includes definition of functions that describe the dynamical evolution of chromosome chain(s) to converge to energy minimum
- **chromoCell.cpp** = includes functions to call functions in files above to simulate the dynamical evolution and convergence 

General parameter set up can be adjusted in **initConfig.hpp** file. **README_parameters_chromo_modelling.xlsx** contains more in-depth description of polymer simulation parameters.



Optional output files and their corresponding parameters are:
•	**_chromoPDB_*.pdb** = PDB file format of chromatin chain with binders/condensin in selected time point


### License
Distributed under **The Francis Crick Institute License**. 
 

### Contact
-	**Xiao Fu** - @foolbirdie, xiao.fu@crick.ac.uk
-	**Tereza Clarence** - @ClarenceTereza, tereza.clarence@crick.ac.uk, clarence.tereza@gmail.com
-	**Paul Bates** - @PaulBatesBMM, paul.bates@crick.ac.uk


### Acknowledgements

This code was developed under the **Biomolecular Modelling Laboratory** https://www.crick.ac.uk/research/labs/paul-bates at Francis Crick Institute (https://www.crick.ac.uk/). 
Please cite our paper (doi:10.21203/rs.3.rs-757454/v1) when using the code.

