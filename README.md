# mpi-bucky
mpi-BUCKy README

Tyler K. Chafin - University of Arkansas
tkchafin@uark.edu

## -----------------------------------

THIS VERSION OF BUCKY IS A MODIFICATION OF THE ORIGINAL SOURCE CODE. PLEASE SEE REFERENCES BELOW FOR THE ORIGINAL PUBLICATION AND PROGRAM DESCRIPTION. 

TO CONTRIBUTE OR FOR QUESTIONS, PLEASE CONTACT ME AT tkchafin@uark.edu

More details on this multithreaded implementation to come soon... 

## -----------------------------------

MPI-BUCKy is a multi-threaded implementation of BUCKy using the OpenMPI message passing interface. BUCKy is an open-source program for species tree estimation using Bayesian concordance analysis of multiple genomic loci, free of assumptions regarding the source of discordance among gene trees. This

Please also refer to the original authors' user manual, website  (http://www.stat.wisc.edu/~ane/bucky/), and please cite the following publications if this program is useful to you, to recognize the original authors: 

C. Ané, B. Larget, D.A. Baum, S.D. Smith, A. Rokas (2007). Bayesian estimation of concordance among gene trees. Molecular Biology and Evolution 24(2), 412-426.

B. Larget, S.K. Kotha, C.N. Dewey, C. Ané (2010). BUCKy: Gene tree / species tree reconciliation with the Bayesian concordance analysis. Bioinformatics 

### PREREQUISITES
You must have the following installed (this may be an incomplete list):

    g++ compiler 
    OpenMPI >1.6 (mpic++)

### COMPILATION
If you have g++ and OpenMPI installed, compile the software with these commands, where $HOME is the directory containing bucky.

   cd $HOME/mpi-bucky/src;
   make

This will compile programs mbsum and bucky.
I suggest making dynamic links in ~/local/bin if this is in your path.

This has only been tested in a Linux environment (Ubuntu 14.04)

### HELP
Type these commands for very brief help messages.

    mbsum --help
    bucky --help

NOTE: Number of chain iterations * number of runs cannot exceed 2^32

### EXAMPLE
Suppose that you have a directory where each file is of the form *.t and is a MrBayes output file.
Use mbsum to summarize each file.  Remove the first 1000 trees of each for burnin.

    for X in *.t; do mbsum -n 1000 $X; done

This will create a file named <filename>.in for each file named <filename>.t .
Warning!  It will overwrite files with the name <filename>.in if they exist.

Next, to run bucky with multiple threads (where $procs is an integer representing the number of threads to use):

    mpirun -np $procs bucky *.in

This will create a bunch of output files of the form run1.* .
You can pick your own root file name.

### YEAST EXAMPLE
To try the yeast example described in the Ane et-al MBE paper (with a much smaller number of updates),
you can try the following.

    cd $BUCKY_HOME/mpi-bucky/data/yeast
    mpirun -np 8 ../../src/bucky y???/*.in

### EXAMINING OUTPUT
Mac users who want to use standard Mac applications to read the output can open any of the output files
using TextEdit.
