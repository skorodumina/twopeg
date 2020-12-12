TWOPEG is an event generator for the pi+pi- electroproduction off protons.

The generator is exclusive and electroproduction (but at low Q2 it uses photoproduction data as a limiting point). It works for all Q2 starting from very small (0.0005 GeV^2) and W from the threshold to 4.5 GeV.

The generator performs weighted event generation: first, it generates phase space (flat) distributions and then for each event multidimensional cross section is applied as a weight (see CLAS12-NOTE-2017-001). 

The generator includes radiative effects according to Mo&Tsai. It also includes Fermi smearing if needed (see CLAS12-NOTE-2017-014).

--------------------------------------------------

Output files:

1) The generator can produce output in lund format with minor changes:  in the Header in field 6* instead of "x" the event number is put, and in field 10* instead of "nu" the cross section value is put (see e.g. this slide https://clasweb.jlab.org/wiki/images/d/dc/Hybrid_Jan15_2016.pdf ),

2) The generator also can make BOS output that is compatible with CLAS6 reconstruction procedure,

3) Two ROOT outputs: "out_hist_test.root" with some histograms for testing purpose and "tree_sigma.root" with the cross section weights for each event (stored as a tree).

--------------------------------------------------

The code incorporates rather large number of files, but the major part of them corresponds to the interpolation/extrapolation procedures for the cross section calculation (all files with names get_xsect* and interpol*).

The key files are:

main_prog.cxx - the main program,

input.cxx - reads input parameters,

read_xsect_files.cxx - reads the files with cross sections from the "data" sub-folder, 

out_file_open.cxx, out_file_fill.cxx, out_file_close.cxx - are needed to open fill, and close output files, respectively,

anti_rot.cxx - contains the subroutine that calculates the 4-vectors of all final particles in the lab frame from  W, Q2 and hadron variables generated in the CMS,

rot.cxx - contains the subroutine that calculates hadron variables in the CMS from their 4-vectors in the lab frame.

radcorr.cxx - performs radiative effects.

--------------------------------------------------

There are two compiling options: 

1) "make nobos" - to compile without BOS libraries, no output in BOS format is possible in this case (suitable for CLAS12),
2) "make bos" - to compile with BOS libraries (suitable for CLAS). BOS output can be created according to the flag in the input file. Some of the libraries (-lc_bos_io -lmapmanager -lfputil -lfpack -lbos -lbankdefs) are taken from "/u/home/skorodum/lib/bos/bos_gcc920". Note that these libraries were compiled for gcc 9.2.0, and so may not match the older compliper versions. For gcc 5.3.0 the libraries can be taken from "/u/home/skorodum/lib/bos/bos_gcc530",the path in the Makefifle should be then changed accordingly.

It should compile on ifarms.

------------------------------------------------

The generator supports two options for taking the input parameters:

1) As a cin stream from a certain input file. To use this option run as "./twopeg < inp1", where "inp1" is in the same directory and contains input parameters with comments. 

2) The command line input is also supported. The EG can then be run as "./twopeg a", where "a" is any charachter. The default values of input parameters are then used unless reset through the cmd line. For cmd line option and default values type "./twopeg --help".

--------------------------------------------------

The generator needs "dat" files with tabulated structure functions and fit parameters. They are located in the "data" subfolder inside the EG directory. If you move it, you need to define the environment variable "TWOPEG_DATA_DIR" that points to the new folder location (e.g. in csh use "setenv TWOPEG_DATA_DIR new_path/").

--------------------------------------------------

See more details in CLAS12-NOTE-2017-001 (arXiv:1703.08081) and CLAS12-NOTE-2017-014 (arXiv:1712.07712).
See also Iu. Skorodumina's wiki page: https://clasweb.jlab.org/wiki/index.php/TWOPEG_event_generator  

--------------------------------------------------

Contact persons: Iuliia Skorodumina (skorodum@jlab.org) and Gleb Fedotov (gleb@jlab.org)
 
