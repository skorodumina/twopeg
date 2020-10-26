Generator is exclusive and electroproduction (but at low Q2 it uses photoproduction data as a limiting point). It works for all Q2 starting from very small (0.005 GeV^2) and W from the threshold to 4.5 GeV.

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

read_xsect_files.cxx - reads the files with cross sections from the corresponding directories, 

out_file_open.cxx, out_file_fill.cxx, out_file_close.cxx - are needed to open fill, and close output files, r espectively,

anti_rot.cxx - contains the subroutine that calculates the 4-vectors of all final particles in the lab frame from  W, Q2 and hadron variables generated in the CMS,

rot.cxx - contains the subroutine that calculates hadron variables in the CMS from their 4-vectors in the lab frame.

radcorr.cxx - performs radiative effects.

--------------------------------------------------

There are two compiling options: 

1) "make nobos" - to compile without BOS libraries, no output in BOS format is possible in this case (suitable for CLAS12),
2) "make bos" - to compile with BOS libraries (suitable for CLAS). BOS output can be created according to the flag in the input file. Some of the libraries are located at "/u/home/gleb/lib/LinuxRHFC8" and sooner or later they will become irrelevant.

It should compile on ifarms.

------------------------------------------------

The generator supports two options for taking the input parameters:

1) As a cin stream from a certain input file. To use this option run as "./twopeg_bos.exe < inp1" (or "./twopeg < inp1"), where inp1 is in the same directory and contains input parameters with comments. 

2) The command line input is also supported. For this option all input parameters are automatically taken from "data/inp_cmd_line". The EG should then be run as "./twopeg_bos.exe a" (or "./twopeg a"), where "a" is any charachter. The number of events to generate can also be reset with "./twopeg_bos.exe --trig Nevents" (or "./twopeg --trig Nevents"), where "Nevents" is the desired number of events.

--------------------------------------------------

The generator needs "dat" files with tabulated structure functions and fit parameters. They are located in the "data" subfolder inside the EG directory. If you move it, you need to define the environment variable "data_dir_2pi" that points to the new folder location (e.g. in csh use "setenv data_dir_2pi new_path/").

--------------------------------------------------

See more details in CLAS12-NOTE-2017-001 (arXiv:1703.08081) and CLAS12-NOTE-2017-014 (arXiv:1712.07712).
See also Iu. Skorodumina's wiki page: https://clasweb.jlab.org/wiki/index.php/TWOPEG_event_generator  

--------------------------------------------------

Contact persons: Iuliia Skorodumina (skorodum@jlab.org) and Gleb Fedotov (gleb@jlab.org)
 
