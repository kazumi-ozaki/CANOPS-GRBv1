This Fortran code was written by Kazumi Ozaki for examining the evolution and stability of O2 levels
in the ocean-atmosphere system. The results are submitted to Geoscientific Model Development.
The basic design of the model is based on the CANOPS model that is originally developed by Kazumi
Ozaki (Ozaki and Tajika, 2013 EPSL; Ozaki et al., 2019 Geobiology). The reader should consult 
these papers for a comprehensive explanation and empirical/theoretical basis. The model's 
framework was expanded for this study by including (1) mass balance calculation of O2 in the 
atmosphere, (2) mass balance calculation of sedimentary reservoirs, and (3) simplified global redox 
(O2) budget. The improvments are fully explained in the main text. 

This model is still undergoing regular development and it is recommended that potential users contact
the corresponding author (Kazumi Ozaki; kazumi.ozaki@sci.toho-u.ac.jp) to obtain the latest version.
This code is provided freely but with the requirement that prospective users contact the corresponding
author with their research plans to avoid parallel projects emerging.

The model is written in modular form: Constants.f90 summarizes the constants used in the model, and
main.f90 is a main code. 

The model can be compiled by typing, as follows:
[for gfortran compiler]
	gfortran Constants.f90 main.f90
[for Intel Fortran compiler]
	ifort Constants.f90 main.f90

You can run the model by './a' (gfortran) or 'main' (ifort)

Input files are located in the 'input' directory, and the output files will 
be placed in the 'output' directory once the model is run. Each simulation will take a few weeks
to complete on a desktop computer, depending on the performance of the CPU.