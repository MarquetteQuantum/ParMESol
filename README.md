Here you can find the preliminary code and some subcodes used for the calculation of the rate coefficient using accurately calculated rotational vibrational states.

The files included and their roles are as follows:

1. barrier.py : This python code is used to interpolate the values of barrier energies from original barrier energy file found in OZONEData which are calculted only for every other J and every other K.
2. barrier_energies_module.f90 :  Uses aforementioned calculated barrier energies and gives them to the main code as barrier energy matrix
3. compile : Make file the produces the executable(named RK here). Should be adjusted to suit environment.
4. parmesol.f90 : Main code which carries out rate coefficient calculations. Here the calculations are done for Ozone at 1 bar and room temperature. Be mindful to change directies of OZONEData and barrier energies file.
5. histogram.py : File for plotting the distribution of partition function vs resonance width(used for analysis of data)
