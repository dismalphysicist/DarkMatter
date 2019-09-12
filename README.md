# DarkMatter
Summer project 2019

A C++ program comprising two classes and a main program. Calculates the velocity-weighted cross section for the process DM DM -> f+ f-, where the final state is approximated as massless. 

# To Use:
Modify the file multichannel_integration.cpp to input the required parameters for dark matter mass, temperature and the dark matter and final state coupling constants. Modify the integration parameters to produce suitable functions in the Expo_fit class. Compile using the command g++ multichannel_integration.cpp integrand.cpp expo_fit.cpp -o [programname]. Run the executable to output the value of the velocity-weighted cross section and its error to the standard outstream (terminal). 
