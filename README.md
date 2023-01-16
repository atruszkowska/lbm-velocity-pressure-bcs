# lbm-velocity-pressure-bcs
Lattice Boltzmann method with velocity and pressure boundary conditions for two-phase flows

This is a serial MATLAB code that can be used to simulate two-phase flows with LBM under velocity/pressure conditions. Current code simulates the flow of two fluids through a t-mixer. The mixer geometries, flow, and fluid properties that this code is capable of generating are found [here](https://ir.library.oregonstate.edu/concern/graduate_thesis_or_dissertations/5q47rt56m):

> Truszkowska, Agnieszka. "Multiscale modeling of two-phase flows in microarchitectures with microfeatures." (2014).

The code was run on Linux and tested with MATLAB 2022b.

To run the code, run the ```two_ph_velocity.m``` - this is also where all main modifications should go. The code ```tp_cont.m``` can be ran if the main code stops prematurely and needs to be continued from the last saved point. Other code are functions supporting the simulation that should not require any case-base changes.

This code was validated, however, the code structure is not sophisticated. The user is recommended to change it with care and encouraged to report any issues or bugs. The code here is one set out of a series of cases, if the user needs something else, reach out and the repository will be updated if possible.


