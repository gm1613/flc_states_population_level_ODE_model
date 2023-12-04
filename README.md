# flc_states_population_level_ODE_model
Matlab code for simulating the ODE model for FLC transcription regulation and histone modification dynamics, developed in 
Nielsen et al., PNAS, 2023, Antisense transcription and PRC2 repression function in parallel during vernalization

Three Matlab files are provided:
flc_model_simulations.m contains code to perform simulations of the ODE model developed in the paper. This file also contains code for processing the simulation output, to plot model predicted changes in histone modification levels over time. The Matlab function ode15s is used to solve the system of ODEs.

flc_state_model_fluct.m is a Matlab function containing code to compute the time derivatives of all the model variables as specified in the model description. This function is designed to be passed to the ODE solver (in this case, ode15s). This function is intended for simulating the wild type case (ColFRI), COOLAIR defective mutants, and H3K27me3 nucleation mutants. The same function also allows simulating the different temperature conditions: constant 5 degrees cold, fluctuating mild, and fluctuating strong conditions (see supplementary information to Nielsen et al., PNAS, 2023, for a description of the model and these conditions)

flc_state_model_fluct_spreading_mutant.m is a Matlab function, similar to flc_state_model_fluct.m, but intended for simulating H3K27me3 spreading mutants
