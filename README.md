# AeroelasticModels
This repository contains all the code used in the doctoral thesis developed at ITA.

The Folders:
-> 00_odeSolver: Contains the runge-kutta solver used.
-> 01_FEM_Software: Contains Finite-Element Models used to generate the wings data for coupled models. The data used for the thesis is included pre-generated.
-> 02_Wagner_linear_aero: Linear aerodynamic model
-> 03_LDVM_nonlinear_aero: Non-linear aerodynamic model LDVM
-> 04_UVLM_nonlinear_aero: Non-linear aerodynamic model UVLM in various forms.

The Functions:
Dynamics_LDVM_flatplate
dynamics_lin
Dynamics_UVLM_flatplate
Dynamics_UVLM_flatplate_LEV
Implement the respective couple models that are solved by the ode4 solver.

The Scripts:
-> UVLMflatplate_Prescribed (Runs UVLM prescribed model)
-> LDVMflatplate_testbench (Runs LDVM and Wagner coupled models)
-> UVLMflatplate_testbench (Runs UVLM coupled models)
Are the start point for running the simulations.

The Scripts:
-> LDVMFlatplate_testbench_state
-> UVLMFlatplate_testbench_state
Allow to resume the respective previous simulation from a saved state.
