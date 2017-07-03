%%Draws figures for GJI paper
%N.J. Mancinelli -- June 2017
%% Sp Isochrons

%% Model Setup (use illustrator)

%% Numerical Kernels

%% Components of analytical kernels (Timing, Geometrical Spreading, Scattering Patterns)

%% Benchmarking figure

%% Layered Kernels (Timing, Geometrical Spreading, Scattering Patterns)
addpath('MATLAB-CODE/FUNCTIONS')
plot_kernels_layered()

%% Migration scheme
curdir=pwd;
cd('OUTPUT/MIGRATION_FOR_KAREN')
plot_model('model_1')
cd(curdir)

%% Inversion scheme
curdir=pwd;
cd('OUTPUT/INVERSION_FOR_KAREN')
plot_model('model_1')
cd(curdir)
