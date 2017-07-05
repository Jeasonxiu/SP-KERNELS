%%Draws figures for GJI paper
%N.J. Mancinelli -- June 2017
addpath('MATLAB-CODE/FUNCTIONS')
addpath('MATLAB-CODE/CLASSES')
addpath('MIGRA')

%% Sp Isochrons
f1=figure1();
save(f1);

%% Model Setup (use illustrator)

%% Numerical Kernels (show timing, angular dependence)
f3=figure3;
save(f3);

%% Components of analytical kernels (Timing, Geometrical Spreading, Scattering Patterns)
f4=figure4;
save(f4);

%% Benchmarking figure
f5=figure5;
save(f5);

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
