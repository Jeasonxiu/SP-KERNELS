%%Draws figures for GJI paper
%N.J. Mancinelli -- June 2017
%% Sp Isochrons


%% Benchmark Test
%

%% Components (Timing, Geometrical Spreading, Scattering Patterns)

%
%% Layered Kernels
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