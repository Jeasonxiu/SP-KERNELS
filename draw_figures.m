%%Draws figures for GJI paper
%N.J. Mancinelli -- June 2017
addpath('MATLAB-CODE/FUNCTIONS')
addpath('MATLAB-CODE/CLASSES')

%% (1) Sp Isochrons
f1=figure1();
save(f1);

%% (2) Model Setup (use illustrator)

%% (3) Numerical Kernels (show timing, angular dependence)
f3=figure3;
save(f3);

%% (4) Components of analytical kernels (Timing, Geometrical Spreading, Scattering Patterns)
f4=figure4;
save(f4);

%% (5) Benchmarking figure
f5=figure5;
save(f5);

%% (6) Migration scheme
curdir=pwd;
cd('OUTPUT/07-Jul-2017/1')
load I
cd(curdir)
plot_model(Inversion.VelocityModel2D,'figure6')


%% (7) Inversion scheme
curdir=pwd;
cd('OUTPUT/07-Jul-2017/2')
load I
cd(curdir)
plot_model(Inversion.VelocityModel2D,'figure7')

