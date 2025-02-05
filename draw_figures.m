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
f4=resize(f4,20,25);
save(f4);

%% (5) Benchmarking figure
f5=figure5;
save(f5);

%% (6) Migration scheme
f6=Figure6;


%% (7) Inversion scheme
f7=Figure7;
