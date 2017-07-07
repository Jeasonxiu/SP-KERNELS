% A master script to perform Sp migration or inversion.
% N. J. Mancinelli -- June 2017
%
clear classes

runname=sprintf('OUTPUT/%s/TESTING',date);

%% Link Dirs
curdir=pwd;
mkdir(runname)
cd(runname)

DirsToLink = {'MATLAB-CODE/FUNCTIONS', 'MATLAB-CODE/CLASSES', 'MIGRA', 'KERNEL-SEM/TEST_MODELS'};

for Dir = DirsToLink;
    cmd=sprintf('ln -sf %s/%s .', curdir, Dir{1});
    system(cmd);
end

addpath('FUNCTIONS');
addpath('CLASSES')
%% Set inversion params and launch inversion

%IP=InversionParams();
%IP.direction=3;
%save IP;
I=Inversion(0);
%save I
%back_project_synthetics(IP);
cd(curdir)