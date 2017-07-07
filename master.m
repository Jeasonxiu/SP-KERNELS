% A master script to perform Sp migration or inversion.
% N. J. Mancinelli -- June 2017
%
function master(label,mode1or2)
%clear classes

runname=sprintf('OUTPUT/%s/%s',date,label);

%% Link Dirs
curdir=pwd;
mkdir(runname)
cd(runname)

DirsToLink = {'MATLAB-CODE/FUNCTIONS', 'MATLAB-CODE/CLASSES', 'KERNEL-SEM/TEST_MODELS'};

for Dir = DirsToLink;
    cmd=sprintf('ln -sf %s/%s .', curdir, Dir{1});
    system(cmd);
end

addpath('FUNCTIONS');
addpath('CLASSES')
%% Launch inversion

I=Inversion();
I=SetDefaultInversionParams(I, mode1or2);
I=SetUpKernel(I);
I=SetUpMatrices(I);
I=RunInversion(I);
Save(I);
%plot_model(I.VelocityModel2D)
cd(curdir)
end
