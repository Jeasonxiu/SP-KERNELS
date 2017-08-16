% A master script to perform Sp migration or inversion.
% N. J. Mancinelli -- June 2017
%
function master(label,mode1or2,lab_amplitude,lab_wavelength,lab_depth,skipsta)
%clear classes

runname=sprintf('OUTPUT/%s/%s',date,label);

%% Link Dirs
curdir=pwd;
mkdir(runname)
cd(runname)

DirsToLink = {'MATLAB-CODE/FUNCTIONS', 'MATLAB-CODE/CLASSES', 'KERNEL-SEM/TEST_MODELS', 'MATLAB-CODE/TOOLS'};

for Dir = DirsToLink;
    cmd=sprintf('ln -sf %s/%s .', curdir, Dir{1});
    system(cmd);
end

addpath('FUNCTIONS');
addpath('CLASSES');
addpath('TOOLS')
%% Launch inversion

I=Inversion();
I=SetDefaultInversionParams(I, mode1or2);
I=SetUpKernel(I);
I=SetDataParams(I,lab_amplitude,lab_wavelength,lab_depth,skipsta);
I=SetUpMatrices(I);
I=RunInversion(I);
%ITest=I.testForwardModel(I);
Save(I);
%plot_model(ITest.VelocityModel2D,'test')
%plot_model(I.VelocityModel2D)
cd(curdir)
end
