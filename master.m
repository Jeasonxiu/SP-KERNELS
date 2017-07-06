% A master script to perform Sp migration or inversion.
% N. J. Mancinelli -- June 2017
%
runname=sprintf('OUTPUT/%s/SPACING_TEST_DBL',date);

InversionParams.skipSta=10;
InversionParams.tDeci=20;
InversionParams.dxin=5;
InversionParams.dzin=2;
InversionParams.Nus=[3.15E+0,3.15E-1,3.15E-2];
InversionParams.Norm_Opts=[2];
InversionParams.nIterMax=500;
InversionParams.direction=3; %3 (synthetics from both directions)
InversionParams.Kernel_Type=3;  %3 (analytical) should be best
InversionParams.ImagingMethod=2; %1 for Back Proj, %2 for CG Inversion
InversionParams.saveFilename='model';

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

back_project_synthetics(InversionParams);

cd(curdir)
copyfile('master.m',[runname, '/Parameters.txt']); 
