% A master script to perform Sp migration or inversion.
% N. J. Mancinelli -- June 2017
%
runname=sprintf('OUTPUT/%s/SPACING_TEST_SNG',date);

InversionParams.skipSta=5;
InversionParams.tDeci=5;
InversionParams.dxin=5;
InversionParams.dzin=2;
InversionParams.Nus=[3.15E+1,3.15E+0,3.15E-1,3.15E-2];
InversionParams.Norm_Opts=[2];
InversionParams.nIterMax=500;
InversionParams.direction=3;
InversionParams.Kernel_Type=3;
InversionParams.ImagingMethod=2;
InversionParams.saveFilename='model';

curdir=pwd;

mkdir(runname)
cd(runname)

DirsToLink = {'MATLAB-CODE/FUNCTIONS','MIGRA','../../TEST_MODELS'};

for Dir = DirsToLink;
    cmd=sprintf('ln -sf %s/%s .', curdir, Dir{1});
    system(cmd);
end

addpath('FUNCTIONS');

back_project_synthetics(InversionParams);

cd(curdir)
copyfile('master.m',[runname, '/Parameters.txt']); 