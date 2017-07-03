function [Kernel,X,Y,Z,KTimes,nTimes] = load_kernel(path,tDeci)

    Kernel=[];
    KTimes=[];
    
    str1=[path '/Kernel_Angles_x2_0.01_2500.mat'];
    str2=[path '/stalocs.txt'];
    
    load(str1)
    stalocs = load(str2)/1000.0 - 1500.0; %convert to km

    %Decimate Kernels to speed up inversions
    %ip=1; %just use first angle calc for now
    
    %Going to time zero appears to work
    it1=1;
    
    Kernel=squeeze(Kernel(:,:,:,it1:tDeci:end));
    KTimes=KTimes(it1:tDeci:end);
    nTimes=length(KTimes);
    %Kernel=cumtrapz(Kernel,3);
    
    [X,Y,Z] = meshgrid(-stalocs,Scat_Depths,1:nTimes);

end