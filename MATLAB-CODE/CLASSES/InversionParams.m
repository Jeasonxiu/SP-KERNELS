classdef InversionParams
    properties
        skipSta=10;
        tDeci=20;
        dxin=25;
        dzin=5;
        Nus=3.15E-1;
        Norm_Opts=2;
        nIterMax=500;
        direction=1; %3 (synthetics from both directions)
        Kernel_Type=3;  %3 (analytical) should be best
        ImagingMethod=1; %1 for Back Proj, %2 for CG Inversion
        saveFilename='model';
    end
    methods
        function obj=InversionParams()
        end
        
    end
    
end