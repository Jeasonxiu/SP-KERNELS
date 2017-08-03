classdef InversionParams
    properties
        tDeci
        dxin
        dzin
        Nus
        Norm_Opts
        nIterMax
        direction      %3 (synthetics from both directions)
        Kernel_Type    %3 (analytical) should be best
        ImagingMethod  %1 for Back Proj, %2 for CG Inversion
        TakeDifferences
        DeconvolveParentWaveform
        saveFilename='model';
    end
    methods
        function obj=InversionParams()
        end
        function obj=SetDefaultParams1(obj)
            obj.dxin=5;
            obj.dzin=2;
            obj.Nus=3.15E-1;
            obj.Norm_Opts=2;
            obj.nIterMax=500;
            obj.direction=3; %3 (synthetics from both directions)
            obj.Kernel_Type=3;  %3 (analytical) should be best
            obj.ImagingMethod=1; %1 for Back Proj, %2 for CG Inversion
            obj.TakeDifferences=false;
            obj.DeconvolveParentWaveform=false;
        end      
        function obj=SetDefaultParams2(obj)
            obj.dxin=5;
            obj.dzin=2;
            obj.Nus=3.15E-1;
            obj.Norm_Opts=2;
            obj.nIterMax=500;
            obj.direction=3; %3 (synthetics from both directions)
            obj.Kernel_Type=3;  %3 (analytical) should be best
            obj.ImagingMethod=2; %1 for Back Proj, %2 for CG Inversion
            obj.TakeDifferences=true;
            obj.DeconvolveParentWaveform=false;
        end
    end
    
end
