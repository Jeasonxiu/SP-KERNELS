classdef DataParams
    properties
        lab_amplitude
        lab_wavelength
        lab_depth
        skipSta
        KTimes
        Angles
        TakeDifferences
        DeconvolveParentWaveform
        direction
    end
   methods
       function obj=DataParams(lab_amplitude,lab_wavelength,lab_depth,skipSta,KTimes,Angles,TakeDifferences,DeconvolveParentWaveform,direction)
           obj.lab_amplitude=lab_amplitude;
           obj.lab_wavelength=lab_wavelength;
           obj.lab_depth=lab_depth;
           obj.skipSta=skipSta;
           obj.KTimes=KTimes;
           obj.Angles=Angles;
           obj.TakeDifferences=TakeDifferences;
           obj.DeconvolveParentWaveform=DeconvolveParentWaveform;
           obj.direction=direction;
       end
       
   end
    
end