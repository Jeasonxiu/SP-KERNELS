classdef DataParams
    properties
        path2seis
        path2ref
        skipSta
        KTimes
        Angles
        TakeDifferences
        DeconvolveParentWaveform
        Direction
    end
   methods
       function obj=DataParams(path2seis,path2ref,skipSta,KTimes,Angles,TakeDifferences,DeconvolveParentWaveform,Direction)
           obj.path2seis=path2seis;
           obj.path2ref=path2ref;
           obj.skipSta=skipSta;
           obj.KTimes=KTimes;
           obj.Angles=Angles;
           obj.TakeDifferences=TakeDifferences;
           obj.DeconvolveParentWaveform=DeconvolveParentWaveform;
           obj.Direction=Direction;
       end
       
   end
    
end