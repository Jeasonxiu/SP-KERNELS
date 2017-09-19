classdef FigureInversion < MyFigure 
    properties
        clabel='Scattering Potential'
        path
        labdepth
        labwavlen
        labamp
    end
    methods
        function obj=FigureInversion(...
                path,labdepth,labwavlen,labamp,savename)
            
            obj.labdepth=labdepth;
            obj.labwavlen=labwavlen;
            obj.labamp=labamp;
            obj.path=path;
            labprops.wavlen=obj.labwavlen;
            labprops.amplitude=obj.labamp;
            labprops.depth=obj.labdepth;
            Inversion=ReadFromDisk(obj);
            obj.fig=plot_model(...
                Inversion.VelocityModel2D(),savename,obj.clabel,labprops);
            
            %[therat,medrat,stdrat] = plot_recovery_performance(Inversion.VelocityModel2D());
            
            %fprintf('%.3f   %.3f+/-%.3f\n',therat, medrat, stdrat)
           
        end

    end
    methods (Access = private)
        function Inversion=ReadFromDisk(obj)
            Inversion=[];
            load(sprintf('%s/I',obj.path));
        end
    end    
end
