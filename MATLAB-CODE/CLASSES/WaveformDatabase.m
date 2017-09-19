classdef WaveformDatabase
    properties
        filelist
        nfiles
        nfiles_max
    end
    
    methods
        function obj=WaveformDatabase(varargin)
            if nargin>0
                obj.nfiles_max=varargin{1};
            end
            
        end
        
        function obj=populate_filelist(obj)
            obj.nfiles=0;
            cdw;
            cd KERNEL-SEM
            directories=dir('OUTPUT_FILES*');
            for idir=1:length(directories);
                directory=directories(idir).name;
                tmpstr=sprintf('%s/OUTPUT_FILES/*semd',directory);
                files=dir(tmpstr);
                for ifil=1:length(files)
                    file=files(ifil).name;
                    fullpath=sprintf('KERNEL-SEM/%s/OUTPUT_FILES/%s',directory,file);
                    obj.nfiles=obj.nfiles+1;
                    obj.filelist{obj.nfiles}=fullpath;
                    if obj.nfiles == obj.nfiles_max;
                        cdw;
                        return
                    end
                end
            end
            cdw;

        end
        
        function obj=decimate_files(obj)
            function Wf=load_file(path)
                seis=load(path);
                tstart=seis(1,1);
                dt=seis(2,1)-tstart;
                Wf=Waveform(seis(:,2),tstart,dt);
            end
            
            function write_file(path,Wf)
                fid=fopen(path,'w');
                for ipt =1:Wf.npts                   
                    fprintf(fid,'%8.3f   %E \n', Wf.t(ipt), Wf.u(ipt));
                end
                fclose(fid);
                
                
            end
            
            for ifil = 1:obj.nfiles;
                Wf=load_file(obj.filelist{ifil});
                if abs(Wf.dt - 0.02)>10^-5; continue; end
                
                Wf_dec=decimate(Wf,10);
                
                
                write_file('dummy',Wf_dec)
                
                
            end
            
            
        end
        
        
        
    end
    
    
end