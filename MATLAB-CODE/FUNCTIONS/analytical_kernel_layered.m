% Calcuates Sp kernels for a layercake background model.
% --N.J. Mancinelli 
%
function [Kernel,X,Y,Z,KTimes,nTimes] = analytical_kernel_layered(PdirectList, xs, zs, tchar, nderiv, model)
    %Input: PdirectList - list of ray parameters (s/km) to use in calculation
    %       tchar - characteristic time (s)
    %       nderiv - the number of derivatives to take on source-time function
    %
    %Outputs:
    
    PdirectList=PdirectList*111.1;
    
    %model.hs=[60.0,60.0,180.0];
    %model.vp=[5.660,7.920,7.280];
    %model.vs=[3.2,4.4,4.05];
    
    %Initialize matrix
    KTimes=-0:-0.4:-30;
    nTimes=length(KTimes);
    Kernel=zeros(length(zs),length(xs),nTimes,length(PdirectList));
    
    iPdirect=0;
    for Pdirect = PdirectList;
        iPdirect=iPdirect+1;
        
        [Tdirect,Tscat,Pscat,xs,zs]=read_in_necessary_files(Pdirect, xs, zs, model);

        %Step 1: Calculate array of times (Sp isochrons)

        % Calculate the src time function
        %tchar=0.8;  % affect amplitude and period of source-time function
        %deriv=3.0;
        %tchar=0.1;
        %nderiv=0.0;

        trial_times=-5:0.1:5;
        amps=zeros(length(trial_times),1);

        [G] = geom_spreading(Pscat,xs,zs,model);
        %[G] = geom_spreading_halfspace(xs,zs,0.5);

        [THE,~,~] = calc_angle(Pdirect,Pscat,xs,zs,model);
        vp=7.920;
        vs=4.400;
        rho=3.300;
        vpert=0.01;
        XSP = scattering_pattern(THE,vp,vs,rho,vpert);
        
        %calculate reference amplitude
        for itime=1:length(trial_times);
           amps(itime)=fractional_deriv(tchar,2.0,trial_times(itime)); 
        end
        
	    AmpRef=max(abs(amps));

        amps=amps/AmpRef;
        
        %calculate actual amplitude
        for itime=1:length(trial_times);
           amps(itime)=fractional_deriv(tchar,nderiv,trial_times(itime)); 
        end
        
        amps=amps/AmpRef;
        
        STFun.amplitude=amps;
        STFun.time=trial_times;
        save('STFun','STFun');
        
        [T] = timeshifts(Tdirect,Tscat,xs,zs);
        
        for itime = 1:length(KTimes);

            time=KTimes(itime);    

            %Step 2: Subtract target time from array of times

            DT=T-time;

            %Step 3: Function maps array of times to array of amplitudes (from src time function)

            A = time2amp(DT,trial_times,amps);

            %Step 4: Add factors for geometrical spreading / scattering pattern

            K=A .* G .* XSP;
            %K=A .* XSP;

            %Rotate onto daughter component
            [Krot] = rotate_kernel_to_daughter(K,Pdirect,Pscat,model);

            Kernel(:,:,itime,iPdirect)=Krot';

        end
    
    end
        
    [X,Y,Z] = meshgrid(xs,zs,1:nTimes);
    Kernel(isnan(Kernel))=0.0;
end

function [Krot] = rotate_kernel_to_daughter(K,Pdirect,Pscat,model)

Pdirect = Pdirect / 111.1;
Pscat = Pscat / 111.1;

%[alpha,beta]=get_velocity_from_profile('MIGRA/myvmod.nd',0);
[alpha,beta]=get_v(model,0);

THE=asind(Pscat*alpha);

Kz=K.*cosd(THE);  %vertical
Kr=K.*sind(THE);  %radial

%rotate to daughter using free-surface transform

%ray param of incident S wave
p=Pdirect;
%disp(p);
%incidence angle at 300 km
%[~,betatmp]=get_velocity_from_profile('MIGRA/myvmod.nd',300);
%incang=asind(p*betatmp);
%fprintf('Debug 2: %f %f %f\n',incang,betatmp,p)

qa0=sqrt(alpha^-2 - p.^2);

Vpz=-(1.-2.*beta^2*p.^2) ./ (2.*alpha*qa0);
Vpr=p*beta^2 / alpha;

M11=Vpz;
M12=Vpr;

%fprintf('alpha, beta, incangle, p = %f, %f, %f, %f\n',alpha,beta, incangle, p);
%fprintf('Ms: %f, %f\n',M11,M12)

%rotate to P component (daughter for S-to-P)
Krot = -M11 * Kz + M12 * Kr;

end