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
    
    iPdirect=0;
    for Pdirect = PdirectList;
        iPdirect=iPdirect+1;
        
        [Tdirect,Tscat,Pscat,xs,zs]=read_in_necessary_files(Pdirect, xs, zs, model);

        %Step 1: Calculate array of times (Sp isochrons)

        KTimes=-0:-0.2:-30;
        nTimes=length(KTimes);

        Angles=999;

        % Calculate the src time function
        %tchar=0.8;  % affect amplitude and period of source-time function
        %deriv=3.0;
        %tchar=0.1;
        %nderiv=0.0;

        trial_times=-10:0.1:10;
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
        
        for ithe = 1:length(Angles);
            [T,xs,zs] = calc_isochrons(Tdirect,Tscat,xs,zs);
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

function [XSP] = scattering_pattern(THE,alpha,beta,rho,dbeta_over_beta)
%See Sato et al. (2012), Eqn. 4.62 or Rondenay (2009) or Bostock and
%Rondenay (1999)

%gam0=sqrt(3);
%a=5;
%T=2.0;
%el0=2*pi/(beta*T);

%fprintf(' CALCULATING SCATTERING PATTERN\n\n')
%fprintf(' dbeta/beta, beta, alpha, rho = %.3f  %.3f  %.3f  %.3f\n\n',dbeta_over_beta, beta, alpha, rho)

fac = -rho * (dbeta_over_beta * 2 * beta/alpha);

%fprintf(' fac = %f\n\n', fac)

%ARG = el0/gam0 * sqrt(1+gam0^2-2*gam0*cos(THE));

%Pterm = dbeta_over_beta * sqrt(exp(-ARG.^2*a^2 / 4));

XSP = fac * sind(2*THE); %.* Pterm;

%Inspired by eqn 7 of Rondenay, Bostock, and Fischer
%1/(2pi) rather than 1/(4pi) because this is 2-D migration rather than 3d

XSP = XSP * 1/2/pi;

end

function [A] = time2amp(DT,trial_times,amps)

%Simple gaussian 
extrapval=0.0;

A=interp1(trial_times,amps,DT,'pchip',extrapval);

%A = exp(-DT.^2./(2.0*tchar^2));

%Ricker wavelet
%A = (1-(DT./tchar).^2) .* A;

%Deriv of Ricker wavelet
%A = -2*(tchar-DT)/tchar^2 .* A + (1-(DT/tchar).^2) .* (-DT/tchar^2) .* A;

%A = -A .* (DT.^3 - 3*tchar^2 *DT)/tchar^4;

end

function [T,xs,zs] = calc_isochrons(Tdirect,Tscat,xs,zs)
%
[T,xs,zs] = timeshifts_layercake(Tdirect,Tscat,xs,zs);
end

function [T,xs0,zs0] = timeshifts_layercake(Tdirect,Tscat,xs0,zs0)
T10=Tdirect;
T20=Tscat;

[XS0, ZS0] = meshgrid(xs0,zs0);
%[XS, ZS] = meshgrid(xs,zs);


%T1=interp2(XS0,ZS0,T10',XS,ZS);
%T2=interp2(XS0,ZS0,T20',XS,ZS);

tref=interp2(XS0,ZS0,T10,0,0);

T = T10'+T20 - tref;
end

function [THE,xs,zs] = calc_angle(Pdirect,Pscat,xs,zs,model)
    P1=Pdirect / 111.1;
    P2=Pscat / 111.1;
    
    %u
    %[vps,vss]=get_velocity_from_profile('MIGRA/myvmod.nd',zs);
    
    [vps,vss]=get_v(model,zs);
    
    [XS,~]=meshgrid(xs,zs);
    [VPS,~]=meshgrid(vps,xs);
    [VSS,~]=meshgrid(vss,xs);
    
    THE1=asind(P1.*VSS);
    THE2=asind(P2.*VPS);

    THE=THE1-THE2;
    
end

function [amp] = geom_spreading(Pscat,xs,zs,model)
%%
P = Pscat / 111.1;

%[vps,~]=get_velocity_from_profile('MIGRA/myvmod.nd',zs);
[vps,~]=get_v(model,zs);

[VPS,~]=meshgrid(vps,xs);

U1=1./VPS;
THE1=asind(P./U1);
THE2=asind(P/U1(1,1));

%take derivative and resample
tmp=diff(P);
xmid=(xs(1:(end-1))+xs(2:end))/2;
dpdx=interp2(xmid,zs,tmp',xs,zs);


%This is an ad hoc expression modified from Eq 6.23 in Shearers Intro to
%Seismoogy
energy=(1/2/pi)*abs(dpdx)'./cosd(THE1)./cosd(THE2)./U1;

amp=(sqrt(energy));

%figure(1)
%subplot(2,1,1)
%
%contourf(xs,-zs,log10(amp')); shading flat; colorbar
%
%subplot(2,1,2)
%
%[XS,ZS]=meshgrid(xs,zs);
%D=sqrt(XS.^2 + ZS.^2);
%contourf(xs,-zs,log10(1./sqrt(D))); shading flat; colorbar
%%
end

function [G] = geom_spreading_halfspace(xs,zs,n)
%Amp. goes as 1/sqrt(r) in 2-D
%  1/r in 3-D 

[XS,ZS]=meshgrid(xs,zs);

D=sqrt(XS.^2 + ZS.^2);
G=D.^-n;
G=G';
end
