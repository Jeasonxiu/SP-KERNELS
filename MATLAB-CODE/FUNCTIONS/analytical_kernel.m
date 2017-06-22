function [Kernel,zs,xs,Angles,KTimes] = analytical_kernel()
    %Step 1: Calculate array of times (Sp isochrons)
    
    %incidence angle and halfspace velocities
    vp=7.92;
    vs=4.4;
    rho=3.3;
    vpert=0.01;
    
    %exponent for geometrical spreading (0.5 for 2-D)
    exp=0.5;
    
    KTimes=-0:-0.2:-30;
    Angles=[15,20,25];
    
    % set up nodes
    x1=-600.0;
    x2=150.0;
    dx=5.0335;
  
    z1=0.0;
    z2=300.0;
    dz=1.0;

    xs=x1:dx:x2;
    zs=z1:dz:z2;
    
    Kernel=zeros(length(zs),length(xs),length(Angles),length(KTimes));
    
    % Calculate the src time function
    tchar=1.6;  % affect amplitude and period of source-time function
    nderiv=3.0;

    trial_times=-10:0.1:10;
    amps=zeros(length(trial_times),1);

    for itime=1:length(trial_times);
       amps(itime)=fractional_deriv(tchar,nderiv,trial_times(itime)); 
    end
    
    for ithe = 1:length(Angles);
        theinc=Angles(ithe);
        [T,xs,zs] = calc_isochrons(xs,zs,theinc,vp,vs);
        for itime = 1:length(KTimes);
            
            time=KTimes(itime);    

            %Step 2: Subtract target time from array of times

            DT=T-time;

            %Step 3: Function maps array of times to array of amplitudes (from src time function)

            A = time2amp(DT,trial_times,amps);

            %Step 4: Add factors for geometrical spreading / scattering pattern

            G = geom_spreading(xs,zs,exp);

            THE = calc_angle(xs,zs,theinc);
            XSP = scattering_pattern(THE,vp,vs,rho,vpert);

            K=A.*G.* XSP;
           
            %Rotate onto daughter component
            [Krot] = rotate_kernel_to_daughter(K,xs,zs,vp,vs,theinc);

            Kernel(:,:,ithe,itime)=Krot;

        end
    end
    
end

function [Krot] = rotate_kernel_to_daughter(K,xs,zs,alpha,beta,incangle)

[THE] = calc_angle(xs,zs,0.0);

Kz=K.*cos(THE);  %vertical
Kr=K.*sin(THE);  %radial

%rotate to daughter using free-surface transform

%ray param of incident S wave
p=sind(incangle)/beta;

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

fac = rho * (dbeta_over_beta * 2 * beta/alpha);

%ARG = el0/gam0 * sqrt(1+gam0^2-2*gam0*cos(THE));

%Pterm = dbeta_over_beta * sqrt(exp(-ARG.^2*a^2 / 4));

XSP = fac * sin(2*THE); %.* Pterm;

end

function [THE] = calc_angle(xs,zs,theinc)

[XS,ZS]=meshgrid(xs,zs);

THE = -atan2(XS,ZS);

THE=THE-theinc*pi/180.0;

end 

function [G] = geom_spreading(xs,zs,n)
%Amp. goes as 1/sqrt(r) in 2-D
%  1/r in 3-D 

[XS,ZS]=meshgrid(xs,zs);

D=sqrt(XS.^2 + ZS.^2);
G=D.^-n;
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

function [T,xs,zs] = calc_isochrons(xs,zs,the,vp,vs)
%
% Plots isochron on axiss
%
[T,xs,zs]=timeshifts_halfspace(xs,zs,0,0,the,vp,vs);
%[T,xs,zs] = timeshifts_layercake(xs,zs);
end

function [T,xs,zs] = timeshifts_layercake(xs,zs)
[T10,xs0,zs0] = read_table('MIGRA/table_1_time.txt');
[T20,~,~] = read_table('MIGRA/table_2_time.txt');

[XS0, ZS0] = meshgrid(xs0,zs0);
[XS, ZS] = meshgrid(xs,zs);


T1=interp2(XS0,ZS0,T10',XS,ZS);
T2=interp2(XS0,ZS0,T20',XS,ZS);

tref=interp2(XS0,ZS0,T10',0,0);

T = T1+T2 - tref;
end

%%Above funcitons for halfspace

function [A,xs,zs] = read_table(tablename)

fid=fopen(tablename);
[~] = fgetl(fid);
dims =str2num(fgetl(fid));

nz=dims(1);
nx=dims(2);

tmp=fgetl(fid);
zs=str2num(tmp);

xs=zeros(nx,1);
A=zeros(nx,nz);

for iline=1:nx;
    tmp=fgetl(fid);
    nfo=str2num(tmp);
    xs(iline)=nfo(1);
    A(iline,:)=nfo(2:end);
end

fclose(fid);

end

function [T,xs,zs] = timeshifts_halfspace(xs,zs,xrec,zrec,the,vp,vs)

  [XS,ZS]=meshgrid(xs,zs);
  
% calc traveltime for S-wave
  P1=[0.0 0.0];
  P2=[cosd(the) sind(the)];
  NUM= (P2(2)-P1(2))*XS - (P2(1)-P1(1))*ZS + P2(1)*P1(2) - P2(2)*P1(1);
  DEN=((P2(2)-P1(2)).^2 + (P2(1)-P1(1)).^2 )^0.5;
  
  T1=(NUM/DEN)/vs;
  
% calc time from node pt to scatterer

  T2=(sqrt((XS-xrec).^2+(ZS-zrec).^2))/vp;
  
  trec=interp2(XS,ZS,T1,xrec,zrec);
  
  T=T1+T2 - trec;

end

function [THE,xs,zs] = compute_scang()
    [P1,~,~] = read_table('MIGRA/table_1_rayp.txt');
    [P2,xs,zs] = read_table('MIGRA/table_2_rayp.txt');
    
    P1=P1 / 111.1;
    P2=P2 / 111.1;
    
    %u
    [vps,vss]=get_velocity_from_profile('MIGRA/myvmod.nd',zs);
    
    [XS,~]=meshgrid(xs,zs);
    [VPS,~]=meshgrid(vps,xs);
    [VSS,~]=meshgrid(vss,xs);
    
    THE1=asind(P1.*VSS);
    THE2=asind(P2.*VPS).*sign(-XS)';

    THE=THE1-THE2;
    
end

function [THE,xs,zs] = compute_dpdx()
%%
[P,xs,zs] = read_table('MIGRA/table_2_rayp.txt');
P = P / 111.1;

[vps,vss]=get_velocity_from_profile('MIGRA/myvmod.nd',zs);

[VPS,XS]=meshgrid(vps,xs);

U1=1./VPS;
THE1=asind(P./U1);
THE2=asind(P/U1(1,1));

%take derivative and resample
tmp=diff(P.*sign(XS));
xmid=(xs(1:(end-1))+xs(2:end))/2;
dpdx=interp2(xmid,zs,tmp',xs,zs);


%This is an ad hoc expression modified from Eq 6.23 in Shearers Intro to
%Seismoogy
energy=(1/2/pi)*abs(dpdx)'./cosd(THE1)./cosd(THE2)./U1;

amp=(sqrt(energy));

figure(1)
subplot(2,1,1)

contourf(xs,-zs,log10(amp')); shading flat; colorbar

subplot(2,1,2)

[XS,ZS]=meshgrid(xs,zs);
D=sqrt(XS.^2 + ZS.^2);
contourf(xs,-zs,log10(1./sqrt(D))); shading flat; colorbar
%%
end