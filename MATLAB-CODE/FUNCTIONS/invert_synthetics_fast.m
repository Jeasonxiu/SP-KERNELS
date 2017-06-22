% call invert_synthetics(1.0e-6,5,0.00005) to get a good result inverting
%  for kernels. Nick Mancinelli - September 2016.
%
%
% This code loads synthetic waveforms from a file,
% and invers it for structure using a strategy outlined in
% Chapter 3 of Bob Parker's inverse theory book.
%
% Inputs:
%  nu0 -- initial value of the tradeoff parameter.  If subsequent
%  iterations produce a negative nu, try making nu0 smaller and start
%  again.
%  nIter -- number of interations
%  sigma_n -- noise level (e.g., 0.0001)
%  tdeci -- factor by which to decimate time dimension (e.g., 10)
%  dxin,dzin -- spacing between model nodes (in km, e.g. 5, 5)
%
%  Nick Mancinelli - September, 2016
function invert_synthetics_fast(nu0,nIter,sigma_n,tDeci,dxin,dzin,norm_opt,label)

    load('Kernel_Angles_x2_0.02_2500.mat')
    stalocs = load('stalocs.txt')/1000.0; %convert to km

    %Decimate Kernels to speed up inversions
    ip=1; %just use first angle calc for now
    Kernel=squeeze(Kernel(:,:,ip,1:tDeci:end));
    KTimes=KTimes(1:tDeci:end);
    nTimes=length(KTimes);
    
    Sta1=1;
    Sta2=150;
    skipSta=1;
    Stas=Sta1:skipSta:Sta2;
    
    nSta=length(Stas);
    
    d=[];
    for iseis = Stas;
        Daughter=load(sprintf('~/data2/nmancine/nmancine/SP-RF-SYN/KERNEL_GMSH/KERNEL/TWOLAY.4s/OUTPUT_FILES_0.05-23/OUTPUT_FILES/AA.S0%03d.BXP.semd',iseis));
        Parent=load(sprintf('~/data2/nmancine/nmancine/SP-RF-SYN/KERNEL_GMSH/KERNEL/TWOLAY.4s/OUTPUT_FILES_0.05-23/OUTPUT_FILES/AA.S0%03d.BXS.semd',iseis));
        [~,imax]=max(abs(Parent(:,2)));
        time=Parent(:,1) - Parent(imax,1);
        
        dtmp=interp1(time,Daughter(:,2),KTimes)';
        d=[dtmp; d];
    
    end
    
    !cat ~/data2/nmancine/nmancine/SP-RF-SYN/KERNEL_GMSH/KERNEL/TWOLAY.4s/OUTPUT_FILES_0.05-23/DATA/STATIONS | awk {'print $3'} > tmp.txt
    stalocs2=load('tmp.txt')/1000.0; % convert to km

    %nDep=length(Scat_Depths);
    %nOff=length(Stations)*2;  %Extend width of model
    
    %Set model parameters
    x1=min(stalocs2)*0.85;
    x2=max(stalocs2)*1.02;
    dx=dxin;
    
    z1=0;
    z2=150;
    dz=dzin;
    
    xs=x1:dx:x2;
    zs=z1:dz:z2;
    
    nDep=length(zs);
    nOff=length(xs);
    
    G=zeros(nTimes*nSta,nOff,nDep);
    %d=zeros(nTimes,1);
    
    %Here we make a matrix C nOff by nDep which store icol values
    C=zeros(nOff,nDep);
    icol=0;
    for ioff = 1:1:nOff;
       for idep = 1:1:nDep;
            icol=icol+1;
            C(ioff,idep)=icol;
       end
    end
    
    %R is the regularization matrix
    R=zeros(nOff*nDep,nOff*nDep);
    [X,Y,Z] = meshgrid(-stalocs+1500,Scat_Depths,1:nTimes);
    %volume_input=zeros(nDep,nOff);
    fprintf('Generating G and R matrices:\n');
    for ista = 1:nSta;
        
        fprintf('%.2f pct complete... \n',(ista-1)/nSta*100.0);   
                
        [Xq,Yq,Zq]=meshgrid(xs-stalocs2(ista),zs,1:nTimes);
        tmp=interp3(X,Y,Z,Kernel(:,:,1:nTimes),Xq,Yq,Zq,'spline',0.0);
        %tmp=interp3(X,Z,Y,Kernel(:,1:nTimes,:),Xq,Zq,Yq,'spline',0.0);
        G((1:nTimes)+(ista-1)*nTimes,:,:)=permute(tmp,[3 2 1]);

        for ioff = 1:1:nOff;
           for idep = 1:1:nDep;

            icol=C(ioff,idep);
            %ioff_for_sta=ioff-(Stas(ista)+1);

            %Construct regularization matrix
            if idep > 1 && idep < nDep && ioff > 1 && ioff < nOff;        
                icolu=C(ioff,idep+1);
                icold=C(ioff,idep-1);
                icoll=C(ioff+1,idep);
                icolr=C(ioff-1,idep);
                R(icol,icol)  = 1.0;
                R(icol,icolu) =-0.25;
                R(icol,icold) =-0.25;
                R(icol,icoll) =-0.25;
                R(icol,icolr) =-0.25;
            elseif idep ==1 && ioff > 1 && ioff < nOff;
                icolu=C(ioff,idep+1);
                icoll=C(ioff+1,idep);
                icolr=C(ioff-1,idep);
                R(icol,icol)  = 1.0;
                R(icol,icolu) =-1.0/3.0;
                R(icol,icoll) =-1.0/3.0;
                R(icol,icolr) =-1.0/3.0;
            elseif idep ==nDep && ioff > 1 && ioff < nOff;
                icold=C(ioff,idep-1);
                icoll=C(ioff+1,idep);
                icolr=C(ioff-1,idep);
                R(icol,icol)  = 1.0;
                R(icol,icold) =-1.0/3.0;
                R(icol,icoll) =-1.0/3.0;
                R(icol,icolr) =-1.0/3.0;
            elseif idep > 1 && idep < nDep && ioff == 1;
                icolu=C(ioff,idep+1);
                icold=C(ioff,idep-1);
                icoll=C(ioff+1,idep);
                R(icol,icol)  = 1.0;
                R(icol,icolu) =-1.0/3.0;
                R(icol,icold) =-1.0/3.0;
                R(icol,icoll) =-1.0/3.0;
            elseif idep > 1 && idep < nDep && ioff == nOff;
                icolu=C(ioff,idep+1);
                icold=C(ioff,idep-1);
                icolr=C(ioff-1,idep);
                R(icol,icol)  = 1.0;
                R(icol,icolu) =-1.0/3.0;
                R(icol,icold) =-1.0/3.0;
                R(icol,icolr) =-1.0/3.0;
            else
                %Don't do anything to the corners for now...
                %0.0 for nothing, 1.0 for small
                R(icol,icol) = 1.0;
            end
           end
       end            
    end
    
    G=permute(G,[1,3,2]);
    G=reshape(G,nTimes*nSta,nOff*nDep);
    
    %print_G(G);
    
    %***
    %Uncomment if you want to force it to minimize the model norm
    if norm_opt==1;
       disp('Minimizing model norm');R=diag(ones(nOff*nDep,1));
    else
       disp('Minimizing model roughness');
    end
    %***
    
    %Add noise to data
    %sigma_n=0.0001;
    %d=d+randn(size(d))*sigma_n;
    
    SIGMA_INV=eye(length(d))*sigma_n.^-1;
    
    dhat=SIGMA_INV*d;
    B=SIGMA_INV*G;
    
    %Add constraints of minimal model size
    
    nu=nu0;
    
    %Set T as the expected value of error norm, i.e, sqrt(N)
    TSquared=length(d);
    
    for iter = 1:nIter;
    
        %Invert system one
        Bf=(horzcat(B',R*nu.^-0.5)');
        df=[dhat',zeros(1,nOff*nDep)]';
        
        [dim1,dim2]=size(Bf);
        fprintf('Inverting system one, %d X %d matrix... ',dim1,dim2)
        tic; m=Bf\df; toc;
        fprintf('done!\n')
        
        if iter == nIter; break; end

        E=R*m * nu^(-1.5);

        %Invert system two
        
        TMP=[zeros(1,length(dhat)), E' ]';

        fprintf('Inverting system two, %d X %d matrix... ',dim1,dim2)
        dmdnu=Bf\TMP;
        fprintf('done!\n')

        F=sum((dhat-B*m).^2);
        dFdnu= -2/nu .* (R' * R * m)' * dmdnu;

        if abs(F-TSquared) < 0.1;
            break
        end

        fprintf('nu, F, T^2 = %e %f %f \n',nu , F, TSquared)

        nu=nu - (F-TSquared)/dFdnu;
    
    end
    
    icol=0;
    volume=zeros(nDep,nOff);
    for ioff = 1:1:nOff;
        for idep = 1:1:nDep;
           icol=icol+1;
            volume(idep,ioff)=m(icol); 
        end
    end
    
    %figure(1);clf;subplot(2,2,2);
    %pcolor(1:nOff,-Scat_Depths,fliplr(volume));
    %title(sprintf('Nu = %e ',nu));
    %shading flat
    %title('Inverted Structure')
    %colorbar;
    
    %figure(1);subplot(2,2,1);
    %pcolor(log10(abs(G))); shading flat;
    %shading flat
    %title('Input Structure')
    %colorbar;
    
    %figure(1);subplot(2,2,[3 4]);
    %dt=0.04;
    %time=(0:length(dhat)-1) * dt;
    %time=-time;
    %plot(time,dhat *sigma_n)
    %hold on
    %plot(time,B*m *sigma_n)
    %title('Daughter Waveforms')
    
    tmpstr=sprintf('recover_test_%s',label);
    
    %print('recovery_test','-dpng','-r300');
   
    %Force data to equal something
    %mtest=m*0.0;
    %icol=C(ceil(nOff/2),ceil(nDep/2));
    %mtest(icol)=1.0;
    %dhat_test=B*mtest;
    
    save(tmpstr,'nOff','volume','nu','dhat','sigma_n','B','m','label','nu0','nIter','sigma_n','tDeci','label','Scat_Depths','Stas','xs','zs','stalocs2','-v7.3');
    
end

function print_G(G)
	%
	figure(7192)
	pcolor(G); shading flat; colorbar;
	print('G.png','-dpng','-r300');
    clf;
    disp('G matrix exported to png file.')
	
end
