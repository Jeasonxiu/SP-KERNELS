% call invert_synthetics(1.0e-6,5,0.00005) to get a good result inverting
%  for kernels. Nick Mancinelli - September 2016.
%
%
%

function invert_synthetics(nu0,nIter,sigma_n)
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
    %
    %  Nick Mancinelli - September, 2016
    
    load Kernel_Plus_Weaker

    %nTimes=100;
    tDeci=10;
    Kernel=Kernel(:,:,1:tDeci:end);
    KTimes=KTimes(1:tDeci:end);
    nTimes=length(KTimes);
    
    nSta=75;
    sta1=1;
    
    d=[];
    for iseis = sta1:(sta1+nSta-1);
    
        Daughter=load(sprintf('~/PROJECTS/SP-RF-SYN/KERNEL/OUTPUT_FILES_150000-0.01-1-24/OUTPUT_FILES/AA.S00%02d.BXP.semd',iseis));
        Parent=load(sprintf('~/PROJECTS/SP-RF-SYN/KERNEL/OUTPUT_FILES_150000-0.01-1-24/OUTPUT_FILES/AA.S00%02d.BXS.semd',iseis));
        [~,imax]=max(abs(Parent(:,2)));
        time=Parent(:,1) - Parent(imax,1);
        %dtmp=Daughter(imax-nTimes*tDeci+1:imax,2);
        
        dtmp=interp1(time,Daughter(:,2),KTimes)';
        %dtmp=Daughter(time>=min(KTimes) & time<=max(KTimes),2);
        
        %dtmp=flipud(dtmp(1:end));

        d=[dtmp; d];
    
    end

    nDep=length(Scat_Depths);
    nOff=length(Stations);
    
    G=zeros(nTimes*nSta,nDep*nOff);
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
    
    %volume_input=zeros(nDep,nOff);
    for ista = 1:nSta;
        for ioff = 1:1:nOff;
           for idep = 1:1:nDep;
            
            icol=C(ioff,idep);
            ioff_for_sta=ioff+50-(ista+1);
            if ioff_for_sta<1 || ioff_for_sta > nOff;
                G((1:nTimes)+(ista-1)*nTimes,icol)=zeros(1,nTimes);
            else
                G((1:nTimes)+(ista-1)*nTimes,icol)=Kernel(idep,ioff_for_sta,1:nTimes);
            end
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
                %Ask corners to be small for now...
                R(icol,icol) = 1.0;
            end
           end      
       end
    end
    
    %***
    %Uncomment if you want to force it to minimize the model norm
    disp('Minimizing model norm');R=diag(ones(nOff*nDep,1));
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
        m=Bf\df;

        E=R*m * nu^(-1.5);

        %Invert system two
        TMP=[zeros(1,length(dhat)), E' ]';

        dmdnu=Bf\TMP;

        F=sum((dhat-B*m).^2);
        dFdnu= -2/nu .* (R' * R * m)' * dmdnu;

        icol=0;
        volume=zeros(nDep,nOff);
        for ioff = 1:1:nOff;
            for idep = 1:1:nDep;
               icol=icol+1;
                volume(idep,ioff)=m(icol); 
            end
        end

        if abs(F-TSquared) < 0.1;
            break
        end

        fprintf('nu, F, T^2 = %e %f %f \n',nu , F, TSquared)

        nu=nu - (F-TSquared)/dFdnu;
    
    end
    
    figure(1);clf;subplot(2,2,2);
    pcolor(Stations,-Scat_Depths,fliplr(volume));
    title(sprintf('Nu = %e ',nu));
    shading flat
    title('Inverted Structure')
    colorbar;
    
    %figure(1);subplot(2,2,1);
    %pcolor(log10(abs(G))); shading flat;
    %shading flat
    %title('Input Structure')
    %colorbar;
    
    figure(1);subplot(2,2,[3 4]);
    dt=0.04;
    time=(0:length(dhat)-1) * dt;
    time=-time;
    plot(time,dhat *sigma_n)
    hold on
    plot(time,B*m *sigma_n)
    title('Daughter Waveforms')
    
    
    
   print('checkerboard','-dpng','-r300');
   
    %Force data to equal something
    %mtest=m*0.0;
    %icol=C(ceil(nOff/2),ceil(nDep/2));
    %mtest(icol)=1.0;
    %dhat_test=B*mtest;
    
end


