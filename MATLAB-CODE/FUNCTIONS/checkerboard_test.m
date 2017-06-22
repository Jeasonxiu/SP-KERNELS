function checkerboard_test(nu0,nIter,sigma_n)
    %
    % This code generates a synthetic waveform associated with a
    % checkerboard structure by loading in a set of kernels from a Matlab
    % file.  Random gaussian noise is then added to each time point, and
    % the result is inverted for structure using a strategy outlined in
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
    
    %nTimes=length(KTimes);
    nTimes=1000;
    
    nDep=length(Scat_Depths);
    nSta=length(Stations);
    
    G=zeros(nTimes,nDep*nSta);
    d=zeros(nTimes,1);
    
    %Here we make a matrix C nSta by nDep which store icol values
    C=zeros(nSta,nDep);
    icol=0;
    for ista = 1:1:nSta;
       for idep = 1:1:nDep;
            icol=icol+1;
            C(ista,idep)=icol;
       end
    end
    
    %R is the regularization matrix
    R=zeros(nSta*nDep,nSta*nDep);
    
    volume_input=zeros(nDep,nSta);
    for ista = 1:1:nSta;
       for idep = 1:1:nDep;
            icol=C(ista,idep);
            G(1:nTimes,icol)=Kernel(idep,ista,1:nTimes);
            %Make checkerboard
            %if (mod(ista,10) >= 5 && mod(idep,10) >= 5) || ...
            %   (mod(ista,10) < 5 && mod(idep,10) < 5);
            %ss
            %Make interface
            if (idep<nDep/5.0)
                tmp=squeeze(Kernel(idep,ista,1:nTimes)); %
                d(1:nTimes,1)=d+tmp;
                volume_input(idep,ista)=1;
            else 
                %tmp=squeeze(Kernel(idep,ista,:));
                %d(:,1)=d+tmp;
                volume_input(idep,ista)=0.0;
            end
            
            %Construct regularization matrix
            if idep > 1 && idep < nDep && ista > 1 && ista < nSta;        
                icolu=C(ista,idep+1);
                icold=C(ista,idep-1);
                icoll=C(ista+1,idep);
                icolr=C(ista-1,idep);
                R(icol,icol)  = 1.0;
                R(icol,icolu) =-0.25;
                R(icol,icold) =-0.25;
                R(icol,icoll) =-0.25;
                R(icol,icolr) =-0.25;
            elseif idep ==1 && ista > 1 && ista < nSta;
                icolu=C(ista,idep+1);
                icoll=C(ista+1,idep);
                icolr=C(ista-1,idep);
                R(icol,icol)  = 1.0;
                R(icol,icolu) =-1.0/3.0;
                R(icol,icoll) =-1.0/3.0;
                R(icol,icolr) =-1.0/3.0;
            elseif idep ==nDep && ista > 1 && ista < nSta;
                icold=C(ista,idep-1);
                icoll=C(ista+1,idep);
                icolr=C(ista-1,idep);
                R(icol,icol)  = 1.0;
                R(icol,icold) =-1.0/3.0;
                R(icol,icoll) =-1.0/3.0;
                R(icol,icolr) =-1.0/3.0;
            elseif idep > 1 && idep < nDep && ista == 1;        
                icolu=C(ista,idep+1);
                icold=C(ista,idep-1);
                icoll=C(ista+1,idep);
                R(icol,icol)  = 1.0;
                R(icol,icolu) =-1.0/3.0;
                R(icol,icold) =-1.0/3.0;
                R(icol,icoll) =-1.0/3.0;
            elseif idep > 1 && idep < nDep && ista == nSta;        
                icolu=C(ista,idep+1);
                icold=C(ista,idep-1);
                icolr=C(ista-1,idep);
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
    
    %***
    %Uncomment if you want to force it to minimize the model norm
    %disp('Minimizing model norm')
    %R=diag(ones(nSta*nDep,1));
    %***
    
    %Add noise to data
    %sigma_n=0.0001;
    d=d+randn(size(d))*sigma_n;
    
    SIGMA_INV=eye(length(d))*sigma_n.^-1;
    
    dhat=SIGMA_INV*d;
    B=SIGMA_INV*G;
    
    %Add constraints of minimal model size
    
    nu=nu0;
    
    %Set T as the expected value of error norm, i.e, sqrt(N)
    TSquared=length(d);
    
    for iter = 1:nIter;
    
    %Invert system one
    Bf=horzcat(B',R*nu.^-0.5)';
    df=[dhat',zeros(1,nSta*nDep)]';
    m=Bf\df;

    E=R*m * nu^(-1.5);
    
    %Invert system two
    TMP=[zeros(1,length(dhat)), E' ]';
    
    dmdnu=Bf\TMP;
    
    F=sum((dhat-B*m).^2);
    dFdnu= -2/nu .* (R' * R * m)' * dmdnu;
    
    icol=0;
    volume=zeros(nDep,nSta);
    for ista = 1:1:nSta;
        for idep = 1:1:nDep;
           icol=icol+1;
            volume(idep,ista)=m(icol); 
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
    
    figure(1);subplot(2,2,1);
    pcolor(Stations,-Scat_Depths,fliplr(volume_input));
    shading flat
    title('Input Structure')
    colorbar;
    
    figure(1);subplot(2,2,[3 4]);
    dt=0.04;
    time=(0:length(dhat)-1) * dt;
    time=-time;
    plot(time,dhat *sigma_n)
    hold on
    plot(time,B*m *sigma_n)
    title('Daughter Waveforms')
    
   print('checkerboard','-dpng','-r300');
    
end


