% call invert_synthetics(1.0e-6,5,0.00005) to get a good result inverting
%  for kernels. Nick Mancinelli - September 2016.
%
%
%

function invert_synthetics_l1(tDeci,label)
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
    %
    %  Nick Mancinelli - September, 2016
    
    load Kernel_Plus_Weaker

    %Decimate Kernels to speed up inversions
    Kernel=Kernel(:,:,1:tDeci:end);
    KTimes=KTimes(1:tDeci:end);
    nTimes=length(KTimes);
    
    nSta=99;
    sta1=1;
    
    d=[];
    for iseis = sta1:(sta1+nSta-1);
        Daughter=load(sprintf('~/PROJECTS/SP-RF-SYN/CASE4/OUTPUT_FILES_50000-40000-1/OUTPUT_FILES/AA.S00%02d.BXP.semd',iseis));
        Parent=load(sprintf('~/PROJECTS/SP-RF-SYN/CASE4/OUTPUT_FILES_50000-40000-1/OUTPUT_FILES/AA.S00%02d.BXS.semd',iseis));
        %Daughter=load(sprintf('~/PROJECTS/SP-RF-SYN/CASE1/OUTPUT_FILES_45000-160000-1/OUTPUT_FILES/AA.S00%02d.BXP.semd',iseis));
        %Parent=load(sprintf('~/PROJECTS/SP-RF-SYN/CASE1/OUTPUT_FILES_45000-160000-1/OUTPUT_FILES/AA.S00%02d.BXS.semd',iseis));
        %Daughter=load(sprintf('~/PROJECTS/SP-RF-SYN/KERNEL/OUTPUT_FILES_150000-0.01-1-24/OUTPUT_FILES/AA.S00%02d.BXP.semd',30));
        %Parent=load(sprintf('~/PROJECTS/SP-RF-SYN/KERNEL/OUTPUT_FILES_150000-0.01-1-24/OUTPUT_FILES/AA.S00%02d.BXS.semd',30));
        [~,imax]=max(abs(Parent(:,2)));
        time=Parent(:,1) - Parent(imax,1);
        %dtmp=Daughter(imax-nTimes*tDeci+1:imax,2);
        
        dtmp=interp1(time,Daughter(:,2),KTimes)';
        %dtmp=Daughter(time>=min(KTimes) & time<=max(KTimes),2);
        
        %dtmp=flipud(dtmp(1:end));

        d=[dtmp; d];
    
    end

    nDep=length(Scat_Depths);
    nOff=length(Stations)*3;  %Triple width of model
    
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
    
    %volume_input=zeros(nDep,nOff);
    for ista = 1:nSta;
        for ioff = 1:1:nOff;
           for idep = 1:1:nDep;
            
            icol=C(ioff,idep);
            ioff_for_sta=ioff-(ista+1);
            if ioff_for_sta<1 || ioff_for_sta > length(Stations);
                G((1:nTimes)+(ista-1)*nTimes,icol)=zeros(1,nTimes);
            else
                G((1:nTimes)+(ista-1)*nTimes,icol)=Kernel(idep,ioff_for_sta,1:nTimes);
            end
         
           end      
       end
    end

    %d=G*m
    
    epsilon=0.1;
    
    m0=G'*d;
    m = l1qc_logbarrier(m0, G, [], d, epsilon, 1e-3);
    
    %Generate image from m
    icol=0;
    volume=zeros(nDep,nOff);
    for ioff = 1:1:nOff;
        for idep = 1:1:nDep;
           icol=icol+1;
            volume(idep,ioff)=m(icol); 
        end
    end
    
    
    figure(1);clf;subplot(2,2,2);
    pcolor(1:nOff,-Scat_Depths,fliplr(volume));
    %title(sprintf('Nu = %e ',nu));
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
    time=(0:length(d)-1) * dt;
    time=-time;
    plot(time,d)
    hold on
    plot(time,G*m)
    title('Daughter Waveforms')
    
    tmpstr=sprintf('recovery_test_%s',label);
    
    print(tmpstr,'-dpng','-r300');
   
    %Force data to equal something
    %mtest=m*0.0;
    %icol=C(ceil(nOff/2),ceil(nDep/2));
    %mtest(icol)=1.0;
    %dhat_test=B*mtest;
    
end


