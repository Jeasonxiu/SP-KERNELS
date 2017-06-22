function linearity_test()
    close all
    [ktimes,kernel1]=load_kernel('Kernel_Angles_x2_0.01_2500.mat');
    [~,kernel2]=load_kernel('Kernel_Angles_x2_0.02_2500.mat');
  
    slice1 = slice_kernel(kernel1,10,50,1);
    slice2 = slice_kernel(kernel2,10,50,1);
    
    l1=plot(ktimes,slice1,'-k','LineWidth',2);
    hold on;
    l2=plot(ktimes,slice2*0.5,'--r','LineWidth',2);
    l3=plot(ktimes,slice2*0.5-slice1,'LineWidth',1);
    
    legend([l1,l2],'1%','2%')
    title('Slice through kernels at A, B, and C')
    xlabel('Time (s)')
    ylabel('Daughter Amplitude (Scaled)')
    
end

function slice=slice_kernel(kernel,ii,jj,kk)
    slice=squeeze(kernel(ii,jj,kk,:));
end

function [KTimes,Kernel]=load_kernel(filename)
        load(filename)
end