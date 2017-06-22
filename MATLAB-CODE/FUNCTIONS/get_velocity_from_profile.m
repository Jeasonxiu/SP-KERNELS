function [vptargs,vstargs]=get_velocity_from_profile(path2file,deptargs)
    fid=fopen(path2file);
    i=0;
    dp=zeros(100,1);
    vp=zeros(100,1);
    vs=zeros(100,1);
    
    while 1 
      line=fgetl(fid);
      if ~ischar(line), break, end
      nfo=str2num(line);
      if isempty(nfo), continue, end
      dp(i+1)=nfo(1);
      vp(i+1)=nfo(2);
      vs(i+1)=nfo(3);
      i=i+1;
    end
    
    fclose(fid);
    
    %resample at 1 km intervals
    dp2=min(deptargs):max(deptargs);
    vp2=min(deptargs):max(deptargs)*0.0;
    vs2=min(deptargs):max(deptargs)*0.0;
    
    for jj =1:length(dp2)
    for ii =1:length(dp);
        if dp2(jj) == dp(ii);
            vp2(jj)=vp(ii);
            vs2(jj)=vs(ii);
            break
        elseif dp2(jj) < dp(ii);
            ind1=ii-1;
            ind2=ii;
            vp2(jj) = interp1([dp(ind1),dp(ind2)],[vp(ind1),vp(ind2)],dp2(jj));
            vs2(jj) = interp1([dp(ind1),dp(ind2)],[vs(ind1),vs(ind2)],dp2(jj));
            break
        end
    end
    end
    
    if length(deptargs)==1;
       vptargs=vp2;
       vstargs=vs2;
       
       return
    end
    
    vptargs=interp1(dp2,vp2,deptargs);
    vstargs=interp1(dp2,vs2,deptargs);

    return
    
end