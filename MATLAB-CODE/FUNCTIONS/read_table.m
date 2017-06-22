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