function plot_geom_spreading(kernel,xs,zs)
[XS,ZS]=meshgrid(xs,zs);
D=sqrt(XS.^2+ZS.^2);

figure(2)
fun = @(x) x.^-0.5 /180;
ezplot(fun,1,1000); hold on
plot(D(:),abs(kernel(:)),'.');



end
