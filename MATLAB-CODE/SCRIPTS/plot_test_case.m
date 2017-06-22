xs=-450:5:150; zs=1:5:300;
[XS,ZS]=meshgrid(xs,zs);
XS=reshape(XS,1,[]);
ZS=reshape(ZS,1,[]);
scatter(XS,ZS,5,ZS>100,'filled')
axis equal
set(gca,'Ydir','reverse')
ylabel('Depth (km)')
xlabel('Lateral Offset (km)')

set(gcf, 'PaperPositionMode', 'auto');
print -depsc2 test_case.eps
close;

!gs test_case.eps

