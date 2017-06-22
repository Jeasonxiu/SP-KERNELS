function plot_fresnel()

z=0:300;

c=5.0;
T=4.0;

for den = [2.0 3.0 4.0];

    d = c*T / den;

    [FZHW] = Fresnel_Zone_Half_Width(d,z);
    
    plot(FZHW,-z)
    
    hold on

end


end

function [FZHW] = Fresnel_Zone_Half_Width(d,z)

    A=(d+z).^2 - z.^2;
    FZHW=sqrt(A);

end