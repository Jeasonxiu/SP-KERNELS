function [XSP] = scattering_pattern(THE,alpha,beta,rho,dbeta_over_beta)
%See Sato et al. (2012), Eqn. 4.62 or Rondenay (2009) or Bostock and
%Rondenay (1999)

%gam0=sqrt(3);
%a=5;
%T=2.0;
%el0=2*pi/(beta*T);

%fprintf(' CALCULATING SCATTERING PATTERN\n\n')
%fprintf(' dbeta/beta, beta, alpha, rho = %.3f  %.3f  %.3f  %.3f\n\n',dbeta_over_beta, beta, alpha, rho)

fac = -rho * (dbeta_over_beta * 2 * beta/alpha);

%fprintf(' fac = %f\n\n', fac)

%ARG = el0/gam0 * sqrt(1+gam0^2-2*gam0*cos(THE));

%Pterm = dbeta_over_beta * sqrt(exp(-ARG.^2*a^2 / 4));

XSP = fac * sind(2*THE); %.* Pterm;

%Inspired by eqn 7 of Rondenay, Bostock, and Fischer
%1/(2pi) rather than 1/(4pi) because this is 2-D migration rather than 3d

XSP = XSP * 1/2/pi;

end
