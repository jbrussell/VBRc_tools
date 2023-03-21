function [ Te,presG,rho ] = calc_HSC( Tp,age,spr_rate_cmyr,depth_m )
% HSC from Turcotte and Schubert as implimented by Faul's scripts.
% Assumes spreading rate of 5 cm/year by default
%
% INPUT
% Tp  : Mantle potential temperature in Kelvin 
% age : Seafloor age in Myr
% spr_rate_cmyr: Spreading rate in cm/yr
% depth_m : depth below seafloor in m
%
% OUTPUT
% Te    : temperature in Kelvin
% presG : pressure in GPa
% rho   : density in kg/m3
% 
% JBR 10/24/19

% age = 70; % Ma
% spr_rate_cmyr = 5; % cm/yr
dist = age *spr_rate_cmyr*1e4;
% Tp = 1623; % mantle potential temperature for geotherm, K

%depth sampling, m
% depth_m = [(5000:2000:197000),(200000:5000:400000)]; 
% dk = -depth_m(2:length(depth_m))/1000; lz = length(depth_m);

% Tr = 1173; Pr = 0.2; %reference temperature and pressure, part of fit
% parameters for T and density calculation (Turcotte and Schubert)
grav = 9.98; kappa = 1E-6; 
rho0 = 3310; 
betaa = 6E-12;

% adiabatic T and rho increase, Turcotte and Schubert, Ch. 4-16 and 4-27
Tad = Tead(depth_m,Tp);
% add conductive cooling
vel = spr_rate_cmyr/100/(60*60*24*365); age = (dist/vel)/(60*60*24*365*1E6);
Te = 273 + (Tad - 273) .* erf(depth_m/(2*sqrt(kappa*dist/vel)));
ff = find(Te+5 > Tad);

% calculate density and pressure
rr = 1 - rho0*grav*betaa*depth_m;
pres = -log(rr) / betaa;
presG = pres/1E9;
rho = rho0./rr;

end

function [Tad] = Tead(depth_m,Tp)
% Calculates the adiabatic temperature gradient for a given Tp
cp = 1350; alv = 2.9E-5; grav = 9.98;
Tad = alv * grav * Tp * depth_m/cp + Tp;
end

