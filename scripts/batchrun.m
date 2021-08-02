clc;
clear;
close all;

%   t: time points at which a solution is requested
%   l: Thickness (e.g. 483nm)
%   df: Diffusivity (e.g. 2.2e-10 cm^2/s)
%   sc: Surface Concentration (e.g. 3.708e-3 mol/cm^3)
%   pc: Accessible Polymer Concentration (e.g. 5.758e-3 mol/cm^3)
%   hd: Hindering Factor (e.g. 1200 cm^3/mol)
%   k: Reaction Rate (e.g. 1 cm^3/mol s)

t1 = linspace(0,5000,2001);
t2 = linspace(5000,120000,4601);
t = [t1 t2(2:end)];
writematrix(t', '../results/timeIndex.txt');
l = 4.83E-5*2;

D = readtable('../results/maxpro/maxpro_design.csv');
for i = 1:size(D,1)
   disp(i);
   df = D.diffusivity(i);
   sc = D.surface_conc(i);
   pc = D.polymer_conc(i);
   hd = D.hindering(i);
   k = D.reaction_rate(i);
   mass = TMA_PMMA(t, l, df, sc, pc, hd, k);
   writematrix(mass, sprintf('../results/maxpro/run%d.txt', i));
end
