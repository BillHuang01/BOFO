clc;
clear;
close all;

%   t: time points at which a solution is requested
%   l: Thickness (e.g. 480nm)
%   df: Diffusivity (e.g. 2.2e-10 cm^2/s)
%   sc: Surface Concentration (e.g. 3.708e-3 mol/cm^3)
%   pc: Accessible Polymer Concentration (e.g. 5.758e-3 mol/cm^3)
%   hd: Hindering Factor (e.g. 1200 cm^3/mol)
%   k: Reaction Rate (e.g. 1 cm^3/mol s)

t1 = linspace(0,5000,2001);
t2 = linspace(5000,120000,4601);
t = [t1 t2(2:end)];
l = 4.8E-5*2;

% load other variables
fileinfo = dir("*.txt");
filename = fileinfo(1).name;
var = readmatrix(filename);
noRun = var(1);
df = var(2);
sc = var(3);
pc = var(4);
hd = var(5);
k = var(6);
mass = TMA_PMMA(t, l, df, sc, pc, hd, k);

runinfo = split(filename,'_');
runinfo = string(runinfo(1));
writematrix(mass, sprintf('../results/%s/run%d.txt', runinfo, noRun));