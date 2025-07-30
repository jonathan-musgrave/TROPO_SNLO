close all;
clear all;
addpath(fullfile(cd,'utils'))

N.w  = 4000; % Time slices 
N.rt = 500*4; % Round trip number
N.z  = 31; % LNL/N.z = dz;
lam0.p = 2100e-9;
lam0.s = 2658.2e-9;
lam0.i = 1/(1/lam0.p-1/lam0.s);

I0.p = 1.885; I0.s = 0; I0.i = 0; % CW pump Powers
I0.p_seed = I0.p; I0.s_seed = I0.p; I0.i_seed = 0*I0.p; % Pulse seeding power
I0.t_seed = 1e-12; % FWHM of pulse seed

% Resonator construction
% Resonator loss due to crystal absorptions and/or diffraction
F = 160;
res.alpha_pc = pi/F; 
res.alpha_sc = pi/F; 
res.alpha_ic = pi/F;
% Resonator coupling
res.theta_p  = pi/F;
res.theta_s  = pi/F;
res.theta_i  = 1;

% Decide if simulation will consider GVM
res.delta_s = -6;
res.delta_p = lam0.s./lam0.p*res.delta_s;
res.GVM_on = 0;
res.delta_i = 0;

% Decide if simulation will consider exteneral GVM and GDD compensation
res.GVM_comp_on = 0;
res.GDD_comp_on = 0;

% If so, decide on the total % compensation that will be placed externally
% For example if GDDp = 1ps^2/mm and Lnl = 1mm than a choise of
% res.GVMp_comp = 1 will apply a compensation of -1ps^2 at the B.C. of the
% cavity
res.GVMp_comp = 1;
res.GVMs_comp = 1;
res.GVMi_comp = 1;

res.GDDp_comp = 1;
res.GDDs_comp = 1;
res.GDDi_comp = 1;


crystal.dk = 1*pi; % normalized to crystal length such that dk = crystal.dk/Lnl
crystal.mode_radius = 3e-6;
crystal.length = 3e-3;
crystal.name = 'BGGSE';
crystal.temperature = 300;
crystal.QPM = 'no';
crystal.PM_type = 'type-1';
crystal.deff = 27.6e-12;
crystal.qpm_currentPolValue = 1;
ng_pump = 2.2152;
ng_signal = 2.2152;
ng_idler = 2.2152;

GDD_pump = 281;
GDD_signal = 213;
GDD_idler = -1096;
c =  2.99792458E8;           % velocity of light in vacuum, m/s
eps_0 = 8.854E-12;           % dielectric constant in vacuum, F/m
crystal.beta1_p = (c./ng_pump).^-1; crystal.beta1_s = (c./ng_signal).^-1; crystal.beta1_i = (c./ng_idler).^-1;
crystal.beta2_p = GDD_pump.*1e-15*1e-15./1e-3; crystal.beta2_s = GDD_signal.*1e-15*1e-15./1e-3; crystal.beta2_i = GDD_idler.*1e-15*1e-15./1e-3;
    

%% Run different cases PLotting
% Detuning negative
% Cavity detuning (normalized to linewidth of pump resonance)
ff = figure(2);
fig_num = ff.Number;
TROPO_Generalized(lam0,crystal,I0,res,N)
ff = gcf;
run('tilefigs.m')