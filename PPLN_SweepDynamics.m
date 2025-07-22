close all;
clear all;
addpath(fullfile(cd,'utils'))

N.w  = 400*1./4; % Time slices 
N.z  = 40; % LNL/N.z = dz;
lam0.p = 786e-9;
lam0.s = 1572e-9;
I0.p = 4.5; I0.s = 0; I0.i = 0; % CW pump Powers
I0.p_seed = 0; I0.s_seed = 0; I0.i_seed = 0; % Pulse seeding power
I0.t_seed = 1e-12; % FWHM of pulse seed

% Resonator construction
% Resonator loss due to crystal absorptions and/or diffraction
F = 10;
res.alpha_pc = pi/F; 
res.alpha_sc = pi/F; 
res.alpha_ic = pi/F;
% Resonator coupling
res.theta_p  = pi/F;
res.theta_s  = pi/F;
res.theta_i  = pi/F;
% Cavity detuning (normalized to linewidth of pump resonance)
res.delta_p  = 6; 
res.delta_s  = 3;
res.delta_i  = res.delta_p-res.delta_s;

SweepSpeed  = 1/1e13*3e8./10e-3; % in units of linewidth/s
Delta_delta = 12;

% Decide if simulation will consider GVM
res.GVM_on = 0;

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


crystal.dk =0*pi; % normalized to crystal length such that dk = crystal.dk/Lnl
crystal.mode_radius = 3e-6;
crystal.length = 10e-3;
crystal.name = 'LNB_M';
crystal.temperature = 300;
crystal.QPM = 'yes';
crystal.PM_type = 'type-1';
crystal.qpm_currentPolValue = 3;
[crystal,lam0] = Crystal_Calculate(crystal,lam0);

% options = [];
% options.SweepSpeed  =; % In units of linewidth/s... Default Sweep speed = 0 means does not sweep

TROPO_Generalized(lam0,crystal,I0,res,N,'SweepSpeed', SweepSpeed)
tilefigs
