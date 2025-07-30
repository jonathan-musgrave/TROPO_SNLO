clear all

% bright-bright

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%-----general definition-----%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

c =  2.99792458E8;           % velocity of light in vacuum, m/s
eps_0 = 8.854E-12;           % dielectric constant in vacuum, F/m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%----------GVM and GDD of PPLN--------%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%---------GV-------%%%%%%%%%
beta1_p = 7389e-12 + 0E-12;                              
beta1_s = 7389e-12 + 0E-12;
beta1_i = 7389e-12 + 0E-12;
beta_offset1 = beta1_p-beta1_s;
beta_offset2 = beta1_s-beta1_s;
beta_offset3 = beta1_i-beta1_s;
%%%%%%%---------GVD-------%%%%%%%%%
beta2_p = 281e-27;
beta2_s = 213e-27;
beta2_i = -87e-27*0+1096e-27;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% BGGSE OPO %%%%%%%%
Ac = pi*(3e-6)^2;                % m^2, mode area in PPLN1
L = 3E-3;                       % length of PPLN1, m
deff1 = 27.6E-12;                  % nonlinear coefficient for PPLN, m/V

lamp0 = 2100e-9;                  % wavelength of pump laser for OPO, m
lams0 = 2658.2e-9;                 % wavelength of signal laser for OPO, m
lami0 = 1/(1/lamp0-1/lams0);     % wavelength of idler laser for OPO, m
w_p = 2*pi*c/lamp0; 
w_s = 2*pi*c/lams0; 
w_i = 2*pi*c/lami0;

n_p = beta1_p*c;
n_s = beta1_s*c;
n_i = beta1_i*c;

kappa_p = sqrt(2)*w_p*deff1/sqrt(n_p*n_s*n_i*c^3*eps_0*Ac); 
kappa_s = sqrt(2)*w_s*deff1/sqrt(n_p*n_s*n_i*c^3*eps_0*Ac);
kappa_i = sqrt(2)*w_i*deff1/sqrt(n_p*n_s*n_i*c^3*eps_0*Ac);


%     kappa_p = 156.6925;
%     kappa_s = 65.587;
%     kappa_i = 91.1055;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%-----------------Sim parameters---------%%%%%%%%%%%%%%%%%%%%%%%%
Nw = 1000;             % slicing number for time
Nzfi = 31;             % slicing number for PPLN1

tww= beta1_s*L*2;      % time window, related to the roundtrip time
t = linspace(-tww,tww,Nw);    
dt = mean(diff(t));                                                     
w = 2*pi*linspace(-1/2/dt,1/2/dt,Nw); % frequency 
dw = mean(diff(w));                                                                                     
z = linspace(0,L,Nzfi);  % PPLN1                                
dz = mean(diff(z));      % PPLN1                                                           
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%-----------FFT transfer -----------%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
alphap= pi/160/L; alphas = alphap; alphai = alphap; %loss 1/m
D_p = exp(-alphap/2.*(dz/2)).*exp(1j*(beta_offset1.*w + beta2_p/2.*w.^2)*dz/2); 
D_s = exp(-alphas/2.*(dz/2)).*exp(1j*(beta_offset2.*w + beta2_s/2.*w.^2)*dz/2);
D_i = exp(-alphai/2.*(dz/2)).*exp(1j*(beta_offset3.*w + beta2_i/2.*w.^2)*dz/2); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%---------Initial pulse--------%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ip0 = 18e-3*3.4*1.1; 
tseed = 1.0E-12
Ap0 = sqrt(Ip0)*ones(1,Nw); 
Ap00 = sqrt(Ip0)*ones(1,Nw) + sqrt(Ip0*0)*exp(-2*sqrt(2)*(t/tseed).^2);
As0 = sqrt(Ip0*2)*exp(-2*sqrt(2)*(t/tseed).^2);  
Ai0 = sqrt(Ip0*2)*exp(-2*sqrt(2)*(t/tseed).^2);  

As = As0;
Ap = Ap0;

Rs = 1-pi/160;
Rp = 1-pi/160;


detunes = -alphas*L*6;
detunep = lams0/lamp0*detunes;
   %  detunes =  0.2258;
% detunei = lams0/lami0*detunes;

Nrt = 4*500;                     % roundtrip number
dk  = 1*pi/L;


close all;
clear all;
addpath(fullfile(cd,'utils'))

N.w  = 400; % Time slices 
N.rt = 500*4; % Round trip number
N.z  = 31; % LNL/N.z = dz;
lam0.p = 2100e-9;
lam0.s = 2658.2e-9;
lam0.i = 1/(1/lam0.p-1/lam0.s);
I0.p = 0.1; I0.s = 0; I0.i = 0; % CW pump Powers
I0.p_seed = 0; I0.s_seed = 2*I0.p; I0.i_seed = 0*I0.p; % Pulse seeding power
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
res.theta_i  = 0;

% Decide if simulation will consider GVM
res.delta_s = -6;
res.delta_p = lam0.s/lam0.p*res.delta_s;
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