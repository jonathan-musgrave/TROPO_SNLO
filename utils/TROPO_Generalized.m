
function TROPO_Generalized(lam0,crystal,I0,res,N,options)
arguments
   lam0 struct
   crystal struct
   I0 struct
   res struct
   N struct
   
   options.SweepSpeed  = 0; % Default Sweep speed = 0 means does not sweep
   options.SweepResonance ='p'; % Sweep p for pump, s for signal, i for idler
   options.SweepDir = -sign(res.delta_p); % if pos sweep towards negative if negative sweep towards positive
   options.Delta_delta = 2*abs([res.delta_p]); % Default delta is from delta_p0 to -delta_p0 if sweep speed is enabled 
   options.offset_delta = 0;
   options.name =  strcat(crystal.name,'lamp_',num2str(round(lam0.p*1e9)),'lams_',num2str(round(lam0.s*1e9)),'lami_',num2str(round(lam0.i*1e9)),'dk_',num2str(crystal.dk./pi),'pi','Ip0',num2str(I0.p));
   options.save = 0;
   
end
% bright-bright
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%-----general definition-----%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% global c eps_0;  
noise = 1e-12;
c =  2.99792458E8;           % velocity of light in vacuum, m/s
eps_0 = 8.854E-12;           % dielectric constant in vacuum, F/m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% OPO Configuration %%%%%%%%
%%%% Nonlinearity and geometry... %%%%
L = crystal.length; Ac = pi.*(crystal.mode_radius).^2;                   % length of PPLN1, m
deff1 = crystal.deff;                  % nonlinear coefficient for PPLN, m/V
%%%% Coupling and losses
alpha_pc = res.alpha_pc; alpha_sc = res.alpha_sc; alpha_ic = res.alpha_ic; % Total cavity loss from diffraction or other
alphap= alpha_pc/L;   alphas = alpha_sc/L;    alphai = alpha_ic/L; %loss 1/m
% Coupling Transmission set too 0 for non resonant condition
Rp = 1-res.theta_p; % 1-Output coupling
Rs = 1-res.theta_s;
Ri = 1-res.theta_i;

alpha = (alpha_pc+res.theta_p)./2;

%  Initial Detuning and phase-mismatch 
detunep = res.delta_p.*alpha;
detunes = res.delta_s.*alpha;
detunei = res.delta_i.*alpha;

dk  = crystal.dk./crystal.length; % Phase mismatch

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%-----------------Sim parameters---------%%%%%%%%%%%%%%%%%%%%%%%%
Nw = N.w;             % slicing number for time
Nzfi = N.z;             % slicing number for PPLN1
if ~options.SweepSpeed
    Nrt = N.rt;            % Round trips
    detunings = ones(3,Nrt).*[res.delta_p,res.delta_s,res.delta_i]'.*alpha;
else
    Nrt = floor(options.Delta_delta./options.SweepSpeed);
    delta_0 = [res.delta_p,res.delta_s,res.delta_i]*strcmpi({'p','s','i'}',repmat(options.SweepResonance,[3,1]));
    delta_t = delta_0 + options.SweepDir.*(options.Delta_delta);
    detune_ar = linspace(delta_0,delta_t,Nrt).*alpha; 
    
    detunings = repmat(strcmpi({'p','s','i'}',repmat(options.SweepResonance,[3,1])),[1,Nrt]).*detune_ar;
    [~,resPumped] = max(detunings(:,1) == 1);
    % Following from delta_p = delta_s + delta_i + dkL + offset the offset
    % is calculated from the initial delta_p, delta_s, delta_i, and dk
    
    if resPumped == 1
        detunings(2,:) = detunings(1,:).*(res.delta_s+options.offset_delta)./res.delta_p-options.offset_delta;
        detunings(3,:) = detunings(1,:) - detunings(2,:)-dk*L-2*options.offset_delta;
    elseif resPumped == 2
        detunings(1,:) = detunings(2,:).*res.delta_p./res.delta_s;
        detunings(3,:) = detunings(1,:) - detunings(2,:)-dk*L;  
    elseif resPumped == 3
        detunings(1,:) = detunings(3,:).*res.delta_p./res.delta_i;
        detunings(2,:) = detunings(1,:) - detunings(3,:)-dk*L;  
        
    end
    for q = 1:length(resPumped)
        if q == 1 || resPumped(q)
            
        elseif q == 2 || resPumped(q)
            
        elseif q == 3 || resPumped(q)
            
        end
        
    end
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%----------GVM and GDD of PPLN--------%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%---------GV-------%%%%%%%%% %%%%%%%%%---------GVM Compensation -------%%%%%%%%%
beta1_p = crystal.beta1_p;
beta1_s = crystal.beta1_s;
beta1_i = crystal.beta1_i;         % Compensation given in units of s for GVM and s^2 for GVD
beta_offset1 = (beta1_p-beta1_p).*res.GVM_on;      beta1_comp_p = res.GVM_comp_on.*res.GVMp_comp.*(-beta_offset1).*L;
beta_offset2 = (beta1_s-beta1_p).*res.GVM_on;      beta1_comp_s = res.GVM_comp_on.*res.GVMs_comp.*(-beta_offset2).*L;
beta_offset3 = (beta1_i-beta1_p).*res.GVM_on;      beta1_comp_i = res.GVM_comp_on.*res.GVMi_comp.*(-beta_offset3).*L;
%%%%%%%---------GVD-------%%%%%%%%% %%%%%%%---------GVD Compensation-------%%%%%%%%%
beta2_p = crystal.beta2_p;    beta2_comp_p = -res.GDD_comp_on.*res.GDDp_comp.*beta2_p.*L;
beta2_s = crystal.beta2_s;    beta2_comp_s = -res.GDD_comp_on.*res.GDDs_comp.*beta2_s.*L;
beta2_i = crystal.beta2_i;    beta2_comp_i = -res.GDD_comp_on.*res.GDDi_comp.*beta2_i.*L;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


w_p = 2*pi*c/lam0.p; 
w_s = 2*pi*c/lam0.s; 
w_i = 2*pi*c/lam0.i;

n_p = beta1_p*c;
n_s = beta1_s*c;
n_i = beta1_i*c;

kappa_p = sqrt(2)*w_p*deff1/sqrt(n_p*n_s*n_i*c^3*eps_0*Ac); 
kappa_s = sqrt(2)*w_s*deff1/sqrt(n_p*n_s*n_i*c^3*eps_0*Ac);
kappa_i = sqrt(2)*w_i*deff1/sqrt(n_p*n_s*n_i*c^3*eps_0*Ac);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tww= beta1_p*L*2;      % time window, related to the roundtrip time of the pump 
h = [-Nw/2:Nw/2-1];
dt = tww./Nw;     
t = h.*dt;
w = 2*pi/(dt*Nw)*h;
dw = 2*pi/(dt*Nw);                                                                                     
z = linspace(0,L,Nzfi);  % PPLN1                                
dz = mean(diff(z));      % PPLN1                                                           
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%-----------FFT transfer -----------%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

D_p = exp(-alphap/2.*(dz/2)).*exp(1j*(beta_offset1.*w + beta2_p/2.*w.^2)*dz/2); 
D_s = exp(-alphas/2.*(dz/2)).*exp(1j*(beta_offset2.*w + beta2_s/2.*w.^2)*dz/2);
D_i = exp(-alphai/2.*(dz/2)).*exp(1j*(beta_offset3.*w + beta2_i/2.*w.^2)*dz/2); 
% Dispersion Compensation
Dcomp_p = exp(1j*(beta1_comp_p.*w + beta2_comp_p/2.*w.^2)); 
Dcomp_s = exp(1j*(beta1_comp_s.*w + beta2_comp_s/2.*w.^2));
Dcomp_i = exp(1j*(beta1_comp_i.*w + beta2_comp_i/2.*w.^2)); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% CW Pumping variables
Ap0 = sqrt(I0.p)*ones(1,Nw); 
As0 = sqrt(I0.s)*ones(1,Nw); 
Ai0 = sqrt(I0.i)*ones(1,Nw); 

% Pulse Pump Variables
tseed = I0.t_seed;
Ap00 = sqrt(I0.p)*ones(1,Nw) + sqrt(I0.p_seed)*exp(-2*sqrt(2)*(t/tseed).^2);
As00 = sqrt(I0.s)*ones(1,Nw) + sqrt(I0.s_seed)*exp(-2*sqrt(2)*(t/tseed).^2);  
Ai00 = sqrt(I0.i)*ones(1,Nw) + sqrt(I0.i_seed)*exp(-2*sqrt(2)*(t/tseed).^2);  

% Initialize
Ap = noise*exp(2*pi*1j*rand(size(Ap0)))+Ap00;
As = noise*exp(2*pi*1j*rand(size(As0)))+As00;
Ai = noise*exp(2*pi*1j*rand(size(Ai0)))+Ai00;

LP = zeros(Nrt,Nw);
LS = zeros(Nrt,Nw);
LI = zeros(Nrt,Nw);
AZp = zeros(Nzfi,Nw);
AZs = zeros(Nzfi,Nw);
AZi = zeros(Nzfi,Nw);

wbd = 2*pi*9E12;
BD = fftshift(exp(-(w/wbd).^8));
BDt = ones(size(BD));


% fftshift before stepping
D_p = fftshift(D_p);
D_s = fftshift(D_s);
D_i = fftshift(D_i);

Dcomp_p = fftshift(Dcomp_p);
Dcomp_s = fftshift(Dcomp_s);
Dcomp_i = fftshift(Dcomp_i);

for indrt = 1:Nrt
    Nrt - indrt
    detunep = detunings(1,indrt);
    detunes = detunings(2,indrt);
    detunei = detunings(3,indrt);
    if indrt<50
        Ap = (sqrt(1-Rp)).*Ap00 + sqrt(Rp)*Ap*exp(-1j*detunep) + noise*exp(2*pi*1j*rand(size(Ap)));
        As = (sqrt(1-Rs)).*As00 + sqrt(Rs)*As*exp(-1j*detunes) + noise*exp(2*pi*1j*rand(size(As)));
        Ai = (sqrt(1-Ri)).*Ai00 + sqrt(Ri)*Ai*exp(-1j*detunei) + noise*exp(2*pi*1j*rand(size(Ai)));
    else
        Ap = (sqrt(1-Rp)).*Ap0 + sqrt(Rp)*Ap*exp(-1j*detunep) + noise*exp(2*pi*1j*rand(size(Ap)));
        As = (sqrt(1-Rs)).*As0 + sqrt(Rs)*As*exp(-1j*detunes) + noise*exp(2*pi*1j*rand(size(As)));
        Ai = (sqrt(1-Ri)).*Ai0 + sqrt(Ri)*Ai*exp(-1j*detunei) + noise*exp(2*pi*1j*rand(size(Ai)));
    end

    % Propagation (1st half), split-step Fourier method and
    % DispersionCompensation
    sAp = BD.*(ifft(ifftshift(Ap)));
    sAs = BD.*(ifft(ifftshift(As)));
    sAi = BD.*(ifft(ifftshift(Ai)));
    Ap = fftshift(fft((D_p.*sAp.*Dcomp_p))).*BDt;
    As = fftshift(fft((D_s.*sAs.*Dcomp_s))).*BDt;
    Ai = fftshift(fft((D_i.*sAi.*Dcomp_i))).*BDt;
    jj = 0;
    
    for indz = z

        jj = jj + 1;


        % nonlinear step using Runga-Kutta 4th order  
        [Ap, As, Ai] = OPOsol(Ap,As,Ai,kappa_p,kappa_s,kappa_i,dk,indz,dz);
        % Propagation (1st half), split-step Fourier method 
        sAp = BD.*(ifft(ifftshift(Ap)));
        sAs = BD.*(ifft(ifftshift(As)));
        sAi = BD.*(ifft(ifftshift(Ai)));
     % Propagation (2st half), split-step Fourier method
     if indz == z(end)
         Ap = fftshift(fft((D_p.^1.*sAp)));
         As = fftshift(fft((D_s.^1.*sAs)));
         Ai = fftshift(fft((D_i.^1.*sAi)));
     else
         Ap = fftshift(fft((D_p.^2.*sAp)));
         As = fftshift(fft((D_s.^2.*sAs)));
         Ai = fftshift(fft((D_i.^2.*sAi)));
     end

         AZp(jj,:) = Ap;
         AZs(jj,:) = As;
         AZi(jj,:) = Ai;

    end

    
    LP(indrt,:) = Ap;  
    LS(indrt,:) = As; 
    LI(indrt,:) = Ai;          
end
ff = gcf;
fig_num = ff.Number;
if ~options.SweepSpeed
    run('run_plots.m')
else
    run('run_plots_sweep.m')
end
% f1 = figure(1);
% f2 = figure(2);
% f6 = figure(6);
% f6.Position = [2 463 1776 420];
% savefig(f6,fullfile('Sim_Results',strcat(name,'_Time_Evolution.fig')))
% savefig(f1,fullfile('Sim_Results',strcat(name,'_PeakPower_Evolution.fig')))
% savefig(f2,fullfile('Sim_Results',strcat(name,'_OutputLzEvolution.fig')))

end
