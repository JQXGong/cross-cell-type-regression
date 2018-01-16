% MATLAB Implementation of the Shannon et al. 2004 version (shannon04) model
% for the rabbit ventricular action potential and calcium transient
%
% The shannon04 model is described in the article "A Mathematical Treatment
% of Integrated Ca Dynamics within the Ventricular Myocyte"
% by Thomas R. Shannon, Fei Wang, José Puglisi, Christopher Weber 
% and Donald M. Bers
%
% The article and supplemental materails are available in the
% Biophysical Journal
% Link to Article:
% http://www.cell.com/biophysj/fulltext/S0006-3495(04)73802-3

clear

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Step 1:  Define constants that will not be varied randomly
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% In order to comply with potential paralell computing, reverse the usual
%%% global variable declearation for all 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % Ionic current maximal conductances

Fx_Ks_SL = 0.89;   % dimensionless (in IKs)
Fx_Ks_jct = 0.11;   % dimensionless (in IKs)
pKNa = 0.01833;   % dimensionless (in IKs)

Fx_NaBk_SL = 0.89;   % dimensionless (in INab)
Fx_NaBk_jct = 0.11;   % dimensionless (in INab)
G_NaBk = 0.297e-3;   % milliS_per_microF (in INab)
Fx_Na_SL = 0.89;   % dimensionless (in INa)
Fx_Na_jct = 0.11;   % dimensionless (in INa)
G_INa = 16.0;   % milliS_per_microF (in INa)
G_tof = 0.02;   % milliS_per_microF (in Itof)
G_tos = 0.06;   % milliS_per_microF (in Itos)

Fx_Cl_SL = 0.89;   % dimensionless (in ICl_Ca)
Fx_Cl_jct = 0.11;   % dimensionless (in ICl_Ca)
G_Cl = 0.109625;   % milliS_per_microF (in ICl_Ca)
Kd_ClCa = 0.1;   % millimolar (in ICl_Ca)
G_ClBk = 0.009;   % milliS_per_microF (in IClb)

% % ICaL parameters
Fx_ICaL_SL = 0.1;   % dimensionless (in ICaL)
Fx_ICaL_jct = 0.9;   % dimensionless (in ICaL)
PCa = 5.4e-4;   % litre_per_farad_millisecond (in ICaL)
PK = 2.7e-7;   % litre_per_farad_millisecond (in ICaL)
PNa = 1.5e-8;   % litre_per_farad_millisecond (in ICaL)
Q10_CaL = 1.8;   % dimensionless (in ICaL)
gamma_Cai = 0.341;   % dimensionless (in ICaL)
gamma_Cao = 0.341;   % dimensionless (in ICaL)
gamma_Ki = 0.75;   % dimensionless (in ICaL)
gamma_Ko = 0.75;   % dimensionless (in ICaL)
gamma_Nai = 0.75;   % dimensionless (in ICaL)
gamma_Nao = 0.75;   % dimensionless (in ICaL)

GKr_ = 0.03 ;
GKs_ = 0.07 ;
GK1_ = 0.9 ;

% % SR Ca release parameters

KSRleak = 5.348e-6;   % per_millisecond (in Jleak_SR)
H2 = 1.787;   % dimensionless (H in Jpump_SR)
Kmf = 0.000246;   % millimolar (in Jpump_SR)
Kmr = 1.7;   % millimolar (in Jpump_SR)
Q10_SRCaP = 2.6;   % dimensionless (in Jpump_SR)
V_max_3 = 286.0e-6;   % millimolar_per_millisecond (V_max in Jpump_SR)
EC50_SR = 0.45;   % millimolar (in Jrel_SR)
HSR = 2.5;   % dimensionless (in Jrel_SR)
Max_SR = 15.0;   % dimensionless (in Jrel_SR)
Min_SR = 1.0;   % dimensionless (in Jrel_SR)
kiCa = 0.5;   % per_millimolar_per_millisecond (in Jrel_SR)
kim = 0.005;   % per_millisecond (in Jrel_SR)
koCa = 10.0;   % per_millimolar2_per_millisecond (in Jrel_SR)
kom = 0.06;   % per_millisecond (in Jrel_SR)
ks = 25.0;   % per_millisecond (in Jrel_SR)

% % Na-K pump paramters
Fx_NaK_SL = 0.89;   % dimensionless (in INaK)
Fx_NaK_jct = 0.11;   % dimensionless (in INaK)
H_NaK = 4.0;   % dimensionless (in INaK)
I_NaK_max = 1.91;   % microA_per_microF (in INaK)
Km_Ko = 1.5;   % millimolar (in INaK)
Km_Nai = 11.0;   % millimolar (in INaK)
Q10_Km_Nai = 1.49;   % dimensionless (in INaK)
Q10_NaK = 1.63;   % dimensionless (in INaK)

% % NCX parameters
Fx_NCX_SL = 0.89;   % dimensionless (in INaCa)
Fx_NCX_jct = 0.11;   % dimensionless (in INaCa)
HNa = 3.0;   % dimensionless (in INaCa)
K_mCai = 0.00359;   % millimolar (in INaCa)
K_mCao = 1.3;   % millimolar (in INaCa)
K_mNai = 12.29;   % millimolar (in INaCa)
K_mNao = 87.5;   % millimolar (in INaCa)
Kd_act = 0.000256;   % millimolar (in INaCa)
Q10_NCX = 1.57;   % dimensionless (in INaCa)
V_max_2 = 9.0;   % microA_per_microF (V_max in INaCa)
eta = 0.35;   % dimensionless (in INaCa)
ksat = 0.27;   % dimensionless (in INaCa)


% % SL Ca pump and background currents
Fx_CaBk_SL = 0.89;   % dimensionless (in ICab)
Fx_CaBk_jct = 0.11;   % dimensionless (in ICab)
G_CaBk = 0.0002513;   % milliS_per_microF (in ICab)
Fx_SLCaP_SL = 0.89;   % dimensionless (in ICap)
Fx_SLCaP_jct = 0.11;   % dimensionless (in ICap)
H1 = 1.6;   % dimensionless (H in ICap)
Km = 0.0005;   % millimolar (in ICap)
Q10_SLCaP = 2.35;   % dimensionless (in ICap)
V_maxAF = 0.0673;   % microA_per_microF (in ICap)  %%%%%%%%%%%%%%% in the paper different value and unit,
% probably being converted

% % Buffering parameters
Bmax_Calsequestrin = 0.14;   % millimolar (in Ca_buffer)
Bmax_SLB_SL = 0.0374;   % millimolar (in Ca_buffer)
Bmax_SLB_jct = 0.0046;   % millimolar (in Ca_buffer)
Bmax_SLHigh_SL = 0.0134;   % millimolar (in Ca_buffer)
Bmax_SLHigh_jct = 0.00165;   % millimolar (in Ca_buffer)
koff_Calsequestrin = 65.0;   % per_millisecond (in Ca_buffer)
koff_SLB = 1.3;   % per_millisecond (in Ca_buffer)
koff_SLHigh = 30.0e-3;   % per_millisecond (in Ca_buffer)
kon_Calsequestrin = 100.0;   % per_millimolar_per_millisecond (in Ca_buffer)
kon_SL = 100.0;   % per_millimolar_per_millisecond (in Ca_buffer)

Bmax_SL = 1.65;   % millimolar (in Na_buffer)
Bmax_jct = 7.561;   % millimolar (in Na_buffer)
koff = 1.0e-3;   % per_millisecond (in Na_buffer)
kon = 0.0001;   % per_millimolar_per_millisecond (in Na_buffer)
Bmax_Calmodulin = 0.024;   % millimolar (in cytosolic_Ca_buffer)
Bmax_Myosin_Ca = 0.14;   % millimolar (in cytosolic_Ca_buffer)
Bmax_Myosin_Mg = 0.14;   % millimolar (in cytosolic_Ca_buffer)
Bmax_SRB = 0.0171;   % millimolar (in cytosolic_Ca_buffer)
Bmax_TroponinC = 0.07;   % millimolar (in cytosolic_Ca_buffer)
Bmax_TroponinC_Ca_Mg_Ca = 0.14;   % millimolar (in cytosolic_Ca_buffer)
Bmax_TroponinC_Ca_Mg_Mg = 0.14;   % millimolar (in cytosolic_Ca_buffer)
koff_Calmodulin = 238.0e-3;   % per_millisecond (in cytosolic_Ca_buffer)
koff_Myosin_Ca = 0.46e-3;   % per_millisecond (in cytosolic_Ca_buffer)
koff_Myosin_Mg = 0.057e-3;   % per_millisecond (in cytosolic_Ca_buffer)
koff_SRB = 60.0e-3;   % per_millisecond (in cytosolic_Ca_buffer)
koff_TroponinC = 19.6e-3;   % per_millisecond (in cytosolic_Ca_buffer)
koff_TroponinC_Ca_Mg_Ca = 0.032e-3;   % per_millisecond (in cytosolic_Ca_buffer)
koff_TroponinC_Ca_Mg_Mg = 3.33e-3;   % per_millisecond (in cytosolic_Ca_buffer)
kon_Calmodulin = 34.0;   % per_millimolar_per_millisecond (in cytosolic_Ca_buffer)
kon_Myosin_Ca = 13.8;   % per_millimolar_per_millisecond (in cytosolic_Ca_buffer)
kon_Myosin_Mg = 15.7e-3;   % per_millimolar_per_millisecond (in cytosolic_Ca_buffer)
kon_SRB = 100.0;   % per_millimolar_per_millisecond (in cytosolic_Ca_buffer)
kon_TroponinC = 32.7;   % per_millimolar_per_millisecond (in cytosolic_Ca_buffer)
kon_TroponinC_Ca_Mg_Ca = 2.37;   % per_millimolar_per_millisecond (in cytosolic_Ca_buffer)
kon_TroponinC_Ca_Mg_Mg = 3.0e-3;   % per_millimolar_per_millisecond (in cytosolic_Ca_buffer)

% % Geometric and diffusional parameters
A_SL_cytosol = 1.3e-4;   % cm2 (in ion_diffusion)
A_jct_SL = 3.01e-6;   % cm2 (in ion_diffusion)
D_Ca_SL_cytosol = 1.22e-6;   % dm2_per_second (in ion_diffusion)
D_Ca_jct_SL = 1.64e-6;   % dm2_per_second (in ion_diffusion)
D_Na_SL_cytosol = 1.79e-5;   % dm2_per_second (in ion_diffusion)
D_Na_jct_SL = 1.09e-5;   % dm2_per_second (in ion_diffusion)
x_SL_cytosol = 0.45;   % micrometre (in ion_diffusion)
x_jct_SL = 0.5;   % micrometre (in ion_diffusion)

% % Ionic concentrations 
Cao = 1.8;   % millimolar (in model_parameters)
Cli = 15.0;   % millimolar (in model_parameters)
Clo = 150.0;   % millimolar (in model_parameters)
Ki = 135.0;   % millimolar (in model_parameters)
Ko = 5.4;   % millimolar (in model_parameters)
Mgi = 1.0;   % millimolar (in model_parameters)
Nao = 140.0;   % millimolar (in model_parameters)

% % Physical parameters
F = 96486.7;   % coulomb_per_mole (in model_parameters)
R2 = 8314.3;   % joule_per_kilomole_kelvin (R in model_parameters)
T = 310.0;   % kelvin (in model_parameters)

% % Cell geometry
Cm_per_area = 2.0e-6;   % farad_per_cm2 (in model_parameters)
cell_length = 100.0;   % micrometre (in model_parameters)
cell_radius = 10.25;   % micrometre (in model_parameters)

% ps, default value = 1, scaling the time constant
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% entry created for modifying gating time constants 
%%% baseline values set equal to 1, indicating no modification
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p_d = 1 ;
p_f = 1 ;
p_Xr = 1 ;
p_Xs = 1 ;
p_h = 1 ;
p_j = 1 ;
p_m = 1 ;
p_X_tof = 1 ;
p_Y_tof = 1 ;
p_R_tos = 1 ;
p_X_tos = 1 ;
p_Y_tos = 1 ;

n_ps = 12 ;
baseline_ps = ones(n_ps,1) ;

% Vs, default value = 0, voltage dependense of gating
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% entry created for modifying gating V1/2 shift 
%%% baseline values set equal to 0, indicating no modification
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Vshift_d = 0 ;
Vshift_f = 0 ;
Vshift_K1 = 0 ;
Vshift_Xr = 0 ;
Vshift_Xs = 0 ;
Vshift_h = 0 ;
Vshift_j = 0 ;
Vshift_m = 0 ;
Vshift_X_tof = 0 ;
Vshift_Y_tof = 0 ;
Vshift_R_tos = 0 ;
Vshift_X_tos = 0 ;
Vshift_Y_tos = 0 ;

n_Vs = 13 ;
baseline_Vs = zeros(n_Vs,1) ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Step 2:  Define simulation, stimulus, and recording parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

stiminterval = 1000 ; % interval bewteen stimuli, ms
n_stimuli = 120 ; % number of stimuli
tend =  n_stimuli*stiminterval ; % end of simulation, ms
stimdelay = 20 ; % delay between start of each cycle and beginning of current injection, ms
stimdur = 2 ; % length of current injection, ms
stim_amp = 20 ; % amplitude of current injection, nA

stim_starts = stimdelay + stiminterval*(0:n_stimuli-1)  ;
stim_ends = stim_starts + stimdur ;

simints = 3*n_stimuli ;
for i=1:n_stimuli
  intervals(3*i-2,:) = [stiminterval*(i-1),stim_starts(i)] ;
  intervals(3*i-1,:) = [stim_starts(i),stim_ends(i)] ;
  intervals(3*i,:) = [stim_ends(i),stiminterval*i] ;
end
intervals(end,:) = [stim_ends(end),tend] ;

Istim = zeros(simints,1) ;
stimindices = 3*(1:n_stimuli) - 1 ;
Istim(stimindices) = -stim_amp ;


numbertokeep = 1 ; % number of stimulus to record state variable values for

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Step 3:  Set initial conditions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Original valus as specified in the Shannon et al. article
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

% % Ca_Calsequestrin = 1.1865 ;
% % Ca_SL = 0.0001064 ;
% % Ca_SLB_SL = 0.0098686 ;
% % Ca_SLB_jct = 0.0077808 ;
% % Ca_SLHigh_SL = 0.11444 ;
% % Ca_SLHigh_jct = 0.077504 ;
% % Ca_SR = 0.54561 ;
% % Ca_jct = 0.00017484 ;
% % Cai = 8.735e-005 ;
% % d = 6.9975e-006 ;
% % fCaB_SL = 0.015353 ;
% % fCaB_jct = 0.024609 ;
% % f = 1.0007 ;
% % Xr = 0.0084716 ;
% % Xs = 0.006874 ;
% % h = 0.98714 ;
% % j = 0.99182 ;
% % m = 0.0013707 ;
% % X_tof = 0.0040112 ;
% % Y_tof = 0.99463 ;
% % R_tos = 0.38343 ;
% % X_tos = 0.0040113 ;
% % Y_tos = 0.29352 ;
% % I = 9.272e-008 ;
% % O = 7.1126e-007 ;
% % R1 = 0.88467 ;
% % Na_SL = 8.8741 ;
% % Na_SL_buf = 0.77612 ;
% % Na_jct = 8.8728 ;
% % Na_jct_buf = 3.5571 ;
% % Nai = 8.8745 ;
% % V = -85.7197 ;
% % Ca_Calmodulin = 0.00029596 ;
% % Ca_Myosin = 0.0019847 ;
% % Ca_SRB = 0.0021771 ;
% % Ca_TroponinC = 0.0089637 ;
% % Ca_TroponinC_Ca_Mg = 0.118 ;
% % Mg_Myosin = 0.1375 ;
% % Mg_TroponinC_Ca_Mg = 0.010338 ;

load InitialCondition_4models_1Hz.mat
statevar_start = InitialCondition.shannon04 ;
%%% To improve simulation efficiency and ensure model steady state arrival,
%%% proper initial conditions are pre-identified as the last time step 
%%% state variable values following a 1000-second simulation started with
%%% originally specified values, under the same simulation protocol

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Step 4:  ODE Solver
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            
statevar_i = statevar_start ;            

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ORI_globals = [ Bmax_Calsequestrin Bmax_SLB_SL Bmax_SLB_jct Bmax_SLHigh_SL Bmax_SLHigh_jct...
    koff_Calsequestrin koff_SLB koff_SLHigh kon_Calsequestrin kon_SL Fx_ICaL_SL...
    Fx_ICaL_jct PCa PK PNa Q10_CaL gamma_Cai gamma_Cao gamma_Ki...
    gamma_Ko gamma_Nai gamma_Nao Fx_CaBk_SL Fx_CaBk_jct G_CaBk Fx_SLCaP_SL...
    Fx_SLCaP_jct H1 Km Q10_SLCaP V_maxAF Fx_Cl_SL Fx_Cl_jct...
    G_Cl Kd_ClCa G_ClBk Fx_Ks_SL Fx_Ks_jct pKNa Fx_NCX_SL Fx_NCX_jct...
    HNa K_mCai K_mCao K_mNai K_mNao Kd_act Q10_NCX V_max_2 eta ksat Fx_NaK_SL...
    Fx_NaK_jct H_NaK I_NaK_max Km_Ko Km_Nai Q10_Km_Nai Q10_NaK Fx_NaBk_SL...
    Fx_NaBk_jct G_NaBk Fx_Na_SL Fx_Na_jct G_INa G_tof G_tos KSRleak...
    H2 Kmf Kmr Q10_SRCaP V_max_3 EC50_SR HSR Max_SR Min_SR kiCa kim koCa...
    kom ks Bmax_SL Bmax_jct koff kon Bmax_Calmodulin Bmax_Myosin_Ca...
    Bmax_Myosin_Mg Bmax_SRB Bmax_TroponinC Bmax_TroponinC_Ca_Mg_Ca...
    Bmax_TroponinC_Ca_Mg_Mg koff_Calmodulin koff_Myosin_Ca koff_Myosin_Mg...
    koff_SRB koff_TroponinC koff_TroponinC_Ca_Mg_Ca koff_TroponinC_Ca_Mg_Mg...
    kon_Calmodulin kon_Myosin_Ca kon_Myosin_Mg kon_SRB kon_TroponinC...
    kon_TroponinC_Ca_Mg_Ca kon_TroponinC_Ca_Mg_Mg...
    Cao Cli Clo F Cm_per_area Ki Ko Mgi Nao R2 T cell_length cell_radius...
    GKr_ GKs_ GK1_...
    p_d p_f p_Xr p_Xs p_h p_j p_m p_X_tof p_Y_tof p_R_tos p_X_tos p_Y_tos...
    Vshift_d Vshift_f Vshift_K1 Vshift_Xr Vshift_Xs Vshift_h...
    Vshift_j Vshift_m Vshift_X_tof Vshift_Y_tof Vshift_R_tos...
    Vshift_X_tos, Vshift_Y_tos ] ;

options = [] ;  
if (n_stimuli > numbertokeep)
    for i=1:simints-3*numbertokeep
        [post,posstatevars] = ode15s(@dydt_shannon04,intervals(i,:),statevar_i,options,Istim(i),ORI_globals) ;
        statevar_i = posstatevars(end,:) ;
        t = post(end) ;
    end
    statevars = statevar_i;
    for i=simints-3*numbertokeep+1:simints
        [post,posstatevars] = ode15s(@dydt_shannon04,intervals(i,:),statevar_i,options,Istim(i),ORI_globals) ;
        t = [t;post(2:end)] ;
        statevars = [statevars;posstatevars(2:end,:)] ;
        statevar_i = posstatevars(end,:) ;
    end
        
else % n_stimuli = numbertokeep
    t = 0 ;
    statevars = statevar_i;
    for i=1:simints
        [post,posstatevars] = ode15s(@dydt_shannon04,intervals(i,:),statevar_i,options,Istim(i),ORI_globals) ;
        t = [t;post(2:end)] ;
        statevars = [statevars;posstatevars(2:end,:)] ;
        statevar_i = posstatevars(end,:) ;
    end
end

t = t - min(t) ;

outputcell = num2cell(statevars,1) ;

[Ca_Calsequestrin,Ca_SL,Ca_SLB_SL,Ca_SLB_jct,Ca_SLHigh_SL,Ca_SLHigh_jct, ...
  Ca_SR,Ca_jct,Cai,d,fCaB_SL,fCaB_jct,f,Xr,Xs,h,j,m,X_tof,Y_tof,R_tos,X_tos, ...
  Y_tos,I,O,R1,Na_SL,Na_SL_buf,Na_jct,Na_jct_buf,Nai,V,Ca_Calmodulin, ...
  Ca_Myosin,Ca_SRB,Ca_TroponinC,Ca_TroponinC_Ca_Mg,Mg_Myosin, ...
  Mg_TroponinC_Ca_Mg] ... 
  = deal(outputcell{:}) ;

figure
subplot(1,2,1)
hold on
title('Action Potential')
plot(t,V)
xlabel('time (ms)')
ylabel('V_m (mV)')
hold off

subplot(1,2,2)
hold on
title('Calcium Transient')
plot(t,Cai)
xlabel('time (ms)')
ylabel('[Cai]_i (mM)')
hold off






