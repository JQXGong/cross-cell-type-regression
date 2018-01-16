% MATLAB Implementation of the Luo-Rudy 2009 version(lr09) model for 
% guinea pig ventricular cell action potential and calcium transient
%
% The lr09 model is described in the article "Uniqueness and Stability 
% of Action Potential Models during Rest, Pacing, and Conduction Using
% Problem-Solving Environment"
% by Leonid Livshitz and Yoram Rudy

% The article and supplemental materails are available in the
% Biophysical Journal
% Link to Article:
% http://www.cell.com/biophysj/fulltext/S0006-3495(09)01159-X
%    
%    t   time variable
%    V   membrane potantial

clear

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Step 1:  Define constants that will not be varied randomly
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% In order to comply with potential paralell computing, reverse the usual
%%% global variable declearation for all 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% global F R T RTF FRT ;
% Physical constants

F = 96485 ;                  % Faraday's constant, C/mol
R = 8314 ;                   % gas constant, mJ/K
T = 273+37 ;                 % absolute temperature, K 
RTF = R*T/F ;
FRT = 1/RTF ;

%%% global Acap Vmyo Vss VJSR VNSR ;
length_cell = 0.01;       % Length of the cell (cm)
radius = 0.0011;     % Radius of the cell (cm)
Vcell = 1000*pi*radius*radius*length_cell ;     %   3.801e-5 uL   % Cell volume (uL)
Ageo = 2*pi*radius*radius+2*pi*radius*length_cell ;  %   7.671e-5 cm^2    % Geometric membrane area (cm^2)
Acap = 2*Ageo ;             %   1.534e-4 cm^2    % Capacitive membrane area (cm^2)
Vmyo = Vcell*0.68 ;    % Myoplasm volume (uL)
Vmito = Vcell*0.24 ;  % Mitochondria volume (uL)
%% data.vsr = Vcell*0.06;    % SR volume (uL)
VNSR = Vcell*0.0552 ;   % NSR volume (uL)
VJSR =  Vcell*0.0048 ;   % JSR volume (uL)
Vss = Vcell*0.02 ;

% % These numbers from Hund-Rudy model, should compare
% % % Cell geometry constants
% % Acap = 1.534e-4 ;            % cm^2
% % Vmyo = 25.84e-6 ;            % uL
% % Vss = 0.76e-6 ;              % uL
% % VJSR = 0.182e-6 ;            % uL
% % VNSR = 2.098e-6 ;            % uL

% Fixed ionic concentrations
% Initial conditions of others listed below
%%% global Ko Nao Cao ;
Ko = 4.5 ;                  % uM
Nao = 140 ;               % uM
Cao = 1.8 ;                 % uM

%%% global GNa_ GNab 
% % Na current
GNa_= 16 ;                          % mS/cm^2
GNab = 0.004 ;

% % Ca current
%%% global GCaT GCab PCa PCa_Na PCa_K gamma_Cao gamma_Cai 
%%% global gamma_Nao gamma_Nai gamma_Ko gamma_Ki KmCa
PCa = 5.4e-4 ;                        % cm/s
PCa_Na = 6.75e-7 ;                    % cm/s
PCa_K = 1.93e-7 ;                     % cm/s
PCab = 1.995084e-7 ;                  % cm/s
gamma_Cao = 0.341 ;                   % dimensionless
gamma_Cai = 1 ;                       % dimensionless
gamma_Nao = 0.75 ;                    % dimensionless
gamma_Nai = 0.75 ;                    % dimensionless
gamma_Ko = 0.75 ;                     % dimensionless
gamma_Ki = 0.75 ;                     % dimensionless
% hLca = 1 ;                            % dimensionless, Hill coefficient
% % Need to make sure this variable not used elsewhere
KmCa = 6e-4 ;                         % Half saturation constant, mM

% % T-type & background currents
GCaT = 0.05 ;
GCab = 0.003016 ;

% % % % % % Not quite sure where these fit
% % % data.IKsCa_max=0.6;
% % % data.IKsCa_Kd=38e-6;

%%% global GK1_ GKr_ GKs_ pKNa GKp_ 
GK1_ = 0.75;
GKr_ = 0.02614 ;
GKs_ = 0.433 ;
pKNa = 0.01833 ;                  % relative permeability of IKs, Na to K
GKp_ = 5.52e-3 ;

%%% global INaK_ KmNa_NaK KmK_NaK
INaK_ = 2.25 ;             % Max. current through Na-K pump (uA/uF)
KmNa_NaK = 10 ;             % Half-saturation concentration of NaK pump (mM)
KmK_NaK = 1.5 ;             % Half-saturation concentration of NaK pump (mM)

%%% global kNCX ksat eta
kNCX = 0.00025 ;
ksat = 0.0001 ;
eta = 0.15 ;

%%% global alpha_rel Krel_inf hrel beta_tau Krel_tau
alpha_rel = 0.125 ;
Krel_inf = 1 ;
hrel = 9 ;
beta_tau = 4.75 ;
Krel_tau = 0.0123 ;

%%% global IpCa_ KmpCa
IpCa_ = 1.15 ; % Max. Ca current through sarcolemmal Ca pump (uA/uF)
KmpCa = 5e-4 ; % Half-saturation concentration of sarcolemmal Ca pump (mM)

%%% global Vserca Kmserca CaNSR_max tau_transfer
Vserca = 8.75e-3 ;              % mM/ms 
Kmserca = 9.0e-4 ;              % mM
CaNSR_max = 15.0 ;
tau_transfer = 120 ;

%%% global TRPNtot KmTRPN CMDNtot KmCMDN CSQNtot KmCSQN 
% % Buffers in cytosol
TRPNtot = 70e-3 ;    
KmTRPN = 0.5e-3 ;    

% % % % troponin buffering might be more complicated now, need to check this
% trpnf = 40;   % forward  buffered in TRPN (mM)
% trpnb = 0.02;   % backward  TRPN (mM)
% cmdnf = 100;   % forward  buffered in TRPN (mM)
% cmdnb = 0.238;   % backward  TRPN (mM)

CMDNtot = 50e-3 ;
KmCMDN = 2.38e-3 ;    

% % Buffers in JSR
CSQNtot = 10 ;
KmCSQN = 0.8 ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% entry created for modifying gating time constants 
%%% baseline values set equal to 1, indicating no modification
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% global p_m p_h p_j p_d p_f p_b p_g p_xKr p_xs p_release 
p_m = 1 ;
p_h = 1 ;
p_j = 1 ;
p_d = 1 ;
p_f = 1 ;
p_b = 1 ;
p_g = 1 ;
p_xKr = 1 ;
p_xs = 1 ;
p_release = 1 ;

n_ps = 10 ;
baseline_ps = ones(n_ps,1) ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% entry created for modifying gating V1/2 shift 
%%% baseline values set equal to 0, indicating no modification
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% global Vshift_m Vshift_h Vshift_j Vshift_d Vshift_f 
%%% global Vshift_b Vshift_g Vshift_xKr Vshift_xs  
%%% global Vshift_K1 Vshift_RKr Vshift_Kp Vshift_NaK Vshift_NCX 
Vshift_m = 0 ;
Vshift_h = 0 ;
Vshift_j = 0 ;
Vshift_d = 0 ;
Vshift_f = 0 ;
Vshift_b = 0 ;
Vshift_g = 0 ;
Vshift_xKr = 0 ;
Vshift_xs = 0 ;
Vshift_K1 = 0 ;
Vshift_RKr = 0 ;
Vshift_Kp = 0 ;
Vshift_NaK = 0 ;
Vshift_NCX = 0 ;

n_Vs = 14 ;
baseline_Vs = zeros(n_Vs,1) ;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Step 2:  Define simulation, stimulus, and recording parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

stiminterval = 1000 ; % interval bewteen stimuli, ms
n_stimuli = 120; % number of stimuli
tend =  n_stimuli*stiminterval ; % end of simulation, ms
stimdelay = 20 ; % delay between start of each cycle and beginning of current injection, ms
stimdur = 1 ; % length of current injection, ms
stim_amp = 40 ; % amplitude of current injection, nA


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
%%% Original valus as specified in the Livshitz and Rudy article
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% %     V = -84.7 ;
% %     Cai = 0.0822e-3 ;
% %     CaNSR = 1.25 ;
% %     CaJSR = 1.25 ;
% %     Nai = 9.71 ;
% %     Ki = 142.82 ;
% %     m = 2.46e-4 ;
% %     h = 0.99869 ;
% %     j = 0.99887 ;
% %     d = 1e-4 ;
% %     f = 0.983 ;
% %     b = 1e-4 ;
% %     g = 0.983 ;
% %     xKr = 0.229 ;
% %     xs1 = 1e-4 ;
% %     xs2 = 1e-4 ;
% %     Jrel = 1e-4 ;

load InitialCondition_4models_1Hz.mat
statevar_start = InitialCondition.lr09 ;
%%% To improve simulation efficiency and ensure model steady state arrival,
%%% proper initial conditions are pre-identified as the last time step 
%%% state variable values following a 1000-second simulation started with
%%% originally specified values, under the same simulation protocol

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Step 4:  ODE Solver
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            
statevar_i = statevar_start ;
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 ORI_globals = [F R T RTF FRT Acap Vmyo Vss VJSR VNSR Ko Nao Cao...
     GNa_ GNab GCaT GCab PCa PCa_Na PCa_K gamma_Cao gamma_Cai...
     gamma_Nao gamma_Nai gamma_Ko gamma_Ki KmCa GK1_ GKr_ GKs_ pKNa GKp_...
     INaK_ KmNa_NaK KmK_NaK kNCX ksat eta alpha_rel Krel_inf hrel beta_tau Krel_tau...
     IpCa_ KmpCa Vserca Kmserca CaNSR_max tau_transfer TRPNtot KmTRPN CMDNtot KmCMDN CSQNtot KmCSQN...
     p_m p_h p_j p_d p_f p_b p_g p_xKr p_xs p_release...
     Vshift_m Vshift_h Vshift_j Vshift_d Vshift_f Vshift_b Vshift_g Vshift_xKr Vshift_xs...
     Vshift_K1 Vshift_RKr Vshift_Kp Vshift_NaK Vshift_NCX] ;

 
options = [] ;
if (n_stimuli > numbertokeep)
    for i=1:simints-3*numbertokeep
        [post,posstatevars] = ode15s(@dydt_lr09,intervals(i,:),statevar_i,options,Istim(i),ORI_globals) ;
        statevar_i = posstatevars(end,:) ;
        t = post(end) ;
    end
    statevars = statevar_i;
    for i=simints-3*numbertokeep+1:simints
        [post,posstatevars] = ode15s(@dydt_lr09,intervals(i,:),statevar_i,options,Istim(i),ORI_globals) ;
        t = [t;post(2:end)] ;
        statevars = [statevars;posstatevars(2:end,:)] ;
        statevar_i = posstatevars(end,:) ;
    end
        
else % n_stimuli = numbertokeep
    t = 0 ;
    statevars = statevar_i;
    for i=1:simints
        [post,posstatevars] = ode15s(@dydt_lr09,intervals(i,:),statevar_i,options,Istim(i),ORI_globals) ;
        t = [t;post(2:end)] ;
        statevars = [statevars;posstatevars(2:end,:)] ;
        statevar_i = posstatevars(end,:) ;
    end
end   

  t = t - min(t) ;

  outputcell = num2cell(statevars,1) ;

  [V,Cai,CaNSR,CaJSR,Nai,Ki,m,h,j,d,f,b,g,xKr,xs1,xs2,Jrel] = ...
    deal(outputcell{:}) ;

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
