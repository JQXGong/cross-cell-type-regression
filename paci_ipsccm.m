% MATLAB Implementation of the Paci model for human induced pluripotent stem cell
% derived cardiomyocytes action potential and calcium transient
%
% The Paci model is described in the article "Computational Models of 
% Ventricular- and Atrial-Like Human Induced Pluripotent Stem Cell Derived 
% Cardiomyocytes"
% by Michelangelo Paci, Jari Hyttinen, Katriina Aalto-Setälä, and Stefano Severi
%
% The article and supplemental materails are available in the
% journal Annals of Biomedical Engineering
% Link to article:
% https://link.springer.com/article/10.1007%2Fs10439-013-0833-3


clear
celltype = 1; % 1 for ventricular 2 for atrial

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Step 1:  Define constants that will not be varied randomly
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% In order to comply with potential paralell computing, reverse the usual
%%% global variable declearation for all 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% extracellular and intracellular ionic concentrations
%%% global Nao Ko Cao celltype
Nao = 151; % mM
Ko = 5.4; %mM   
Cao = 1.8; %mM

%%% global Ki 
Ki = 150; %mM

% cell size and dimensions
%%% global Cm Vc Vsr

if celltype==1
        Cm = 98.7109e-12 ; % F
        Vc = 8800 ; % um^3
        Vsr = 583.73 ; % um^3
else
        %Cm = 78.6671 ; % pF
        Cm = 78.6671e-12 ; % F
        Vc = 7012 ; % um^3
        Vsr = 465.20 ; % um^3
end

% computed quantities that do not change during simulation
% maximum conductances and currents
% As always, adjusted for capacitance

%%% global GCaL GKr GKs Gf arel brel crel Vleak GpCa GbNa GbCa G_RyR
GCaL = 8.635702e-5 ; % m^3 /(Fxs)
GKr = 29.8667 ; % S/F
GKs = 2.041 ; % S/F
Gf = 30.10312 ; % S/F
arel = 16.464 ; % mM/s
brel = 0.25 ; % mM
crel = 8.232 ; % mM/s
Vleak = 4.4444e-4 ; % (1/s)
GpCa = 0.4125 ; % A/F
GbNa = 0.9 ; % S/F
GbCa = 0.69264 ; % S/F
G_RyR = 1 ;

%%% global GNa Gto GK1 PNaK KNaCA Vmaxup
if celltype==1
        GNa = 3.6712302e3 ; % S/F
        Gto = 29.9038 ; % S/F
        GK1 = 28.1492 ; % S/F
        PNaK = 1.841424 ; % A/F
        KNaCA = 4900 ; % A/F
        Vmaxup = 0.56064 ; % mM/s
else
        GNa = 6.646185e3 ; % S/F
        Gto = 59.8077 ; % S/F
        GK1 = 19.1925 ; % S/F
        PNaK = 1.4731392 ; % A/F
        KNaCA = 2450 ; % A/F
        Vmaxup = 0.22 ; % mM/s     
end

% other constants

%%% global Bufc Bufsr Kbufc Kbufsr Kup KpCa F R T L0 Pkna Ksat KmCa KmNai alpha gamma KmNa Kmk 
Bufc = 0.25 ; %mM
Bufsr = 10 ; %mM
Kbufc = 0.001 ; %mM
Kbufsr = 0.3 ; %mM
Kup = 0.00025 ; % mM
KpCa = 0.0005 ; %mM
F = 96485.3415 ; % (C/M)
R = 8.314472 ; %(J/(M x K))
T = 310 ; % K
L0 = 0.025 ; % dimensionless
Pkna = 0.03 ; % dimensionless
Ksat = 0.1 ; % dimensionless
KmCa = 1.38 ; % mM
KmNai = 87.5 ; %mM
alpha = 2.8571432 ; % dimensionless
gamma = 0.35 ; % dimensionless
KmNa = 40 ; % mM
Kmk = 1 ; % mM

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Step 2:  Define simulation, stimulus, and recording parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

STIM = 1 ; % 0 - spontaneous, 1 - electrical stimulating

%%% stimuli conditions 
stiminterval = 1 ; % interval bewteen stimuli, s 
n_stimuli = 120 ; % number of stimuli
tend =  n_stimuli*stiminterval ; % end of simulation, s
stimdelay = 0.1 ;% delay between start of each cycle and beginning of current injection, s
stimdur = 0.001 ; % length of current injection, s
stim_amp = 40 ; % amplitude of current injection, nA

stim_starts = stimdelay + stiminterval*(0:n_stimuli-1)  ;
stim_ends = stim_starts + stimdur ;

simints = 3*n_stimuli ; %for every stimulus there are 3 states, delay, during stimulus and after stimulus
for i=1:n_stimuli
  intervals(3*i-2,:) = [stiminterval*(i-1),stim_starts(i)] ; % time during delay
  intervals(3*i-1,:) = [stim_starts(i),stim_ends(i)] ;       % time during stimulus   
  intervals(3*i,:) = [stim_ends(i),stiminterval*i] ;         % after stimulus   
end
intervals(end,:) = [stim_ends(end),tend] ;
% when stimuli is applied, intervals define the relevant timestep

Istim = zeros(simints,1) ;
stimindices = 3*(1:n_stimuli) - 1 ; % where there needs a stimulus
Istim(stimindices) = -stim_amp ;


if STIM == 1 % apply electrical stimuli
    stimulitokeep = 1 ; % number of stimulus to record state variable values for
    numbertokeep = stimulitokeep ;
else % spontaneous 
    tend = 120 ; % end of simulation, s
    secondtokeep = 30 ;% length of time to record state variable values for
end
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Step 3:  Set initial conditions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if STIM == 0 % spontaneous
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Original valus as specified in the Paci et al. article
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    h0 = 0.75 ; % dimensionless
    j0 = 0.75 ; % dimensionless
    m0 = 0 ; % dimensionless
    d0 = 0 ; % dimensionless
    fCa0 = 1 ; % dimensionless
    f1comma0 = 1 ; % dimensionless
    f2comma0 = 1 ; % dimensionless
    r0 = 0 ; % dimensionless
    q0 = 1 ; % dimensionless
    Xr10 = 0 ; % dimensionless
    Xr20 = 1 ; % dimensionless
    Xs0 = 0 ; % dimensionless
    Xf0 = 0.1 ; % dimensionless
    g0 = 1 ; % dimensionless
    V0 = -70e-3 ; % V (volts)
    if celltype==1
        Nai = 10; %mM
    else
        Nai = 14.1; %mM
    end
    Cai = 0.0002; %mM
    Casr = 0.3; %mM
    
else % electrical stimulating

    load InitialCondition_4models_1Hz.mat
    statevar_start = InitialCondition.Paci ;
%%% To improve simulation efficiency and ensure model steady state arrival,
%%% proper initial conditions are pre-identified as the last time step 
%%% state variable values following a 1000-second simulation started with
%%% originally specified values, under the same simulation protocol

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Step 4:  ODE Solver
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
statevar_i = statevar_start ;
         
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ORI_globals = [Cm Vc Vsr Nao Ko Cao Ki GCaL GKr GKs Gf arel brel crel Vleak GpCa GbNa GbCa...
    Bufc Bufsr Kbufc Kbufsr Kup KpCa F R T L0 Pkna Ksat KmCa KmNai alpha gamma KmNa Kmk...
    GNa Gto GK1 PNaK KNaCA Vmaxup celltype G_RyR ] ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
options = odeset('MaxStep',1e-3,'InitialStep',2e-5) ;                                         

switch STIM
    case 1 % apply electrical stimuli

        if (n_stimuli > stimulitokeep)
            for i=1:simints-3*stimulitokeep
                [post,posstatevars] = ode15s(@dydt_paci,intervals(i,:),statevar_i,options,Istim(i),ORI_globals) ;
                statevar_i = posstatevars(end,:) ;
                t = post(end) ;
            end
            statevars = statevar_i;
            for i=simints-3*stimulitokeep+1:simints
                [post,posstatevars] = ode15s(@dydt_paci,intervals(i,:),statevar_i,options,Istim(i),ORI_globals) ;
                t = [t;post(2:end)] ;
                statevars = [statevars;posstatevars(2:end,:)] ;
                statevar_i = posstatevars(end,:) ;
            end
            
        else % % when numbertokeep = n_stimuli

            t = 0 ;
            statevars = statevar_i;
            for i=1:simints
                [post,posstatevars] = ode15s(@dydt_paci,intervals(i,:),statevar_i,options,Istim(i),ORI_globals) ;
                t = [t;post(2:end)] ;
                statevars = [statevars;posstatevars(2:end,:)] ;
                statevar_i = posstatevars(end,:) ;
            end
        end
        
    case 0 % spontaneous
        
        [post,posstatevars] = ode15s(@dydt_paci,[0,(tend-secondtokeep)],statevar_i,options,0,ORI_globals) ;
        statevar_i = posstatevars(end,:) ;
        
        [t,statevars] = ode15s(@dydt_paci,[0,secondtokeep],statevar_i,options,0,ORI_globals) ;
    
end
 
t = t*1000 ; % convert into ms
t = t - min(t) ;
V = statevars(:,15) ;   % but need to align all the action potentials together
V = V*1000; %convert the unit to mV
Cai = statevars(:,17);
    
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
    

