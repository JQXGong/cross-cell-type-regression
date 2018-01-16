% MATLAB Implementation of the O'Hara-Rudy dynamic (ORd) model for the
% undiseased human ventricular action potential and calcium transient
%
% The ORd model is described in the article "Simulation of the Undiseased
% Human Cardiac Ventricular Action Potential: Model Formulation and
% Experimental Validation"
% by Thomas O'Hara, Laszlo Virag, Andras Varro, and Yoram Rudy
%
% The article and supplemental materails are freely available in the
% Open Access journal PLoS Computational Biology
% Link to Article:
% http://www.ploscompbiol.org/article/info:doi/10.1371/journal.pcbi.1002061
% 
% Email: tom.ohara@gmail.com / rudy@wustl.edu
% Web: http://rudylab.wustl.edu
% 
% The ORd model is free software: you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the
% Free Software Foundation, either version 3 of the License, or (at your
% option) any later version. The ORd model is distributed in the hope that
% it will be useful, but WITHOUT ANY WARRANTY; without even the implied
% warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details (www.gnu.org/licenses)

clear

celltype=0; %endo = 0, epi = 1, M = 2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Step 1:  Define constants that will not be varied randomly
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% In order to comply with potential paralell computing, reverse the usual
%%% global variable declearation for all 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%extracellular ionic concentrations
%%% global Nao Cao Ko;
Nao=140.0;
Cao=1.8;
Ko=5.4;

%physical constants
%%% global R T F Cm;
R=8314.0;
T=310.0;
F=96485.0;
Cm=1.0; %uF

%cell geometry
%%% global L rad vcell Ageo Acap vmyo vnsr vjsr vss;
L=0.01;
rad=0.0011;
vcell=1000*3.14*rad*rad*L;
Ageo=2*3.14*rad*rad+2*3.14*rad*L;
Acap=2*Ageo;
vmyo=0.68*vcell;
vnsr=0.0552*vcell;
vjsr=0.0048*vcell;
vss=0.02*vcell;
 
%CaMK constants
%%% global aCaMK bCaMK CaMKo;
aCaMK=0.05;
bCaMK=0.00068;
CaMKo=0.05;

%reverse potential constant
%%% global PKNa;
PKNa=0.01833;

%INa constants
%%% global Ahf Ahs;
Ahf=0.99;
Ahs=1.0-Ahf;

%INaL constants
%%% global thL;
thL=200.0;

%Ito constants
%%% global delta_epi;
if celltype==1
    delta_epi=1.0-(0.95/(1.0+exp((V+70.0)/5.0)));
else
    delta_epi=1.0;
end

%ICaL, ICaNa, ICaK constants
%%% global Aff tjca Kmn k2n zca;
Aff=0.6;
tjca=75.0;
Kmn=0.002;
k2n=1000.0;
zca=2.0;

%INaCa_i constants
%%% global kna1 kna2 kna3 kasymm wna wca wnaca kcaon kcaoff qna qca zna;
kna1=15.0;
kna2=5.0;
kna3=88.12;
kasymm=12.5;
wna=6.0e4;
wca=6.0e4;
wnaca=5.0e3;
kcaon=1.5e6;
kcaoff=5.0e3;
qna=0.5224;
qca=0.1670;
zna=1.0;

%INaCa_ss constants
%%% global KmCaAct;
KmCaAct=150.0e-6;

%INaK constants
%%% global k1p k1m k2p k2m k3p k3m k4p k4m Knai0 Knao0 delta;
k1p=949.5;
k1m=182.4;
k2p=687.2;
k2m=39.4;
k3p=1899.0;
k3m=79300.0;
k4p=639.0;
k4m=40.0;
Knai0=9.073;
Knao0=27.78;
delta=-0.1550;
%%% global Kki Kko MgADP MgATP Kmgatp H eP Khp Knap Kxkur; 
Kki=0.5;
Kko=0.3582;
MgADP=0.05;
MgATP=9.8;
Kmgatp=1.698e-7;
H=1.0e-7;
eP=4.2;
Khp=1.698e-7;
Knap=224.0;
Kxkur=292.0;
%%% global zk;
zk=1.0;

%calcium buffer constants
%%% global cmdnmax kmcmdn trpnmax kmtrpn BSRmax KmBSR BSLmax KmBSL csqnmax kmcsqn;
cmdnmax=0.05;
if celltype==1
    cmdnmax=cmdnmax*1.3;
end
kmcmdn=0.00238;
trpnmax=0.07;
kmtrpn=0.0005;
BSRmax=0.047;
KmBSR=0.00087;
BSLmax=1.124;
KmBSL=0.0087;
csqnmax=10.0;
kmcsqn=0.8;

%jsr constants
%%% global bt a_rel;
bt=4.75;
a_rel=0.5*bt;

% computed quantities that do not change during simulation
%%% global GNa GNaL Gto GKr GKs GK1 Gncx GKb GpCa PCa Pnak PNab PCab KmCaMK KmCaM;

GNa=75;

GNaL=0.0075;
if celltype==1
    GNaL=GNaL*0.6;
end

Gto=0.02;
if celltype==1
    Gto=Gto*4.0;
elseif celltype==2
    Gto=Gto*4.0;
end

GKr=0.046;
if celltype==1
    GKr=GKr*1.3;
elseif celltype==2
    GKr=GKr*0.8;
end

GKs=0.0034;
GK1=0.1908;
if celltype==1
    GKs=GKs*1.4;
end
if celltype==1
    GK1=GK1*1.2;
elseif celltype==2
    GK1=GK1*1.3;
end

Gncx=0.0008;                %GNaCa
if celltype==1
    Gncx=Gncx*1.1;
elseif celltype==2
    Gncx=Gncx*1.4;
end

GKb=0.003;
if celltype==1
    GKb=GKb*0.6;
end

GpCa=0.0005;

PCa=0.0001;
if celltype==1
    PCa=PCa*1.2;
elseif celltype==2
    PCa=PCa*2.5;
end

Pnak=30;
if celltype==1
    Pnak=Pnak*0.9;
elseif celltype==2
    Pnak=Pnak*0.7;
end

PNab=3.75e-10;
PCab=2.5e-8;

KmCaMK=0.15;
KmCaM=0.0015;

%%% global SERCA_total RyR_total
SERCA_total = 1 ;
RyR_total = 1 ;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Step 2:  Define simulation, stimulus, and recording parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

stiminterval = 1000 ;  % interval bewteen stimuli, ms
n_stimuli = 120 ; % number of stimuli
tend =  n_stimuli*stiminterval ; % end of simulation, ms
stimdelay = 100 ; % delay between start of each cycle and beginning of current injection, ms
stimdur = 0.5 ; % length of current injection, ms
stim_amp = 80 ; % amplitude of current injection, nA

stim_starts = stimdelay + stiminterval*(0:n_stimuli-1)  ;
stim_ends = stim_starts + stimdur ;

simints = 3*n_stimuli ;
for i=1:n_stimuli
  intervals(3*i-2,:) = [stiminterval*(i-1),stim_starts(i)] ; % time during delay
  intervals(3*i-1,:) = [stim_starts(i),stim_ends(i)] ;       % time during stimulus   
  intervals(3*i,:) = [stim_ends(i),stiminterval*i] ;         % after stimulus   
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
%%% Original valus as specified in the O'Hara et al. article
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% %     V=-87;
% %     Nai=7;
% %     Nass=Nai;
% %     Ki=145;
% %     Kss=Ki;
% %     Cai=1.0e-4;
% %     Cass=Cai;
% %     Cansr=1.2;
% %     Cajsr=Cansr;
% %     m=0;
% %     hf=1;
% %     hs=1;
% %     j=1;
% %     hsp=1;
% %     jp=1;
% %     mL=0;
% %     hL=1;
% %     hLp=1;
% %     a=0;
% %     iF=1;
% %     iS=1;
% %     ap=0;
% %     iFp=1;
% %     iSp=1;
% %     d=0;
% %     ff=1;
% %     fs=1;
% %     fcaf=1;
% %     fcas=1;
% %     jca=1;
% %     nca=0;
% %     ffp=1;
% %     fcafp=1;
% %     xrf=0;
% %     xrs=0;
% %     xs1=0;
% %     xs2=0;
% %     xk1=1;
% %     Jrelnp=0;
% %     Jrelp=0;
% %     CaMKt=0;

load InitialCondition_4models_1Hz.mat
statevar_start = InitialCondition.ORd ;
%%% To improve simulation efficiency and ensure model steady state arrival,
%%% proper initial conditions are pre-identified as the last time step 
%%% state variable values following a 1000-second simulation started with
%%% originally specified values, under the same simulation protocol
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Step 4:  ODE Solver
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            
statevar_i = statevar_start ;
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ORI_globals = [Nao Cao Ko R T F Cm Acap vmyo vnsr vjsr vss aCaMK bCaMK CaMKo...
    PKNa Ahf Ahs thL delta_epi Aff tjca Kmn k2n zca kna1 kna2 kna3...
    kasymm wna wca wnaca kcaon kcaoff qna qca zna KmCaAct...
    k1p k1m k2p k2m k3p k3m k4p k4m Knai0 Knao0 delta...
    Kki Kko MgADP MgATP Kmgatp H eP Khp Knap Kxkur zk...
    cmdnmax kmcmdn trpnmax kmtrpn BSRmax KmBSR BSLmax KmBSL csqnmax kmcsqn...
    bt a_rel...
    GNa GNaL Gto GKr GKs GK1 Gncx GKb GpCa PCa Pnak PNab PCab KmCaMK KmCaM...
    SERCA_total RyR_total celltype]  ;
 

options = [] ;
if (n_stimuli > numbertokeep)
    for i=1:simints-3*numbertokeep
        [post,posstatevars] = ode15s(@dydt_ohara,intervals(i,:),statevar_i,options,Istim(i),ORI_globals) ;
        statevar_i = posstatevars(end,:) ;
        t = post(end) ;
    end
    statevars = statevar_i;
    for i=simints-3*numbertokeep+1:simints
        [post,posstatevars] = ode15s(@dydt_ohara,intervals(i,:),statevar_i,options,Istim(i),ORI_globals) ;
        t = [t;post(2:end)] ;
        statevars = [statevars;posstatevars(2:end,:)] ;
        statevar_i = posstatevars(end,:) ;
    end
        
else % n_stimuli = numbertokeep
    t = 0 ;
    statevars = statevar_i;
    for i=1:simints
        [post,posstatevars] = ode15s(@dydt_ohara,intervals(i,:),statevar_i,options,Istim(i),ORI_globals) ;
        t = [t;post(2:end)] ;
        statevars = [statevars;posstatevars(2:end,:)] ;
        statevar_i = posstatevars(end,:) ;
    end
end
    t = t - min(t) ;

    outputcell = num2cell(statevars,1) ;

    [V Nai Nass Ki Kss Cai Cass Cansr Cajsr m hf hs j hsp jp mL hL...
        hLp a iF iS ap iFp iSp d ff fs fcaf fcas jca nca ffp fcafp...
        xrf xrs xs1 xs2 xk1 Jrelnp Jrelp CaMKt] = deal(outputcell{:}) ;
    
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


