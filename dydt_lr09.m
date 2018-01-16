function deriv = dydt_lr09(t,statevar,Id,ORI_globals)



  outputcell = num2cell(ORI_globals) ;
  [F R T RTF FRT Acap Vmyo Vss VJSR VNSR Ko Nao Cao...
   GNa_ GNab GCaT GCab PCa PCa_Na PCa_K gamma_Cao gamma_Cai... 
   gamma_Nao gamma_Nai gamma_Ko gamma_Ki KmCa GK1_ GKr_ GKs_ pKNa GKp_... 
   INaK_ KmNa_NaK KmK_NaK kNCX ksat eta alpha_rel Krel_inf hrel beta_tau Krel_tau...
   IpCa_ KmpCa Vserca Kmserca CaNSR_max tau_transfer TRPNtot KmTRPN CMDNtot KmCMDN CSQNtot KmCSQN... 
   p_m p_h p_j p_d p_f p_b p_g p_xKr p_xs p_release... 
   Vshift_m Vshift_h Vshift_j Vshift_d Vshift_f Vshift_b Vshift_g Vshift_xKr Vshift_xs...  
   Vshift_K1 Vshift_RKr Vshift_Kp Vshift_NaK Vshift_NCX] = deal(outputcell{:}) ;


 
statevarcell = num2cell(statevar) ;

[V,Cai,CaNSR,CaJSR,Nai,Ki,m,h,j,d,f,b,g,xKr,xs1,xs2,Jrel] = ...
  deal(statevarcell{:}) ;

% Reversal potentials
ENa = RTF*log(Nao/Nai) ;
EK = RTF*log(Ko/Ki) ;
EKs = RTF*log((Ko + pKNa*Nao)/(Ki + pKNa*Nai)) ;
ECa = 0.5*RTF*log(Cao/Cai) ;

% %% compute ionic currents
Vold = V ;

% Na currents
INa = GNa_*m^3*h*j*(V - ENa) ;
INab = GNab*(V - ENa) ;

% % L-type Ca current
ICa_ = PCa*4*F*FRT*V* ...
  (gamma_Cai*Cai*exp(2*V*FRT) - gamma_Cao*Cao)/ ...
  (exp(2*V*FRT) - 1) ;
ICaK_ = PCa_K*F*FRT*V* ...
  (gamma_Ki*Ki*exp(V*FRT) - gamma_Ko*Ko)/ ...
  (exp(V*FRT) - 1) ;
ICaNa_ = PCa_Na*F*FRT*V* ...
  (gamma_Nai*Nai*exp(V*FRT) - gamma_Nao*Nao)/ ...
  (exp(V*FRT) - 1) ;
fCa = 1/(Cai/KmCa + 1) ;

ICaL = ICa_*d*f*fCa ;
ICaL_K = ICaK_*d*f*fCa ;
ICaL_Na = ICaNa_*d*f*fCa ;

ICab = GCab*(V - ECa) ;

IpCa = IpCa_*Cai/(Cai + KmpCa) ;

% % T-type Ca current 
ICaT = GCaT*b^2*g*(V-ECa) ;

% K currents
V = Vold + Vshift_K1 ;
xK1 = 0.004*(1+exp(0.6987*(V-EK+11.724)))/ ...
  (1+exp(0.6168*(V-EK+4.872))) ;
IK1 = GK1_*sqrt(Ko/5.4)*(V-EK)/(1+xK1) ;
V = Vold ;

V = Vold + Vshift_RKr ;
RKr = 1/(exp((V+9)/22.4) + 1) ;
V = Vold ;

IKr = GKr_*sqrt(Ko/5.4)*xKr*RKr*(V - EK) ;

IKs = GKs_*(1+0.6/((3.8e-5/Cai)^(1.4)+1))*xs1*xs2*(V - EKs) ;

V = Vold + Vshift_Kp ;
Kp = 1/(1+exp((7.488-V)/5.98)) ;
IKp = GKp_*Kp*(V - EK) ;
V = Vold ;

% Pumps and transporters
V = Vold + Vshift_NaK ;
sigma_NaK = (exp(Nao/67.3) - 1)/7 ; 
fNaK = 1/(1 + 0.1245*exp(-0.1*V*FRT) + 0.0365*sigma_NaK*exp(-V*FRT)) ;
INaK = INaK_*fNaK*Ko/( (Ko + KmK_NaK)*(1 + (KmNa_NaK/Nai)^2) ) ;
V = Vold ;

V = Vold + Vshift_NCX ;
INCX = kNCX*exp((eta-1)*V*FRT)*(Nai^3*Cao*exp(V*FRT)-Nao^3*Cai) / ...
  ( 1+ksat*exp((eta-1)*V*FRT)*(Nai^3*Cao*exp(V*FRT)+Nao^3*Cai) ) ;
V = Vold ;

Iion = INa + INab + ICaL + ICaL_Na + ICaL_K + ICab + ICaT + ...
  IpCa + IKr + IKs + IK1 + IKp + INCX + INaK ;

% %% Compute rate constants to update gates
V = Vold + Vshift_h ;
lambda_na=1-1/(1+exp(-(V+40)/0.024));
ah= lambda_na*0.135*exp(-(80+V)/6.8);
bh= (1-lambda_na)/(0.13*(1+exp((V+10.66)/(-11.1)))) + ...
  lambda_na*(3.56*exp(0.079*V)+3.1*1e5*exp(0.35*V));
hinf = ah/(ah + bh) ;

V = Vold + Vshift_j ;
lambda_na=1-1/(1+exp(-(V+40)/0.024));
aj =  lambda_na*(-1.2714e5*exp(0.2444*V)-3.474e-5*exp(-0.04391*V))* ...
  (V+37.78)/ ...
  (1+exp(0.311*(V+79.23)));

bj= (1-lambda_na)*(0.3*exp(-2.535e-7*V)/(1+exp(-0.1*(V+32)))) + ...
 lambda_na*(0.1212*exp(-0.01052*V)/(1+exp(-0.1378*(V+40.14))));
jinf = aj/(aj + bj) ;

V = Vold + Vshift_m ;
am = 0.32*(V+47.13)/(1-exp(-0.1*(V+47.13))) ;
bm = 0.08*exp(-V/11) ;
minf = am/(am + bm) ;
V = Vold ;

lambda_na=1-1/(1+exp(-(V+40)/0.024));
ah= lambda_na*0.135*exp(-(80+V)/6.8);
bh= (1-lambda_na)/(0.13*(1+exp((V+10.66)/(-11.1)))) + ...
  lambda_na*(3.56*exp(0.079*V)+3.1*1e5*exp(0.35*V));
tauh = 1/(p_h*(ah + bh)) ;

aj =  lambda_na*(-1.2714e5*exp(0.2444*V)-3.474e-5*exp(-0.04391*V))* ...
  (V+37.78)/ ...
  (1+exp(0.311*(V+79.23)));
bj= (1-lambda_na)*(0.3*exp(-2.535e-7*V)/(1+exp(-0.1*(V+32)))) + ...
 lambda_na*(0.1212*exp(-0.01052*V)/(1+exp(-0.1378*(V+40.14))));
tauj = 1/(p_j*(aj + bj)) ;

am = 0.32*(V+47.13)/(1-exp(-0.1*(V+47.13))) ;
bm = 0.08*exp(-V/11) ;
taum = 1/(p_m*(am + bm)) ;

dmdt = (minf-m)/taum ;
dhdt = (hinf-h)/tauh ;
djdt = (jinf-j)/tauj ;

% % % 2) Ca current
V = Vold + Vshift_d ;
dinf_0 = 1/(1+exp(-(V+10)/6.24)) ;
dinf_1 = 1/(1+exp(-(V+60)/0.024)) ;
dinf = dinf_0*dinf_1 ;
V = Vold ;

taud = (1/p_d)/(1+exp(-(V+10)/6.24))* ...
  (1-exp(-(V+10)/6.24))/(0.035*(V+10)) ;

V = Vold + Vshift_f ;
finf = 1/(1+exp((V+32)/8))+(0.6)/(1+exp((50-V)/20)) ;
V = Vold ;

tauf = (1/p_f)/(0.0197*exp(-(0.0337*(V+10))^2)+0.02) ;

dddt = (dinf - d)/taud ;
dfdt = (finf - f)/tauf ;

% % T Type current
V = Vold + Vshift_b ;
binf = 1./(1+exp(-(V+14.0)/10.8));
V = Vold ;

taub = (1/p_b)*(3.7+6.1/(1+exp((V+25.0)/4.5))) ;

V = Vold + Vshift_g ;
ginf = 1/(1+exp((V+60.0)/5.6));
V = Vold ;

lambda_g=1-1/(1+exp(-V/0.0024));
taug = (1/p_g)*(lambda_g*(-0.875*V+12.0)+12.0*(1-lambda_g)) ;

dbdt = (binf-b)/taub ;
dgdt = (ginf-g)/taug ;

% % IKr
V = Vold + Vshift_xKr ;
xKrinf = 1/(1+exp(-(V+21.5)/7.5)) ;
V = Vold ;
tauxKr = (1/p_xKr)*(1/(0.00138*(V+14.2)/(1-exp(-0.123*(V+14.2)))+ ...
  0.00061*(V+38.9)/(exp(0.145*(V+38.9))-1)) ) ;
dxKrdt = (xKrinf - xKr)/tauxKr ;

% % IKr

V = Vold + Vshift_xs ;
xsinf = 1/(1+exp(-(V-1.5)/16.7));
V = Vold ;

tauxs1 = (10000/p_xs)/(0.719*(V+30)/(1-exp(-0.148*(V+30)))+1.31*(V+30)/ ...
  (exp(0.0687*(V+30))-1));
tauxs2 = 4*tauxs1 ;

dxs1dt = (xsinf - xs1)/tauxs1 ;
dxs2dt = (xsinf - xs2)/tauxs2 ;

% % Intracellular Ca fluxes
% % SR Ca release, uptake, and leak
% 
Jrelinf = alpha_rel*beta_tau*ICaL/((Krel_inf/CaJSR)^hrel + 1) ;
if (abs(Jrelinf) < 1e-5)
  Jrelinf = 0 ;
end
tau_rel = (1/p_release)*beta_tau/(Krel_tau/CaJSR + 1) ;
if (tau_rel < 0.1)
  tau_rel = 0.1 ;
end
  
dJreldt = - (Jrelinf + Jrel)/tau_rel ;
% if (abs(dJreldt) > 100)
%   dJreldt = 0 ;
%   Jrel = -Jrelinf ;
% end

Jserca = Vserca*(Cai/(Cai+Kmserca) - CaNSR/CaNSR_max ) ;

Jtr = (CaNSR-CaJSR)/tau_transfer ;

% % % Buffering factors for rapid buffering approximation

BJSR = (1 + CSQNtot*KmCSQN/(KmCSQN + CaJSR)^2)^-1 ;

Bi = 1/(1+ CMDNtot*KmCMDN/(Cai+KmCMDN)^2+ ...
  TRPNtot*KmTRPN/(Cai+KmTRPN)^2) ;

% % % Derivatives for voltage and ionic concentrations
dVdt = -(Id + Iion) ;

dNai = -(INa + INab + ICaL_Na + 3*INCX + 3*INaK)*Acap/(Vmyo*F) ;
dKi = -(IKr + IKs + IK1 + IKp + ICaL_K - 2*INaK)*Acap/(Vmyo*F) ;

dCai = Bi*( -Jserca*VNSR/Vmyo + Jrel*VJSR/Vmyo - ...
  (ICaL + ICaT + ICab + IpCa - 2*INCX)*Acap/(2*Vmyo*F) ) ;
dCaJSR = BJSR*(Jtr - Jrel) ;
dCaNSR = Jserca - Jtr*VJSR/VNSR ;

deriv = [dVdt;dCai;dCaNSR;dCaJSR;dNai;dKi; ...
  dmdt;dhdt;djdt;dddt;dfdt;dbdt;dgdt; ... 
  dxKrdt;dxs1dt;dxs2dt;dJreldt] ;

return


