
function deriv = dydt_shannon04(t,statevar,Id,ORI_globals)

outputcell = num2cell(ORI_globals) ;
[ Bmax_Calsequestrin Bmax_SLB_SL Bmax_SLB_jct Bmax_SLHigh_SL Bmax_SLHigh_jct... 
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
  Vshift_X_tos, Vshift_Y_tos ] = deal(outputcell{:}) ;


Y = statevar ;

%-------------------------------------------------------------------------------
% State variables
%-------------------------------------------------------------------------------

% 1: Ca_Calsequestrin (millimolar) (in Ca_buffer)
% 2: Ca_SL (millimolar) (in Ca_buffer)
% 3: Ca_SLB_SL (millimolar) (in Ca_buffer)
% 4: Ca_SLB_jct (millimolar) (in Ca_buffer)
% 5: Ca_SLHigh_SL (millimolar) (in Ca_buffer)
% 6: Ca_SLHigh_jct (millimolar) (in Ca_buffer)
% 7: Ca_SR (millimolar) (in Ca_buffer)
% 8: Ca_jct (millimolar) (in Ca_buffer)
% 9: Cai (millimolar) (in Ca_buffer)
% 10: d (dimensionless) (in ICaL_d_gate)
% 11: fCaB_SL (dimensionless) (in ICaL_fCa_gate)
% 12: fCaB_jct (dimensionless) (in ICaL_fCa_gate)
% 13: f (dimensionless) (in ICaL_f_gate)
% 14: Xr (dimensionless) (in IKr_Xr_gate)
% 15: Xs (dimensionless) (in IKs_Xs_gate)
% 16: h (dimensionless) (in INa_h_gate)
% 17: j (dimensionless) (in INa_j_gate)
% 18: m (dimensionless) (in INa_m_gate)
% 19: X_tof (dimensionless) (in Itof_X_gate)
% 20: Y_tof (dimensionless) (in Itof_Y_gate)
% 21: R_tos (dimensionless) (in Itos_R_gate)
% 22: X_tos (dimensionless) (in Itos_X_gate)
% 23: Y_tos (dimensionless) (in Itos_Y_gate)
% 24: I (dimensionless) (in Jrel_SR)
% 25: O (dimensionless) (in Jrel_SR)
% 26: R1 (dimensionless) (R in Jrel_SR)
% 27: Na_SL (millimolar) (in Na_buffer)
% 28: Na_SL_buf (millimolar) (in Na_buffer)
% 29: Na_jct (millimolar) (in Na_buffer)
% 30: Na_jct_buf (millimolar) (in Na_buffer)
% 31: Nai (millimolar) (in Na_buffer)
% 32: V (millivolt) (in cell)
% 33: Ca_Calmodulin (millimolar) (in cytosolic_Ca_buffer)
% 34: Ca_Myosin (millimolar) (in cytosolic_Ca_buffer)
% 35: Ca_SRB (millimolar) (in cytosolic_Ca_buffer)
% 36: Ca_TroponinC (millimolar) (in cytosolic_Ca_buffer)
% 37: Ca_TroponinC_Ca_Mg (millimolar) (in cytosolic_Ca_buffer)
% 38: Mg_Myosin (millimolar) (in cytosolic_Ca_buffer)
% 39: Mg_TroponinC_Ca_Mg (millimolar) (in cytosolic_Ca_buffer)


Vold = Y(32) ;

% geometry
Vol_Cell = 3.141592654*(cell_radius/1000.0)^2.0*cell_length/1000.0^3.0;
Vol_cytosol = 0.65*Vol_Cell;
Vol_SR = 0.035*Vol_Cell;
% Calsequestrin 
dCalsequestrin = kon_Calsequestrin*Y(7)*(Bmax_Calsequestrin*Vol_cytosol/Vol_SR-Y(1))-koff_Calsequestrin*Y(1);
dY(1, 1) = dCalsequestrin; 
% SL Ca buffering
Vol_SL = 0.02*Vol_Cell;
dCa_SLB_SL = kon_SL*Y(2)*(Bmax_SLB_SL*Vol_cytosol/Vol_SL-Y(3))-koff_SLB*Y(3);
Vol_jct = 0.00051*Vol_Cell;
dCa_SLB_jct = kon_SL*Y(8)*(Bmax_SLB_jct*0.1*Vol_cytosol/Vol_jct-Y(4))-koff_SLB*Y(4);
dCa_SLHigh_SL = kon_SL*Y(2)*(Bmax_SLHigh_SL*Vol_cytosol/Vol_SL-Y(5))-koff_SLHigh*Y(5);
dCa_SLHigh_jct = kon_SL*Y(8)*(Bmax_SLHigh_jct*0.1*Vol_cytosol/Vol_jct-Y(6))-koff_SLHigh*Y(6);

dY(3, 1) = dCa_SLB_SL;
dY(4, 1) = dCa_SLB_jct;
dY(5, 1) = dCa_SLHigh_SL;
dY(6, 1) = dCa_SLHigh_jct;

dCa_jct_tot_bound = dCa_SLB_jct+dCa_SLHigh_jct;
dCa_SL_tot_bound = dCa_SLB_SL+dCa_SLHigh_SL;
Q_CaL = Q10_CaL^((T-310.0)/10.0);
temp = 0.45*Y(10)*Y(13)*Q_CaL*Y(32)*F^2.0/(R2*T);
fCa_jct = 1.0-Y(12);
i_CaL_Ca_jct = temp*fCa_jct*Fx_ICaL_jct*PCa*4.0*(gamma_Cai*Y(8)*exp(2.0*Y(32)*F/(R2*T))-gamma_Cao*Cao)/(exp(2.0*Y(32)*F/(R2*T))-1.0);
Ka_jct = 1.0/(1.0+(Kd_act/Y(8))^3.0);
Q_NCX = Q10_NCX^((T-310.0)/10.0);
temp_jct = (exp(eta*Y(32)*F/(R2*T))*Y(29)^HNa*Cao-exp((eta-1.0)*Y(32)*F/(R2*T))*Nao^HNa*Y(8))/(1.0+ksat*exp((eta-1.0)*Y(32)*F/(R2*T)));
i_NaCa_jct = Fx_NCX_jct*V_max_2*Ka_jct*Q_NCX*temp_jct/(K_mCai*Nao^HNa*(1.0+(Y(29)/K_mNai)^HNa)+K_mNao^HNa*Y(8)*(1.0+Y(8)/K_mCai)+K_mCao*Y(29)^HNa+Y(29)^HNa*Cao+Nao^HNa*Y(8));
E_Ca_jct = R2*T/(2.0*F)*log(Cao/Y(8));
i_Cab_jct = G_CaBk*Fx_CaBk_jct*(Y(32)-E_Ca_jct);
Q_SLCaP = Q10_SLCaP^((T-310.0)/10.0);
i_Cap_jct = Q_SLCaP*V_maxAF*Fx_SLCaP_jct/(1.0+(Km/Y(8))^H1);
i_Ca_jct_tot = i_CaL_Ca_jct-2.0*i_NaCa_jct+i_Cab_jct+i_Cap_jct;
fCa_SL = 1.0-Y(11);
i_CaL_Ca_SL = temp*fCa_SL*Fx_ICaL_SL*PCa*4.0*(gamma_Cai*Y(2)*exp(2.0*Y(32)*F/(R2*T))-gamma_Cao*Cao)/(exp(2.0*Y(32)*F/(R2*T))-1.0);
Ka_SL = 1.0/(1.0+(Kd_act/Y(2))^3.0);
temp_SL = (exp(eta*Y(32)*F/(R2*T))*Y(27)^HNa*Cao-exp((eta-1.0)*Y(32)*F/(R2*T))*Nao^HNa*Y(2))/(1.0+ksat*exp((eta-1.0)*Y(32)*F/(R2*T)));
i_NaCa_SL = Fx_NCX_SL*V_max_2*Ka_SL*Q_NCX*temp_SL/(K_mCai*Nao^HNa*(1.0+(Y(27)/K_mNai)^HNa)+K_mNao^HNa*Y(2)*(1.0+Y(2)/K_mCai)+K_mCao*Y(27)^HNa+Y(27)^HNa*Cao+Nao^HNa*Y(2));
E_Ca_SL = R2*T/(2.0*F)*log(Cao/Y(2));
i_Cab_SL = G_CaBk*Fx_CaBk_SL*(Y(32)-E_Ca_SL);
i_Cap_SL = Q_SLCaP*V_maxAF*Fx_SLCaP_SL/(1.0+(Km/Y(2))^H1);
i_Ca_SL_tot = i_CaL_Ca_SL-2.0*i_NaCa_SL+i_Cab_SL+i_Cap_SL;
Q_SRCaP = Q10_SRCaP^((T-310.0)/10.0);
j_pump_SR = Q_SRCaP*V_max_3*Vol_cytosol/Vol_SR*((Y(9)/Kmf)^H2-(Y(7)/Kmr)^H2)/(1.0+(Y(9)/Kmf)^H2+(Y(7)/Kmr)^H2);
j_leak_SR = KSRleak*(Y(7)-Y(8));
j_rel_SR = ks*Y(25)*(Y(7)-Y(8));

dY(7, 1) = j_pump_SR-(j_leak_SR*Vol_cytosol/Vol_SR+j_rel_SR)-dCalsequestrin; % Ca_SR

Cm = Cm_per_area*2.0*cell_radius/10000.0*pi*cell_length/10000.0;
J_Ca_jct_SL = (Y(8)-Y(2))*8.2413e-13;

dY(8, 1) = -i_Ca_jct_tot*Cm/(Vol_jct*2.0*F)-J_Ca_jct_SL/Vol_jct+j_rel_SR*Vol_SR/Vol_jct+j_leak_SR*Vol_cytosol/Vol_jct-1.0*dCa_jct_tot_bound; % Ca_jct

J_Ca_SL_cytosol = (Y(2)-Y(9))*3.7243e-12;

dY(2, 1) = -i_Ca_SL_tot*Cm/(Vol_SL*2.0*F)+(J_Ca_jct_SL-J_Ca_SL_cytosol)/Vol_SL-1.0*dCa_SL_tot_bound; % Ca_SL

dCa_TroponinC = kon_TroponinC*Y(9)*(Bmax_TroponinC-Y(36))-koff_TroponinC*Y(36);
dCa_TroponinC_Ca_Mg = kon_TroponinC_Ca_Mg_Ca*Y(9)*(Bmax_TroponinC_Ca_Mg_Ca-(Y(37)+Y(39)))-koff_TroponinC_Ca_Mg_Ca*Y(37);
dMg_TroponinC_Ca_Mg = kon_TroponinC_Ca_Mg_Mg*Mgi*(Bmax_TroponinC_Ca_Mg_Mg-(Y(37)+Y(39)))-koff_TroponinC_Ca_Mg_Mg*Y(39);
dCa_Calmodulin = kon_Calmodulin*Y(9)*(Bmax_Calmodulin-Y(33))-koff_Calmodulin*Y(33);
dCa_Myosin = kon_Myosin_Ca*Y(9)*(Bmax_Myosin_Ca-(Y(34)+Y(38)))-koff_Myosin_Ca*Y(34);
dMg_Myosin = kon_Myosin_Mg*Mgi*(Bmax_Myosin_Mg-(Y(34)+Y(38)))-koff_Myosin_Mg*Y(38);
dCa_SRB = kon_SRB*Y(9)*(Bmax_SRB-Y(35))-koff_SRB*Y(35);
dCa_cytosol_tot_bound = dCa_TroponinC+dCa_TroponinC_Ca_Mg+dMg_TroponinC_Ca_Mg+dCa_Calmodulin+dCa_Myosin+dMg_Myosin+dCa_SRB;

dY(9, 1) = -j_pump_SR*Vol_SR/Vol_cytosol+J_Ca_SL_cytosol/Vol_cytosol-1.0*dCa_cytosol_tot_bound; % Cai

i_CaL_Na_jct = temp*fCa_jct*Fx_ICaL_jct*PNa*(gamma_Nai*Y(29)*exp(Y(32)*F/(R2*T))-gamma_Nao*Nao)/(exp(Y(32)*F/(R2*T))-1.0);
i_CaL_Na_SL = temp*fCa_SL*Fx_ICaL_SL*PNa*(gamma_Nai*Y(27)*exp(Y(32)*F/(R2*T))-gamma_Nao*Nao)/(exp(Y(32)*F/(R2*T))-1.0);
i_CaL_K = temp*(fCa_SL*Fx_ICaL_SL+fCa_jct*Fx_ICaL_jct)*PK*(gamma_Ki*Ki*exp(Y(32)*F/(R2*T))-gamma_Ko*Ko)/(exp(Y(32)*F/(R2*T))-1.0);
i_CaL = i_CaL_Ca_SL+i_CaL_Ca_jct+i_CaL_Na_SL+i_CaL_Na_jct+i_CaL_K;

% 
Y(32) = Vold + Vshift_d ;
d_infinity = 1.0/(1.0+exp(-(Y(32)+14.5)/6.0));
Y(32) = Vold ;

tau_d = p_d*( 1.0*d_infinity*(1.0-exp(-(Y(32)+14.5)/6.0))/(0.035*(Y(32)+14.5)) ) ;
dY(10, 1) = (d_infinity-Y(10))/tau_d;
dY(11, 1) = 1.7*Y(2)*(1.0-Y(11))-11.9e-3*Y(11);
dY(12, 1) = 1.7*Y(8)*(1.0-Y(12))-11.9e-3*Y(12);

Y(32) = Vold + Vshift_f ;
f_infinity = 1.0/(1.0+exp((Y(32)+35.06)/3.6))+0.6/(1.0+exp((50.0-Y(32))/20.0));
Y(32) = Vold ;

tau_f = p_f*( 1.0/(0.0197*exp(-(0.0337*(Y(32)+14.5))^2.0)+0.02) ) ;

dY(13, 1) = (f_infinity-Y(13))/tau_f;

i_Cab = i_Cab_SL+i_Cab_jct;
i_Cap = i_Cap_jct+i_Cap_SL;
E_Cl = -R2*T/F*log(Clo/Cli);
i_Cl_Ca = G_Cl*(Y(32)-E_Cl)*(Fx_Cl_jct/(1.0+Kd_ClCa/Y(8))+Fx_Cl_SL/(1.0+Kd_ClCa/Y(2)));
i_Clb = G_ClBk*(Y(32)-E_Cl);

G_K1 = GK1_*sqrt(Ko/5.4);
E_K = R2*T/F*log(Ko/Ki);

Y(32) = Vold + Vshift_K1 ;
alpha_K1 = 1.02/(1.0+exp(0.2385*(Y(32)-(E_K+59.215))));
beta_K1 = (0.49124*exp(0.08032*(Y(32)-E_K+5.476))+1.0*exp(0.06175*(Y(32)-(E_K+594.31))))/(1.0+exp(-0.5143*(Y(32)-E_K+4.753)));
K1_infinity = alpha_K1/(alpha_K1+beta_K1);
Y(32) = Vold ;

i_K1 = G_K1*K1_infinity*(Y(32)-E_K);

G_IKr = GKr_*sqrt(Ko/5.4);
Rr = 1.0/(1.0+exp((33.0+Y(32))/22.4));
i_Kr = G_IKr*Y(14)*Rr*(Y(32)-E_K);

Y(32) = Vold + Vshift_Xr ;
Xr_infinity = 1.0/(1.0+exp(-(50.0+Y(32))/7.5));
Y(32) = Vold ;

tau_Xr = p_Xr*( 1.0/(0.00138*(Y(32)+7.0)/(1.0-exp(-0.123*(Y(32)+7.0)))+0.00061*(Y(32)+10.0)/(exp(0.145*(Y(32)+10.0))-1.0)) ) ;

dY(14, 1) = (Xr_infinity-Y(14))/tau_Xr;

pCa_jct = -log10(Y(8)/1.0)+3.0;
pCa_SL = -log10(Y(2)/1.0)+3.0;
G_Ks_jct = GKs_*(0.057+0.19/(1.0+exp((-7.2+pCa_jct)/0.6)));
G_Ks_SL = GKs_*(0.057+0.19/(1.0+exp((-7.2+pCa_SL)/0.6)));
E_Ks_jct = R2*T/F*log((Ko+pKNa*Nao)/(Ki+pKNa*Y(29)));
E_Ks_SL = R2*T/F*log((Ko+pKNa*Nao)/(Ki+pKNa*Y(27)));
E_Ks = R2*T/F*log((Ko+pKNa*Nao)/(Ki+pKNa*Y(31)));
i_Ks_jct = Fx_Ks_jct*G_Ks_jct*Y(15)^2.0*(Y(32)-E_Ks);
i_Ks_SL = Fx_Ks_SL*G_Ks_SL*Y(15)^2.0*(Y(32)-E_Ks);
i_Ks = i_Ks_jct+i_Ks_SL;

Y(32) = Vold + Vshift_Xs ;
Xs_infinity = 1.0/(1.0+exp(-(Y(32)-1.5)/16.7));
Y(32) = Vold ;

tau_Xs = p_Xs*( 1.0/(7.19e-5*(Y(32)+30.0)/(1.0-exp(-0.148*(Y(32)+30.0)))+1.31e-4*(Y(32)+30.0)/(-1.0+exp(0.0687*(Y(32)+30.0)))) ) ;

dY(15, 1) = (Xs_infinity-Y(15))/tau_Xs;

openProb = Y(18)^3.0*Y(16)*Y(17);
E_Na_jct = R2*T/F*log(Nao/Y(29));
i_Na_jct = Fx_Na_jct*G_INa*openProb*(Y(32)-E_Na_jct);
E_Na_SL = R2*T/F*log(Nao/Y(27));
i_Na_SL = Fx_Na_SL*G_INa*openProb*(Y(32)-E_Na_SL);
i_Na = i_Na_jct+i_Na_SL;
i_NaCa = i_NaCa_jct+i_NaCa_SL;
Q_NaK = Q10_NaK^((T-310.0)/10.0);
Q_Km_Nai = Q10_Km_Nai^((T-310.0)/10.0);
sigma = (exp(Nao/67.3)-1.0)/7.0;
f_NaK = 1.0/(1.0+0.1245*exp(-0.1*Y(32)*F/(R2*T))+0.0365*sigma*exp(-Y(32)*F/(R2*T)));
i_NaK_jct = Fx_NaK_jct*Q_NaK*I_NaK_max*f_NaK/(1.0+(Q_Km_Nai*Km_Nai/Y(29))^H_NaK)*Ko/(Ko+Km_Ko);
i_NaK_SL = Fx_NaK_SL*Q_NaK*I_NaK_max*f_NaK/(1.0+(Q_Km_Nai*Km_Nai/Y(27))^H_NaK)*Ko/(Ko+Km_Ko);
i_NaK = i_NaK_jct+i_NaK_SL;

Y(32) = Vold + Vshift_h ;
if (Y(32) < -40.0)
   alpha_h = 0.135*exp((80.0+Y(32))/-6.8);
else
   alpha_h = 0.0;
end;

if (Y(32) < -40.0)
   beta_h = 3.56*exp(0.079*Y(32))+3.1e5*exp(0.35*Y(32));
else
   beta_h = 1.0/(0.13*(1.0+exp((Y(32)+10.66)/-11.1)));
end;

if (Y(32) < -40.0)
   alpha_h = 0.135*exp((80.0+Y(32))/-6.8);
else
   alpha_h = 0.0;
end;
h_infinity = alpha_h/(alpha_h+beta_h);
Y(32) = Vold ;

if (Y(32) < -40.0)
   beta_h = 3.56*exp(0.079*Y(32))+3.1e5*exp(0.35*Y(32));
else
   beta_h = 1.0/(0.13*(1.0+exp((Y(32)+10.66)/-11.1)));
end;
tau_h = p_h*( 1.0/(alpha_h+beta_h) ) ;

dY(16, 1) = (h_infinity-Y(16))/tau_h;

Y(32) = Vold + Vshift_j ;
if (Y(32) < -40.0)
   alpha_j = (-1.2714e5*exp(0.2444*Y(32))-3.474e-5*exp(-0.04391*Y(32)))*(Y(32)+37.78)/1.0/(1.0+exp(0.311*(Y(32)+79.23)));
else
   alpha_j = 0.0;
end;

if (Y(32) < -40.0)
   beta_j = 0.1212*exp(-0.01052*Y(32))/(1.0+exp(-0.1378*(Y(32)+40.14)));
else
   beta_j = 0.3*exp(-2.535e-7*Y(32))/(1.0+exp(-0.1*(Y(32)+32.0)));
end;
j_infinity = alpha_j/(alpha_j+beta_j);
Y(32) = Vold ;

if (Y(32) < -40.0)
   alpha_j = (-1.2714e5*exp(0.2444*Y(32))-3.474e-5*exp(-0.04391*Y(32)))*(Y(32)+37.78)/1.0/(1.0+exp(0.311*(Y(32)+79.23)));
else
   alpha_j = 0.0;
end;

if (Y(32) < -40.0)
   beta_j = 0.1212*exp(-0.01052*Y(32))/(1.0+exp(-0.1378*(Y(32)+40.14)));
else
   beta_j = 0.3*exp(-2.535e-7*Y(32))/(1.0+exp(-0.1*(Y(32)+32.0)));
end;
tau_j = p_j*( 1.0/(alpha_j+beta_j) ) ;


dY(17, 1) = (j_infinity-Y(17))/tau_j;

alpha_m = 0.32*(Y(32)+47.13)/1.0/(1.0-exp(-0.1*(Y(32)+47.13)));
beta_m = 0.08*exp(-Y(32)/11.0);
tau_m = p_m*( 1.0/(alpha_m+beta_m) ) ;

Y(32) = Vold + Vshift_m ;
alpha_m = 0.32*(Y(32)+47.13)/1.0/(1.0-exp(-0.1*(Y(32)+47.13)));
beta_m = 0.08*exp(-Y(32)/11.0);
m_infinity = alpha_m/(alpha_m+beta_m);
Y(32) = Vold ;

dY(18, 1) = (m_infinity-Y(18))/tau_m;

i_Nab_jct = Fx_NaBk_jct*G_NaBk*(Y(32)-E_Na_jct);
i_Nab_SL = Fx_NaBk_SL*G_NaBk*(Y(32)-E_Na_SL);
i_Nab = i_Nab_jct+i_Nab_SL;
i_tof = G_tof*Y(19)*Y(20)*(Y(32)-E_K);

Y(32) = Vold + Vshift_X_tof ;
X_tof_infinity = 1.0/(1.0+exp(-(Y(32)+3.0)/15.0));
Y(32) = Vold ;

tau_X_tof = p_X_tof*( 3.5*exp(-(Y(32)/30.0)^2.0)+1.5 ) ;

dY(19, 1) = (X_tof_infinity-Y(19))/tau_X_tof;

Y(32) = Vold + Vshift_Y_tof ;
Y_tof_infinity = 1.0/(1.0+exp((Y(32)+33.5)/10.0));
Y(32) = Vold ;

tau_Y_tof = p_Y_tof*( 20.0/(1.0+exp((Y(32)+33.5)/10.0))+20.0 ) ;

dY(20, 1) = (Y_tof_infinity-Y(20))/tau_Y_tof;

i_tos = G_tos*Y(22)*(Y(23)+0.5*Y(21))*(Y(32)-E_K);

Y(32) = Vold + Vshift_R_tos ;
R_tos_infinity = 1.0/(1.0+exp((Y(32)+33.5)/10.0));
Y(32) = Vold ;

tau_R_tos = p_R_tos*( 2.8e3/(1.0+exp((Y(32)+60.0)/10.0))+220.0 ) ;
R_tos_other = R_tos_infinity;

dY(21, 1) = (R_tos_infinity-Y(21))/tau_R_tos;

Y(32) = Vold + Vshift_X_tos ;
X_tos_infinity = 1.0/(1.0+exp(-(Y(32)+3.0)/15.0));
Y(32) = Vold ;

tau_X_tos = p_X_tos*( 9.0/(1.0+exp((Y(32)+3.0)/15.0))+0.5 ) ;

dY(22, 1) = (X_tos_infinity-Y(22))/tau_X_tos;

Y(32) = Vold + Vshift_Y_tos ;
Y_tos_infinity = 1.0/(1.0+exp((Y(32)+33.5)/10.0));
Y(32) = Vold ;

tau_Y_tos = p_Y_tos*( 3000.0/(1.0+exp((Y(32)+60.0)/10.0))+30.0 ) ;

dY(23, 1) = (Y_tos_infinity-Y(23))/tau_Y_tos;

kCaSR = Max_SR-(Max_SR-Min_SR)/(1.0+(EC50_SR/Y(7))^HSR);
koSRCa = koCa/kCaSR;
kiSRCa = kiCa*kCaSR;
RI = 1.0-Y(26)-Y(25)-Y(24);

dY(26, 1) = kim*RI-kiSRCa*Y(8)*Y(26)-(koSRCa*Y(8)^2.0*Y(26)-kom*Y(25));
dY(25, 1) = koSRCa*Y(8)^2.0*Y(26)-kom*Y(25)-(kiSRCa*Y(8)*Y(25)-kim*Y(24));
dY(24, 1) = kiSRCa*Y(8)*Y(25)-kim*Y(24)-(kom*Y(24)-koSRCa*Y(8)^2.0*RI);

dNa_jct_buf = kon*Y(29)*(Bmax_jct-Y(30))-koff*Y(30);
dNa_SL_buf = kon*Y(27)*(Bmax_SL-Y(28))-koff*Y(28);

dY(30, 1) = dNa_jct_buf;
dY(28, 1) = dNa_SL_buf;

J_Na_jct_SL = (Y(29)-Y(27))*1.8313e-14;

dY(29, 1) = -Cm*(i_Na_jct+3.0*i_NaCa_jct+i_Nab_jct+3.0*i_NaK_jct+i_CaL_Na_jct)/(Vol_jct*F)-J_Na_jct_SL/Vol_jct-dNa_jct_buf;

J_Na_SL_cytosol = (Y(27)-Y(31))*1.6386e-12;

dY(27, 1) = -Cm*(i_Na_SL+3.0*i_NaCa_SL+i_Nab_SL+3.0*i_NaK_SL+i_CaL_Na_SL)/(Vol_SL*F)+(J_Na_jct_SL-J_Na_SL_cytosol)/Vol_SL-dNa_SL_buf;
dY(31, 1) = J_Na_SL_cytosol/Vol_cytosol;

dY(32, 1) = -(i_Na+i_Nab+i_NaK+i_Kr+i_Ks+i_tos+i_tof+i_K1+i_NaCa+i_Cl_Ca+i_Clb+i_CaL+i_Cab+i_Cap+ Id);
dY(36, 1) = dCa_TroponinC;
dY(37, 1) = dCa_TroponinC_Ca_Mg;
dY(39, 1) = dMg_TroponinC_Ca_Mg;
dY(33, 1) = dCa_Calmodulin;
dY(34, 1) = dCa_Myosin;
dY(38, 1) = dMg_Myosin;
dY(35, 1) = dCa_SRB;

deriv=dY;

return 

