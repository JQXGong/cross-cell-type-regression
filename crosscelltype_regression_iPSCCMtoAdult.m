
%%% Implementation of cross-cell type modeling approach as described in the article 
%%% "Population-based mechanistic modeling allows for quantitative predictions of drug responses across cell types"
%%% by Jingqi Q.X. Gong and Eric A. Sobie

%%% The article and supplementary materials are freely available in the
%%% Open Access journal Nature Systems Biology and Applications
%%% Link to the article :
%%% LINK

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Required scripts and dataset:
%%% PLS_nipals.m -- Performs partial least-squares regression
%%% zscore.m -- Performs z-normalization
%%% Metrics_humanadult.mat -- Metrics of adult human AP and CaT + outputnames
%%% Metrics_iPSCCM.mat -- Metrics of iPSC-CM AP and CaT under 5 conditions + inputnames

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear 

%% Loading and formalizing the dataset

load Metrics_humanadult.mat 
%%% Metrics of adult human AP and CaT + outputnames
load Metrics_iPSCCM.mat
%%% Metrics of iPSC-CM AP and CaT under 5 conditions + inputnames


%%% Simulation results from 5 conditions in iPSC-CM are used, in the orther of following
%%% 1-Caohigh  2-Naolow  3-Naohigh  4-2Hz  5-Spon 

%%% output matrix, metrics of human adult AP and CaT from population simualtions under 1 Hz pacing
Yblock = Metrics_human ;

%%% input matrix, metrics of iPSC-CM AP and CaT from population simulations under 5 different conditions
%%% in the above mentioned order, vertically concatenated together
Xblock = [Metrics_Caohigh,Metrics_Naolow,Metrics_Naohigh,Metrics_2Hz,Metrics_Spon] ;
   
[n_cells,n_parameters] = size(Xblock) ;
[n_cells,n_outputs] = size(Yblock) ;


%%% for both input matrix and output matrix, remove the models with
%%% metrics values differ from the corresponding population mean by more than 3 standard deviations
for i=1:n_outputs
    mu = mean(Yblock(:,i)) ;
    dev = std(Yblock(:,i)) ;
    goodcells = find( Yblock(:,i) > (mu-3*dev) & Yblock(:,i) < (mu+3*dev) ) ;
    goodcount = length(goodcells) ;
    Xblock = Xblock(goodcells,:) ;
    Yblock = Yblock(goodcells,:) ;
end

[n_cells,n_parameters] = size(Xblock) ;
[n_cells,n_outputs] = size(Yblock) ;
    
for i=1:n_parameters
    mu = mean(Xblock(:,i)) ;
    dev = std(Xblock(:,i)) ;
    goodcells = find( Xblock(:,i) > (mu-3*dev) & Xblock(:,i) < (mu+3*dev) ) ;
    goodcount = length(goodcells) ;
    Xblock = Xblock(goodcells,:) ;
    Yblock = Yblock(goodcells,:) ;
end

%%% log transform input matrix and output matrix for better linearity
X = log(abs(Xblock)) ;
Y = log(abs(Yblock)) ;

[n_cells,n_parameters] = size(X) ;
[n_cells,n_outputs] = size(Y) ;
   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Performing the PLS regression

M_X=mean(X);
M_Y=mean(Y);
S_X=std(X);
S_Y=std(Y);

sigmaX = std(X);
sigmaY = std(Y); 

sigma_matrix = [] ;
for i = 1:n_outputs
    sigma_matrix = [sigma_matrix;(ones(1,n_parameters)*sigmaY(i))./sigmaX] ;
end
sigma_matrix = sigma_matrix' ;

%%% number of component to keep for PLS regression
nfact2keep = 40 ; % Adjust according to the column size of input matrix
%%% Supplementary Methods section "Decision on number of components used in PLS regression"


[T,P,W,Wstar,U,b,C,Bpls,Bpls_star,Xhat,Yhat,R2X,R2Y,Fh,SS_Eh,SS_Fh]=...%%%%%%%%%%%%
    PLS_nipals(X,Y,nfact2keep);%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
TSS = sum((Y-ones(n_cells,1)*mean(Y)).^2); % TSS total sum of squares
RSS = sum((Y-Yhat).^2); % RSS resisual sum of squares
R2each = 1 - RSS./TSS;

R2_adj = 1 - (1 - R2each).*((n_cells-1)/(n_cells-n_parameters-1)); % adjusted R^2 
Bpls_scaled = Bpls.*sigma_matrix ;

Bpls_origin = Bpls ;
Bpls_initial = Bpls_origin ;
R2build_initial = R2_adj ;


%% 5-fold cross validation

Yhat_ori = Yhat ;
Yhat_cross = zeros(size(Yhat)) ;

%%% fold of cross validation 
fold_validation = 5 ;

ROUND = round(n_cells/fold_validation) ;
if ROUND > n_cells/fold_validation
   samples_leaveout = ROUND - 2 ;
else
   samples_leaveout = ROUND - 1 ;
end

%%% account for the conditions where n_cells cannot be exactly devided by fold of cross validation
rest = mod(n_cells,fold_validation) ;

M_X=mean(X);
M_Y=mean(Y);
S_X=std(X);
S_Y=std(Y);


for ii=1:fold_validation

  dices_leaveout = (1:samples_leaveout) + (ii-1)*samples_leaveout ;
  
  dices_leavein = setdiff(1:n_cells,dices_leaveout) ;
  X_build = X(dices_leavein,:) ;
  Y_build = Y(dices_leavein,:) ;
  X_test = X(dices_leaveout,:) ;

  X_test_Z = (X_test - ones(samples_leaveout,1)*M_X)./ ...
    (ones(samples_leaveout,1)*S_X) ;

[T,P,W,Wstar,U,b,C,Bpls,Bpls_star,Xhat,Yhat,R2X,R2Y,Fh,SS_Eh,SS_Fh]=...
    PLS_nipals(X_build,Y_build,nfact2keep);
  Yhat_unscale_Z = X_test_Z*Bpls ;
  Yhat_testset = Yhat_unscale_Z.*(ones(samples_leaveout,1)*S_Y) + ...
      ones(samples_leaveout,1)*M_Y ;
  
  Yhat_cross(dices_leaveout,:) = Yhat_testset ;

end  

%%% for the rest few cells
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dices_leaveout = (1:rest+fold_validation*1) + samples_leaveout*fold_validation ;%%%%%%%%%%%%%%
num_test = rest+fold_validation*1 ;

  dices_leavein = setdiff(1:n_cells,dices_leaveout) ;
  X_build = X(dices_leavein,:) ;
  Y_build = Y(dices_leavein,:) ;
  X_test = X(dices_leaveout,:) ;

  X_test_Z = (X_test - ones(num_test,1)*M_X)./ ...
    (ones(num_test,1)*S_X) ;

[T,P,W,Wstar,U,b,C,Bpls,Bpls_star,Xhat,Yhat,R2X,R2Y,Fh,SS_Eh,SS_Fh]=...
    PLS_nipals(X_build,Y_build,nfact2keep);
  Yhat_unscale_Z = X_test_Z*Bpls ;
  Yhat_testset = Yhat_unscale_Z.*(ones(num_test,1)*S_Y) + ...
      ones(num_test,1)*M_Y ;
  
  Yhat_cross(dices_leaveout,:) = Yhat_testset ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

TSS = sum((Y-ones(n_cells,1)*mean(Y)).^2); % TSS total sum of squares
PRESS = sum((Y-Yhat_cross).^2); % RSS resisual sum of squares
R2each_cross_now = 1 - PRESS./TSS;

R2_adj_cross = R2each_cross_now ;
 
Yhat = Yhat_ori ;

%% Scatter plot of cross validation result

figure
handle = gcf ;
set(handle,'Position',[10,10,1500,1000])
set(handle,'PaperPosition',[1 1 12 6]) ;
    
  for i=1:n_outputs
    
    subplot(3,4,i)
    title(outputnames{i})
    hold on
   
    if i == 5 % Vrest usually with negative values
       plot((-exp(Y(:,i))),(-exp(Yhat_cross(:,i))),'bo')
       plotmin = min([(-exp(Y(:,i)));(-exp(Yhat_cross(:,i)))]) ;
       plotmax = max([(-exp(Y(:,i)));(-exp(Yhat_cross(:,i)))]) ;
    else
       plot(exp(Y(:,i)),exp(Yhat_cross(:,i)),'bo')
       plotmin = min([exp(Y(:,i));exp(Yhat_cross(:,i))]) ;
       plotmax = max([exp(Y(:,i));exp(Yhat_cross(:,i))]) ;
    end
    
    xlabel(['True ',outputnames{i}]) 
    ylabel(['Predicted ',outputnames{i}])
    axis([plotmin plotmax plotmin plotmax])
    plot([plotmin,plotmax],[plotmin,plotmax],'k:')
    text((plotmin+(plotmax-plotmin)*0.4),(plotmin+(plotmax-plotmin)*0.2),['R^2_c_r_o_s_s= ',num2str(R2_adj_cross(i))],'FontSize',12)
    set(gca,'TickDir','Out') 
    
  end

