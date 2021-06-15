%% Reset matlab
clear; clc;
%% load data and variables (CCM and FF_Factor)
addpath('D:\OneDrIVe\Freelance\6.CAPM Matlab\');
load('CCM_RFA')
data =xlsread('FF_FACTORS.xls');
MKTRF = data(:,2)/100;      % market excess return
FF3F = data(:,2:4)/100;     % Fama-French factors
RF = data(:,end)/100;       % risk free rate
MKT = MKTRF+RF;             % market return
MKTRF2 = MKTRF.^2           % squared market return
data= array2table(data)
data.Properties.VariableNames([1:5]) = ["mdate", "MktRF", "SMB", "HML", "RF"];
%%
%select only stocks which have full 618 months
[aGroup, permno] = findgroups(CCM.permno);
Count = groupcounts(aGroup);
stock = table(unique(CCM.permno),Count);
stock(stock{:,2}~=618,:) = [];
Data=table()
for i = 1: size(stock,1)
    sec = stock{i,1};
    data= CCM(CCM{:,1}==sec,:);
    Data = [Data;data];
end
%% create variables of selected stocks
Data.mdate           = datetime(Data.year,Data.month,1);
Data.mdate_num_date  = datenum(Data.mdate);
Data = sortrows(Data,{'mdate_num_date','permno'}); 
ret = unstack(Data,'ret','permno', 'GroupingVariables','mdate');
price = unstack(Data,'prc','permno', 'GroupingVariables','mdate');
me_lag = unstack(Data,'me_lag','permno', 'GroupingVariables','mdate'); 
me = unstack(Data,'me','permno', 'GroupingVariables','mdate'); 
me6 = unstack(Data,'me6','permno', 'GroupingVariables','mdate');
me12 = unstack(Data,'me12','permno', 'GroupingVariables','mdate');
be = unstack(Data,'be','permno', 'GroupingVariables','mdate');
exchcd = unstack(Data,'exchcd','permno', 'GroupingVariables','mdate'); 

mdate = table2array(ret(:,1));
ret  = table2array(ret(:,2:end));
price = table2array(price(:,2:end));
me_lag  = table2array(me_lag(:,2:end));
me  = table2array(me(:,2:end));
me6  = table2array(me6(:,2:end));
me12  = table2array(me12(:,2:end));
be  = table2array(be(:,2:end));
exchcd = table2array(exchcd(:,2:end));
ExcessRet= ret-RF

%% Exercise 1:
window = 36;
T = size(MKTRF,1);
%calculate rolling beta and idiosyncratic volatility risk wrt market excess return
for i = 1:T-window
    select = i:i+window-1;
    [B, SEhac, Vhac, RSQ, e] = regressS_HAC(ExcessRet(select,:),[ones(window,1) MKTRF(select,1)],0);
    BetaM(i,:) = B(2,:);
    IV(i,:) = var(e)
end

%%
%calculate rolling beta and idiosyncratic volatility risk wrt squared market excess return
for i = 1:T-window
    select = i:i+window-1;
    [B, SEhac, Vhac, RSQ, e] = regressS_HAC(ExcessRet(select,:),[ones(window,1) MKTRF2(select,1)],0);
    betaV(i,:) = B(2,:);

end

%% drop the first 36 observations from the sample
if size(ExcessRet,1) == 618
    ExcessRet(1:window,:)=[];    
end
if size(mdate,1) == 618
    mdate(1:window,:)=[];   
end
if size(MKTRF,1) == 618
    MKTRF(1:window,:)=[];   
end
if size(MKTRF2,1) == 618
    MKTRF2(1:window,:)=[];   
end
if size(me,1) == 618
    me(1:window,:)=[];   
end
if size(me6,1) == 618
    me6(1:window,:)=[];   
end
if size(be,1) == 618
    be(1:window,:)=[];   
end
if size(ret,1) == 618
    ret(1:window,:)=[];
end
if size(me_lag,1) == 618
    me_lag(1:window,:)=[];
end
if size(RF,1) == 618
    RF(1:window,:)=[];
end
%% Report the periodic statistic BetaV and Signma square
%The betas are estimated overlapping and cannot change too much 
%from month to month, i.e. they are persistent by construction. 
%For that reason, we focus on July Betas and study how their properties 
%change over time:
betaV_jul = betaV(1:12:end,:);
IV_jul = IV(1:12:end,:);
ExcessRet_jul   = ExcessRet(1:12:end,:);
me_jul   = me(1:12:end,:);

% properties of beta in a gIVen year:
for tau = 1:size(betaV_jul,1)
    betaV_jul_tau = betaV_jul(tau,:);
    IV_jul_tau = IV_jul(tau,:);
    ExcessRet_jul_tau   = ExcessRet_jul(tau,:);
    me_jul_tau   = me_jul(tau,:);
    % restrict to the common sample
    isnan1 = isnan(betaV_jul_tau);
    isnan4 = isnan(IV_jul_tau);
    isnan2 = isnan(ExcessRet_jul_tau);
    isnan3 = isnan(me_jul_tau);
    isnanA = isnan1+isnan2+isnan3+isnan4;
    betaV_jul_tau = betaV_jul_tau(isnanA==0);
    IV_jul_tau = IV_jul_tau(isnanA==0);
    ExcessRet_jul_tau = ExcessRet_jul_tau(isnanA==0);
    me_jul_tau = me_jul_tau(isnanA==0);
    
    mu_beta(tau,:)  = mean(betaV_jul_tau);
    sd_beta(tau,:)  = std(betaV_jul_tau'');
    skew_beta(tau,:)  = skewness(betaV_jul_tau');
    [p_beta(tau,:)] = prctile(betaV_jul_tau',[10 50 90]);
    corr_beta(tau,:)  = corr(betaV_jul_tau', ExcessRet_jul_tau');
    
    mu_IV(tau,:)  = mean(IV_jul_tau);
    sd_IV(tau,:)  = std(IV_jul_tau'');
    skew_IV(tau,:)  = skewness(IV_jul_tau');
    [p_IV(tau,:)] = prctile(IV_jul_tau',[10 50 90]);
    corr_IV(tau,:)  = corr(IV_jul_tau', ExcessRet_jul_tau');
    
    mu_me(tau,:)  = mean(me_jul_tau);
    sd_me(tau,:)  = std(me_jul_tau'');
    skew_me(tau,:)  = skewness(me_jul_tau');
    [p_me(tau,:)] = prctile(me_jul_tau',[10 50 90]);
    corr_me(tau,:)  = corr(me_jul_tau', ExcessRet_jul_tau');
end
%Make a table for beta:
tab_BETA = [mu_beta sd_beta skew_beta p_beta corr_beta];
var_txt = {'mu','sd','skew','p10','p50','p90','corr_rx'};
years = [1963+3:2014]
years_txt = num2str(years');
years_txt = cellstr(years_txt)';
TABLE_BETA = array2table(tab_BETA, 'VariableNames',var_txt,'RowNames',years_txt);
TABLE_BETA2 = TABLE_BETA(1:5:end,:) 
disp(TABLE_BETA2)

%Make a table for market equity (me)
tab_ME = [mu_me sd_me skew_me p_me corr_me];
var_txt = {'mu','sd','skew','p10','p50','p90','corr_rx'};
years = [1963+3:2014];
years_txt = num2str(years');
years_txt = cellstr(years_txt)';
TABLE_ME = array2table(tab_ME, 'VariableNames',var_txt,'RowNames',years_txt);
TABLE_ME2 = TABLE_ME(1:5:end,:) 
disp(TABLE_ME2)

%Make a table for IV:
tab_IV = [mu_IV sd_IV skew_IV p_IV corr_IV];
var_txt = {'mu','sd','skew','p10','p50','p90','corr_rx'};
years = [1963+3:2014]
years_txt = num2str(years');
years_txt = cellstr(years_txt)';
TABLE_IV = array2table(tab_IV, 'VariableNames',var_txt,'RowNames',years_txt);
TABLE_IV2 = TABLE_IV(1:5:end,:) 
disp(TABLE_IV2)
%%
%Let us have a look at the persistence of beta over time. How quickly does beta change from year to year?
for tau = 6:size(betaV_jul,1)
    betaV_jul_tau = betaV_jul(tau,:);
    betaV_jul_tau_lag1 = betaV_jul(tau-1,:);
    betaV_jul_tau_lag5 = betaV_jul(tau-5,:);
    isnan1 = isnan(betaV_jul_tau);
    isnan2 = isnan(betaV_jul_tau_lag1);
    isnan3 = isnan(betaV_jul_tau_lag5);
    isnanA = isnan1+isnan2+isnan3;
    betaV_jul_tau = betaV_jul_tau(isnanA==0);
    betaV_jul_tau_lag1 = betaV_jul_tau_lag1(isnanA==0);
    betaV_jul_tau_lag5 = betaV_jul_tau_lag5(isnanA==0);
    
    rho_lag1(tau-5,:) = corr(betaV_jul_tau',betaV_jul_tau_lag1');
    rho_lag5(tau-5,:) = corr(betaV_jul_tau',betaV_jul_tau_lag5');
end

tab_rho_beta = [rho_lag1 rho_lag5];
var_txt = {'rho_beta_lag1','rho__beta_lag5'};
years = [1963+3+5:2014];
years_txt = num2str(years');
years_txt = cellstr(years_txt)';
TABLE_rho_beta = array2table(tab_rho_beta, 'VariableNames',var_txt,'RowNames',years_txt);
TABLE_rho_beta2 = TABLE_rho_beta(1:5:end,:) 
disp(TABLE_rho_beta2)

for tau = 6:size(IV_jul,1)
    IV_jul_tau = IV_jul(tau,:);
    IV_jul_tau_lag1 = IV_jul(tau-1,:);
    IV_jul_tau_lag5 = IV_jul(tau-5,:);
    isnan1 = isnan(IV_jul_tau);
    isnan2 = isnan(IV_jul_tau_lag1);
    isnan3 = isnan(IV_jul_tau_lag5);
    isnanA = isnan1+isnan2+isnan3;
    IV_jul_tau = IV_jul_tau(isnanA==0);
    IV_jul_tau_lag1 = IV_jul_tau_lag1(isnanA==0);
    IV_jul_tau_lag5 = IV_jul_tau_lag5(isnanA==0);
    
    rho_lag1(tau-5,:) = corr(IV_jul_tau',IV_jul_tau_lag1');
    rho_lag5(tau-5,:) = corr(IV_jul_tau',IV_jul_tau_lag5');
end

tab_rho_IV = [rho_lag1 rho_lag5];
var_txt = {'rho_IV_lag1','rho_IV_lag5'};
years = [1963+3+5:2014];
years_txt = num2str(years');
years_txt = cellstr(years_txt)';
TABLE_rho_IV = array2table(tab_rho_IV, 'VariableNames',var_txt,'RowNames',years_txt);
TABLE_rho_IV2 = TABLE_rho_IV(1:5:end,:) 
disp(TABLE_rho_IV2)

%The beta is computed over a 3 year-rolling window, 
%so high persistence to the previous year is still rather mechanical. 
%However, we find that there is not very high persistence 
%at the 4 year horizon when beta estimation no longer overlaps.

    
%% Exercise 2:
%% Test models 1:
% Stacked observations + pooled OLS
ExcessRet_stacked = ExcessRet(:);
beta_stacked = betaV(:);

isnan1 = isnan(ExcessRet_stacked);
isnan2 = isnan(beta_stacked);
isnanA = isnan1+isnan2;
 
Y=ExcessRet_stacked(isnanA==0);
X=[ones(size(Y,1),1) beta_stacked(isnanA==0)];
[Coef, SEhac, Vhac, SEiid, Viid, RSQ] = regress_HAC(Y,X,0);
disp([Coef*100 Coef./SEhac])
disp(RSQ)
% Fama-MacBeth approach:
for t = 1:size(ExcessRet,1)
    ExcessRet_t   = ExcessRet(t,:);
    beta_t = betaV(t,:);
    isnan1 = isnan(ExcessRet_t);
    isnan2 = isnan(beta_t);
    isnanA = isnan1+isnan2;

    Y=ExcessRet_t(isnanA==0)';
    X=[ones(size(Y,1),1) beta_t(isnanA==0)'];
    B = inv(X'*X)*(X'*Y);
    gamma_0_t(t,:) = B(1);
    gamma_1_t(t,:) = B(2);
end

gamma_0 = mean(gamma_0_t);
gamma_1 = mean(gamma_1_t);

se_g0 = std(gamma_0_t)./sqrt(length(gamma_0_t));
se_g1 = std(gamma_1_t)./sqrt(length(gamma_1_t));

disp([[gamma_0; gamma_1]*100 [gamma_0./se_g0; gamma_1./se_g1]])

%% Test models 2:
% Stacked observations + pooled OLS
ExcessRet_stacked = ExcessRet(:);
IV_stacked = IV(:);

isnan1 = isnan(ExcessRet_stacked);
isnan2 = isnan(IV_stacked);
isnanA = isnan1+isnan2;
 
Y=ExcessRet_stacked(isnanA==0);
X=[ones(size(Y,1),1) IV_stacked(isnanA==0)];
[Coef, SEhac, Vhac, SEiid, Viid, RSQ] = regress_HAC(Y,X,0);
disp([Coef*100 Coef./SEhac])
disp(RSQ)
% Fama-MacBeth approach:
for t = 1:size(ExcessRet,1)
    ExcessRet_t   = ExcessRet(t,:);
    IV_t = IV(t,:);
    isnan1 = isnan(ExcessRet_t);
    isnan2 = isnan(IV_t);
    isnanA = isnan1+isnan2;

    Y=ExcessRet_t(isnanA==0)';
    X=[ones(size(Y,1),1) IV_t(isnanA==0)'];
    B = inv(X'*X)*(X'*Y);
    gamma_0_t(t,:) = B(1);
    gamma_1_t(t,:) = B(2);
end

gamma_0 = mean(gamma_0_t);
gamma_1 = mean(gamma_1_t);

se_g0 = std(gamma_0_t)./sqrt(length(gamma_0_t));
se_g1 = std(gamma_1_t)./sqrt(length(gamma_1_t));

disp([[gamma_0; gamma_1]*100 [gamma_0./se_g0; gamma_1./se_g1]])

%% Test models 3:
% Stacked observations + pooled OLS
ExcessRet_stacked = ExcessRet(:);
beta_stacked = betaV(:);
IV_stacked = IV(:);

isnan1 = isnan(ExcessRet_stacked);
isnan2 = isnan(beta_stacked);
isnan3 = isnan(IV_stacked);
isnanA = isnan1+isnan2+isnan3;
 
Y=ExcessRet_stacked(isnanA==0);
X=[ones(size(Y,1),1) beta_stacked(isnanA==0) IV_stacked(isnanA==0)];
[Coef, SEhac, Vhac, SEiid, Viid, RSQ] = regress_HAC(Y,X,0);
disp([Coef*100 Coef./SEhac])
disp(RSQ)
% Fama-MacBeth approach:
for t = 1:size(ExcessRet,1)
    ExcessRet_t   = ExcessRet(t,:);
    beta_t = betaV(t,:);
    IV_t = IV(t,:);
    isnan1 = isnan(ExcessRet_t);
    isnan2 = isnan(beta_t);
    isnan3 = isnan(IV_t);
    isnanA = isnan1+isnan2+isnan3;

    Y=ExcessRet_t(isnanA==0)';
    X=[ones(size(Y,1),1) beta_t(isnanA==0)' IV_t(isnanA==0)'];
    B = inv(X'*X)*(X'*Y);
    gamma_0_t(t,:) = B(1);
    gamma_1_t(t,:) = B(2);
    gamma_2_t(t,:) = B(3);
end

gamma_0 = mean(gamma_0_t);
gamma_1 = mean(gamma_1_t);
gamma_2 = mean(gamma_2_t);
se_g0 = std(gamma_0_t)./sqrt(length(gamma_0_t));
se_g1 = std(gamma_1_t)./sqrt(length(gamma_1_t));
se_g2 = std(gamma_2_t)./sqrt(length(gamma_2_t));
disp([[gamma_0; gamma_1; gamma_2]*100 [gamma_0./se_g0; gamma_1./se_g1; gamma_2./se_g2]])

%% Test models 4:
% Stacked observations + pooled OLS
ExcessRet_stacked = ExcessRet(:);
beta_stacked = betaV(:);
IV_stacked = IV(:);
ln_me = log(me);
ln_beme = log(be.*me);
ln_me_stacked = ln_me(:);
ln_beme_stacked = ln_beme(:);

isnan1 = isnan(ExcessRet_stacked);
isnan2 = isnan(beta_stacked);
isnan3 = isnan(IV_stacked);
isnan4 = isnan(ln_me_stacked);
isnan5 = isnan(ln_beme_stacked);
isnanA = isnan1+isnan2+isnan3+isnan4+isnan5;

Y=ExcessRet_stacked(isnanA==0);
X=[ones(size(Y,1),1) beta_stacked(isnanA==0) IV_stacked(isnanA==0) ln_me(isnanA==0) ln_beme(isnanA==0)];
[Coef, SEhac, Vhac, SEiid, Viid, RSQ] = regress_HAC(Y,X,0);
disp([Coef*100 Coef./SEhac])
disp(RSQ)
% Fama-MacBeth approach:
for t = 1:size(ExcessRet,1)
    ExcessRet_t   = ExcessRet(t,:);
    beta_t = betaV(t,:);
    IV_t = IV(t,:);
    ln_me_t = ln_me(t,:);
    ln_beme_t = ln_beme(t,:);
    isnan1 = isnan(ExcessRet_t);
    isnan2 = isnan(beta_t);
    isnan3 = isnan(IV_t);
    isnan4 = isnan(ln_me_t);
    isnan5 = isnan(ln_beme_t);
    isnanA = isnan1+isnan2+isnan3+isnan4+isnan5;

    Y=ExcessRet_t(isnanA==0)';
    X=[ones(size(Y,1),1) beta_t(isnanA==0)' IV_t(isnanA==0)' ln_me(isnanA==0)' ln_beme(isnanA==0)'];
    B = inv(X'*X)*(X'*Y);
    gamma_0_t(t,:) = B(1);
    gamma_1_t(t,:) = B(2);
    gamma_2_t(t,:) = B(3);
    gamma_3_t(t,:) = B(4);
    gamma_4_t(t,:) = B(5);
end

gamma_0 = mean(gamma_0_t);
gamma_1 = mean(gamma_1_t);
gamma_2 = mean(gamma_2_t);
gamma_3 = mean(gamma_3_t);
gamma_4 = mean(gamma_4_t);
se_g0 = std(gamma_0_t)./sqrt(length(gamma_0_t));
se_g1 = std(gamma_1_t)./sqrt(length(gamma_1_t));
se_g2 = std(gamma_2_t)./sqrt(length(gamma_2_t));
se_g3 = std(gamma_3_t)./sqrt(length(gamma_3_t));
se_g4 = std(gamma_4_t)./sqrt(length(gamma_4_t));
disp([[gamma_0; gamma_1; gamma_2; gamma_3; gamma_4]*100 [gamma_0./se_g0; gamma_1./se_g1; gamma_2./se_g2; gamma_3./se_g3; gamma_4./se_g4]])

%Based on R squared, all 4 models explain very little of the return
%Based on R-squared value, all 4 model have weak explanatory power over the excess return (or return). 
%Market variance βV,I  have very little explanatory while variance of stock IVi can’t predict 13.7% return.
%The increment of variable doesn’t much this explanatory power.

%% Size and value effect
% size
N =size(ExcessRet,2) %we test N assets  
T=size(ExcessRet,1) %We have T months of data
%Under the null, at 5% significance level the critical value for GRS is:
W_crit = chi2inv(1-0.05,N)
%We can convert the statistic to the exact finite sample F-test
GRS = ((T-N-1)/(T*N))*W_crit %W_crit

%the "correct" p-value:
pv_GRS = 1-fcdf(GRS,N,T-N-1)

% value
% assume Sharpe ratio of the market and the alternatIVe Sharpe ratio
mu_m = 0.08/12;
s_m  = 0.20/sqrt(12);
mu_q = 0.085/12;
s_q  = 0.16/sqrt(12);

SRm = mu_m./s_m
SRq = mu_q./s_q
delta= T * (SRq^2-SRm^2)/(1+SRq^2);

GRS_crit = finv(1-0.05,5,T-N-1);

x = (0.01:0.1:10.01)';
p = fcdf(x,5,240-5-1);          % distribution if SRm=SRq
p1 = ncfcdf(x,5,T-N-1,delta); % distribution if SRm not SRq
figure; plot(x,p,'-',x,p1,'-')
legend('F distribution','Noncentral F')

fcdf(GRS_crit,5,240-5-1)
power = 1- ncfcdf(GRS_crit,5,T-N-1,delta) 
% We find a power of about 41.35%, which is not too big.


%% Exercise 3:
% for equally weighted return
Portfolio= pofo_sort_simple(ret,betaV,5,me_lag)
IDX= Portfolio.IDX
ExRet_ew = Portfolio.POFO_ew - RF
ExRet_vw  = Portfolio.POFO_vw - RF

for i = 1:T-window
    select = i:i+window-1;
    [B, SEhac, Vhac, RSQ, e] = regressS_HAC(ExRet_ew(select,:),[ones(window,1) MKTRF(select,1)],0);
    IV_ew(i,:) = var(e)
end
for i = 1:T-window
    select = i:i+window-1;
    [B, SEhac, Vhac, RSQ, e] = regressS_HAC(ExRet_ew(select,:),[ones(window,1) MKTRF2(select,1)],0);
    betaV_ew(i,:) = B(2,:);
end

% drop the first 36 observations from the sample
if size(ExRet_ew,1) == 618-window
    ExRet_ew(1:window,:)=[];    
end
if size(ExRet_vw,1) == 618-window
    ExRet_vw(1:window,:)=[];    
end

%% Test models 1:
% Stacked observations + pooled OLS
ExcessRet_stacked = ExRet_ew(:);
beta_stacked = betaV_ew(:);

isnan1 = isnan(ExcessRet_stacked);
isnan2 = isnan(beta_stacked);
isnanA = isnan1+isnan2;
 
Y=ExcessRet_stacked(isnanA==0);
X=[ones(size(Y,1),1) beta_stacked(isnanA==0)];
[Coef, SEhac, Vhac, SEiid, Viid, RSQ] = regress_HAC(Y,X,0);
disp([Coef*100 Coef./SEhac])
disp(RSQ)
% Fama-MacBeth approach:
for t = 1:size(ExRet_ew,1)
    ExcessRet_t   = ExRet_ew(t,:);
    beta_t = betaV_ew(t,:);
    isnan1 = isnan(ExcessRet_t);
    isnan2 = isnan(beta_t);
    isnanA = isnan1+isnan2;

    Y=ExcessRet_t(isnanA==0)';
    X=[ones(size(Y,1),1) beta_t(isnanA==0)'];
    B = inv(X'*X)*(X'*Y);
    gamma_0_t(t,:) = B(1);
    gamma_1_t(t,:) = B(2);
end

gamma_0 = mean(gamma_0_t);
gamma_1 = mean(gamma_1_t);

se_g0 = std(gamma_0_t)./sqrt(length(gamma_0_t));
se_g1 = std(gamma_1_t)./sqrt(length(gamma_1_t));

disp([[gamma_0; gamma_1]*100 [gamma_0./se_g0; gamma_1./se_g1]])

%% Test models 2:
% Stacked observations + pooled OLS
ExcessRet_stacked = ExRet_ew(:);
IV_stacked = IV_ew(:);

isnan1 = isnan(ExcessRet_stacked);
isnan2 = isnan(IV_stacked);
isnanA = isnan1+isnan2;
 
Y=ExcessRet_stacked(isnanA==0);
X=[ones(size(Y,1),1) IV_stacked(isnanA==0)];
[Coef, SEhac, Vhac, SEiid, Viid, RSQ] = regress_HAC(Y,X,0);
disp([Coef*100 Coef./SEhac])
disp(RSQ)
% Fama-MacBeth approach:
for t = 1:size(ExRet_ew,1)
    ExcessRet_t   = ExRet_ew(t,:);
    IV_t = IV_ew(t,:);
    isnan1 = isnan(ExcessRet_t);
    isnan2 = isnan(IV_t);
    isnanA = isnan1+isnan2;

    Y=ExcessRet_t(isnanA==0)';
    X=[ones(size(Y,1),1) IV_t(isnanA==0)'];
    B = inv(X'*X)*(X'*Y);
    gamma_0_t(t,:) = B(1);
    gamma_1_t(t,:) = B(2);
end

gamma_0 = mean(gamma_0_t);
gamma_1 = mean(gamma_1_t);

se_g0 = std(gamma_0_t)./sqrt(length(gamma_0_t));
se_g1 = std(gamma_1_t)./sqrt(length(gamma_1_t));

disp([[gamma_0; gamma_1]*100 [gamma_0./se_g0; gamma_1./se_g1]])

%% Test models 3:
% Stacked observations + pooled OLS
ExcessRet_stacked = ExRet_ew(:);
beta_stacked = betaV_ew(:);
IV_stacked = IV_ew(:);

isnan1 = isnan(ExcessRet_stacked);
isnan2 = isnan(beta_stacked);
isnan3 = isnan(IV_stacked);
isnanA = isnan1+isnan2+isnan3;
 
Y=ExcessRet_stacked(isnanA==0);
X=[ones(size(Y,1),1) beta_stacked(isnanA==0) IV_stacked(isnanA==0)];
[Coef, SEhac, Vhac, SEiid, Viid, RSQ] = regress_HAC(Y,X,0);
disp([Coef*100 Coef./SEhac])
disp(RSQ)
% Fama-MacBeth approach:
for t = 1:size(ExRet_ew,1)
    ExcessRet_t   = ExRet_ew(t,:);
    beta_t = betaV_ew(t,:);
    IV_t = IV_ew(t,:);
    isnan1 = isnan(ExcessRet_t);
    isnan2 = isnan(beta_t);
    isnan3 = isnan(IV_t);
    isnanA = isnan1+isnan2+isnan3;

    Y=ExcessRet_t(isnanA==0)';
    X=[ones(size(Y,1),1) beta_t(isnanA==0)' IV_t(isnanA==0)'];
    B = inv(X'*X)*(X'*Y);
    gamma_0_t(t,:) = B(1);
    gamma_1_t(t,:) = B(2);
    gamma_2_t(t,:) = B(3);
end

gamma_0 = mean(gamma_0_t);
gamma_1 = mean(gamma_1_t);
gamma_2 = mean(gamma_2_t);
se_g0 = std(gamma_0_t)./sqrt(length(gamma_0_t));
se_g1 = std(gamma_1_t)./sqrt(length(gamma_1_t));
se_g2 = std(gamma_2_t)./sqrt(length(gamma_2_t));
disp([[gamma_0; gamma_1; gamma_2]*100 [gamma_0./se_g0; gamma_1./se_g1; gamma_2./se_g2]])


%% for value weighted return
Portfolio= pofo_sort_simple(ret,betaV,5,me_lag)
IDX= Portfolio.IDX
ExRet_ew = Portfolio.POFO_ew - RF
ExRet_vw  = Portfolio.POFO_vw - RF
for i = 1:T-window
    select = i:i+window-1;
    [B, SEhac, Vhac, RSQ, e] = regressS_HAC(ExRet_vw(select,:),[ones(window,1) MKTRF(select,1)],0);
    IV_vw(i,:) = var(e)
end
for i = 1:T-window
    select = i:i+window-1;
    [B, SEhac, Vhac, RSQ, e] = regressS_HAC(ExRet_vw(select,:),[ones(window,1) MKTRF2(select,1)],0);
    betaV_vw(i,:) = B(2,:);
end

% drop the first 36 observations from the sample
if size(ExRet_vw,1) == 618-window
    ExRet_vw(1:window,:)=[];    
end
if size(ExRet_ew,1) == 618-window
    ExRet_ew(1:window,:)=[];    
end

%% Test models 1:
% Stacked observations + pooled OLS
ExcessRet_stacked = ExRet_vw(:);
beta_stacked = betaV_vw(:);

isnan1 = isnan(ExcessRet_stacked);
isnan2 = isnan(beta_stacked);
isnanA = isnan1+isnan2;
 
Y=ExcessRet_stacked(isnanA==0);
X=[ones(size(Y,1),1) beta_stacked(isnanA==0)];
[Coef, SEhac, Vhac, SEiid, Viid, RSQ] = regress_HAC(Y,X,0);
disp([Coef*100 Coef./SEhac])
disp(RSQ)
% Fama-MacBeth approach:
for t = 1:size(ExRet_vw,1)
    ExcessRet_t   = ExRet_vw(t,:);
    beta_t = betaV_vw(t,:);
    isnan1 = isnan(ExcessRet_t);
    isnan2 = isnan(beta_t);
    isnanA = isnan1+isnan2;

    Y=ExcessRet_t(isnanA==0)';
    X=[ones(size(Y,1),1) beta_t(isnanA==0)'];
    B = inv(X'*X)*(X'*Y);
    gamma_0_t(t,:) = B(1);
    gamma_1_t(t,:) = B(2);
end

gamma_0 = mean(gamma_0_t);
gamma_1 = mean(gamma_1_t);

se_g0 = std(gamma_0_t)./sqrt(length(gamma_0_t));
se_g1 = std(gamma_1_t)./sqrt(length(gamma_1_t));

disp([[gamma_0; gamma_1]*100 [gamma_0./se_g0; gamma_1./se_g1]])

%% Test models 2:
% Stacked observations + pooled OLS
ExcessRet_stacked = ExRet_vw(:);
IV_stacked = IV_vw(:);

isnan1 = isnan(ExcessRet_stacked);
isnan2 = isnan(IV_stacked);
isnanA = isnan1+isnan2;
 
Y=ExcessRet_stacked(isnanA==0);
X=[ones(size(Y,1),1) IV_stacked(isnanA==0)];
[Coef, SEhac, Vhac, SEiid, Viid, RSQ] = regress_HAC(Y,X,0);
disp([Coef*100 Coef./SEhac])
disp(RSQ)
% Fama-MacBeth approach:
for t = 1:size(ExRet_vw,1)
    ExcessRet_t   = ExRet_vw(t,:);
    IV_t = IV_vw(t,:);
    isnan1 = isnan(ExcessRet_t);
    isnan2 = isnan(IV_t);
    isnanA = isnan1+isnan2;

    Y=ExcessRet_t(isnanA==0)';
    X=[ones(size(Y,1),1) IV_t(isnanA==0)'];
    B = inv(X'*X)*(X'*Y);
    gamma_0_t(t,:) = B(1);
    gamma_1_t(t,:) = B(2);
end

gamma_0 = mean(gamma_0_t);
gamma_1 = mean(gamma_1_t);

se_g0 = std(gamma_0_t)./sqrt(length(gamma_0_t));
se_g1 = std(gamma_1_t)./sqrt(length(gamma_1_t));

disp([[gamma_0; gamma_1]*100 [gamma_0./se_g0; gamma_1./se_g1]])

%% Test models 3:
% Stacked observations + pooled OLS
ExcessRet_stacked = ExRet_vw(:);
beta_stacked = betaV_vw(:);
IV_stacked = IV_vw(:);

isnan1 = isnan(ExcessRet_stacked);
isnan2 = isnan(beta_stacked);
isnan3 = isnan(IV_stacked);
isnanA = isnan1+isnan2+isnan3;
 
Y=ExcessRet_stacked(isnanA==0);
X=[ones(size(Y,1),1) beta_stacked(isnanA==0) IV_stacked(isnanA==0)];
[Coef, SEhac, Vhac, SEiid, Viid, RSQ] = regress_HAC(Y,X,0);
disp([Coef*100 Coef./SEhac])
disp(RSQ)
% Fama-MacBeth approach:
for t = 1:size(ExRet_vw,1)
    ExcessRet_t   = ExRet_vw(t,:);
    beta_t = betaV_vw(t,:);
    IV_t = IV_vw(t,:);
    isnan1 = isnan(ExcessRet_t);
    isnan2 = isnan(beta_t);
    isnan3 = isnan(IV_t);
    isnanA = isnan1+isnan2+isnan3;

    Y=ExcessRet_t(isnanA==0)';
    X=[ones(size(Y,1),1) beta_t(isnanA==0)' IV_t(isnanA==0)'];
    B = inv(X'*X)*(X'*Y);
    gamma_0_t(t,:) = B(1);
    gamma_1_t(t,:) = B(2);
    gamma_2_t(t,:) = B(3);
end

gamma_0 = mean(gamma_0_t);
gamma_1 = mean(gamma_1_t);
gamma_2 = mean(gamma_2_t);
se_g0 = std(gamma_0_t)./sqrt(length(gamma_0_t));
se_g1 = std(gamma_1_t)./sqrt(length(gamma_1_t));
se_g2 = std(gamma_2_t)./sqrt(length(gamma_2_t));
disp([[gamma_0; gamma_1; gamma_2]*100 [gamma_0./se_g0; gamma_1./se_g1; gamma_2./se_g2]])


% Using the simple sort portfolio based on βV,I  , 
% we observed the weaker in prediction power of market variance beta and variance of stock. 
%The cross sectional regression show the same result which 
%lower the absolute value of t-statistic of coefficients.
% Based on very low R-squared, we observed that CAPM can’t predict the excess return


%% Fama French 3 factor model

if size(FF3F,1)==618
    FF3F(1:window*2,:)=[]
end

[B, SEhac, Vhac, RSQ, e] = regressS_HAC(ExRet_ew,[ones(size(FF3F,1),1) FF3F(:,1) FF3F(:,2) FF3F(:,3)],0);
disp(B(2:4,:))
disp(RSQ)

[B, SEhac, Vhac, RSQ, e] = regressS_HAC(ExRet_vw,[ones(size(FF3F,1),1) FF3F(:,1) FF3F(:,2) FF3F(:,3)],0);
disp(B(2:4,:))
disp(RSQ)