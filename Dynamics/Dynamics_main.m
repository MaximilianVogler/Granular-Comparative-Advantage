% This code computes the dynamics of the model.

%% Housekeeping
clear all;
close all;
clc;

tstart0 = tic;

addpath('Data');
addpath('Auxiliary Functions');


%% Parameter Vector from Static Model
% load('Estimate_static');
load('GE_Results')

muT=bestParams(1);
sigmaT=bestParams(2);
tau=bestParams(3);
theta=bestParams(4);
F=bestParams(5);

sigma=5;
% theta=kappa*(sigma-1);
% f=f*4.93*.43e-5;
% F=f/sigma;

% New parameters for dynamic model
nu = 0.0542;
mu = -theta*nu^2/2;
rho = 0.2;

Y0=123;
YF0=2*Y0;

%% Setup of model

% Initialize random number generator and set seed
aseed = 1;
rng(aseed);

% Number of shadow firms
MF_large=20000;
MF_small=2800;
MH_large=10000;
MH_small=1400;

% Wages
wR = 1.13;      
w = 1;          
wF = 1/wR; 

% Markups
vMU = 1;

% Competition Structure
BER = 1;

% Scale of model (i.e. number of sectors)
scale = 1;

% Compute number of sectors 
cdshares_init = csvread('cdshares_v3.csv');             % Cobb-Douglas shares from external data source.

S_multiple = 4*scale;                                         
S_init = length(cdshares_init);                         % Number of sectors in data
S = S_init*S_multiple;                                  % Number of sectors used (see footnote 56)

% Assign CD-shares across sectors 
ALPHA = cdshares_init;      

for iloop = 1:S_multiple-1
    ALPHA = [ALPHA;cdshares_init(randperm(S_init))];
end
ALPHA = ALPHA/S_multiple;

% Split sectors into large and small based on CD-share.
split_param = 1.25;
small = (ALPHA<split_param/S);                            
Nsmall = sum(small);                           
Nlarge = S-Nsmall;

% Make random draws for loop
rng(aseed);                                               % Reset random number generator for consistency with old code
rtdraws = randn(1,S);

UH0S = exprnd(1,MH_small,Nsmall);                         % Draw U of most productive small home shadow firm and spacings in each sector from exponential with mean 1
UF0S = exprnd(1,MF_small,Nsmall);                        
UHS = cumsum(UH0S);                                       % Cumulate to get U of each home firm in each sector (see footnote 57)
UFS = cumsum(UF0S);

UH0L = exprnd(1,MH_large,Nlarge);
UF0L = exprnd(1,MF_large,Nlarge);
UHL = cumsum(UH0L);
UFL = cumsum(UF0L);

% Labor Normalization
L0 = 100;

% Given mu_T and sigma_T draw sectoral productivity T_z for each sector z 
RT=exp(muT+sigmaT*rtdraws);
RTS=RT(small);
RTL=RT(~small);

% Draw productivities phi (step 2 of estimation procedure)
ZHS=(UHS./(repmat(RTS,MH_small,1))).^(-1/theta);
ZFS=UFS.^(-1/theta);
ZHL=(UHL./(repmat(RTL,MH_large,1))).^(-1/theta);                        
ZFL=UFL.^(-1/theta);


%% Solve for lower productivity threshold 

% % Compute PE given initial productivities
% [~,Y,YF,LF,~,~,~,~,~,~,LAMBDAFVEC_0,varphi_BB]=GEreplication_vectorized(sigma,theta,F,tau,ALPHA,RTS,RTL,ZHS,ZFS,ZHL,ZFL,w,wF,L0,Y0,YF0,small,vMU,BER,0);

[~,~,~,~,~,~,~,~,~,~,~,~,~,~,LAMBDAFVEC_0,PHIFVEC_0,~,varphi_BB,~,~]=PEreplication_vectorized(sigma,theta,F,tau,ALPHA,RTS,RTL,ZHS,ZFS,ZHL,ZFL,w,wF,Y,YF,small,vMU,BER,0,0);
disp('Loop 0 is finished')
% save('Data/GE_Results','bestParams','Y','YF','LF','varphi_BB')

% Lower bound
ZS = [ZHS;ZFS];
ZL = [ZHL;ZFL];
varphi_low = max([min(ZS),min(ZL)]);

% Upper bound
ZHS_sorted = sort(ZHS,1,'descend');
ZFS_sorted = sort(ZFS,1,'descend');
ZHL_sorted = sort(ZHL,1,'descend');
ZFL_sorted = sort(ZFL,1,'descend');

varphi = zeros(4,S);
RT_sorted = zeros(1,S);
RT_sorted(1,1:Nsmall) = RT(1,small);
RT_sorted(1,Nsmall+1:end) = RT(1,~small);
A = varphi_BB(:,small);
A(A==0) = 1;
B = varphi_BB(:,~small);
B(B==0) = 1;

for i = 1:Nsmall
    varphi(1,i) = ZHS_sorted(A(1,i),i);
    varphi(2,i) = ZHS_sorted(A(2,i),i);
    varphi(3,i) = ZFS_sorted(A(3,i),i);
    varphi(4,i) = ZFS_sorted(A(4,i),i);
end

for i = 1:Nlarge
    varphi(1,Nsmall+i) = ZHL_sorted(B(1,i),i);
    varphi(2,Nsmall+i) = ZHL_sorted(B(2,i),i);
    varphi(3,Nsmall+i) = ZFL_sorted(B(3,i),i);
    varphi(4,Nsmall+i) = ZFL_sorted(B(4,i),i);
end

% Find productivity threshold

% Small Sectors Home
[value_HS,index_HS] = min(abs(log(min(varphi(1:2,1:Nsmall)))-log(min(ZHS))));
[~,zero_index_HS] = min(abs(log(RT(small))));
gross_intercept = log(min(ZHS(:,zero_index_HS)));
% net_intercept = gross_intercept + 0.5*(abs(log(min(varphi(1:2,index_HS)))-log(min(ZHS(:,index_HS)))));
net_intercept_HS = gross_intercept + 0.5*value_HS;
varphi_bar_HS = (1/theta)*log(RT(small))+net_intercept_HS;
% AAA = log(min(varphi(1:2,index_HS)))*ones(1,Nsmall);
% BBB = log(min(ZHS(:,index_HS)))*ones(1,Nsmall);

figure(1);
hold on
scatter(log(RT_sorted(1:Nsmall)),log(min(ZHS)));
scatter(log(RT_sorted(1:Nsmall)),log(min(varphi(1:2,1:Nsmall))));
scatter(log(RT_sorted(1:Nsmall)),varphi_bar_HS);
% % scatter(log(RT_sorted(1:Nsmall)),(1/theta)*log(RT_sorted(1:Nsmall)))
% scatter(log(RT_sorted(1:Nsmall)),AAA);
% scatter(log(RT_sorted(1:Nsmall)),BBB);
lgd = legend('Log of minimum $\varphi_z$','Log of minimum active $\varphi_z$','$\underline{\varphi}_z$');
set(lgd,'location','northwest','interpreter','latex','fontsize',10)
xlabel('$\log(T_z)$','interpreter','latex')
title('Small sectors Home')
saveas(gcf,'Results/Varphi_bar/small_home.png')


% Large Sectors Home
[value_HL,~] = min(abs(log(min(varphi(1:2,Nsmall+1:end)))-log(min(ZHL))));
[~,zero_index_HL] = min(abs(log(RT(~small))));
gross_intercept = log(min(ZHL(:,zero_index_HL)));
net_intercept_HL = gross_intercept + 0.5*value_HL;
varphi_bar_HL = (1/theta)*log(RT(~small))+net_intercept_HL;

figure(2);
hold on
scatter(log(RT_sorted(Nsmall+1:end)),log(min(ZHL)));
scatter(log(RT_sorted(Nsmall+1:end)),log(min(varphi(1:2,Nsmall+1:end))));
scatter(log(RT_sorted(Nsmall+1:end)),varphi_bar_HL);
lgd = legend('Log of minimum $\varphi_z$','Log of minimum active $\varphi_z$','$\underline{\varphi}_z$');
set(lgd,'location','northwest','interpreter','latex','fontsize',10)
xlabel('$\log(T_z)$','interpreter','latex')
title('Large sectors Home')
saveas(gcf,'Results/Varphi_bar/large_home.png')


% Small Sectors Foreign
[value_FS,~] = min(abs(log(min(varphi(3:4,1:Nsmall)))-log(min(ZFS))));
[~,zero_index_FS] = min(abs(log(RT(small))));
gross_intercept = log(min(ZFS(:,zero_index_FS)));
net_intercept_FS = gross_intercept + 0.5*value_FS;
varphi_bar_FS = net_intercept_FS*ones(1,Nsmall);

figure(3);
hold on
scatter(log(RT_sorted(1:Nsmall)),log(min(ZFS)));
scatter(log(RT_sorted(1:Nsmall)),log(min(varphi(3:4,1:Nsmall))));
scatter(log(RT_sorted(1:Nsmall)),varphi_bar_FS);
lgd = legend('Log of minimum $\varphi_z$','Log of minimum active $\varphi_z$','$\underline{\varphi}_z$');
set(lgd,'location','northwest','interpreter','latex','fontsize',10)
xlabel('$\log(T_z)$','interpreter','latex')
title('Small sectors Foreign')
saveas(gcf,'Results/Varphi_bar/small_foreign.png')


% Large Sectors Foreign
[value_FL,~] = min(abs(log(min(varphi(3:4,Nsmall+1:end)))-log(min(ZFL))));
[~,zero_index_FL] = min(abs(log(RT(~small))));
gross_intercept = log(min(ZFL(:,zero_index_FL)));
net_intercept_FL = gross_intercept + 0.5*value_FL;
varphi_bar_FL = net_intercept_FL*ones(1,Nlarge);

figure(4);
hold on
scatter(log(RT_sorted(Nsmall+1:end)),log(min(ZFL)));
scatter(log(RT_sorted(Nsmall+1:end)),log(min(varphi(3:4,Nsmall+1:end))));
scatter(log(RT_sorted(Nsmall+1:end)),varphi_bar_FL);
lgd = legend('Log of minimum $\varphi_z$','Log of minimum active $\varphi_z$','$\underline{\varphi}_z$');
set(lgd,'location','northwest','interpreter','latex','fontsize',10)
xlabel('$\log(T_z)$','interpreter','latex')
title('Large sectors Foreign')
saveas(gcf,'Results/Varphi_bar/large_foreign.png')

% Compute \underline\varphi for each sector for home and foreign, these are
% the relevant vectors that we will use from now on.
varphi_bar_H = zeros(1,S);
varphi_bar_F = zeros(1,S);

varphi_bar_H(1,small) = log(RTS)*(1/theta)+net_intercept_HS;
varphi_bar_H(1,~small) = log(RTL)*(1/theta)+net_intercept_HL;

varphi_bar_F(1,small) = net_intercept_FS;
varphi_bar_F(1,~small) = net_intercept_FL;

%% Simulate Dynamics

% Set all productivities below threshold to zero
inactive_HS = log(ZHS) < varphi_bar_H(1,small);
inactive_HL = log(ZHL) < varphi_bar_H(1,~small);
inactive_FS = log(ZFS) < varphi_bar_F(1,small);
inactive_FL = log(ZFL) < varphi_bar_F(1,~small);

ZHS(inactive_HS) = 1e-10;
ZHL(inactive_HL) = 1e-10;
ZFS(inactive_FS) = 1e-10;
ZFL(inactive_FL) = 1e-10;

% Cut off lower part of matrix (rows that consist only of zeros)
cutoff_HS = min(sum(inactive_HS==1));
cutoff_HL = min(sum(inactive_HL==1));
cutoff_FS = min(sum(inactive_FS==1));
cutoff_FL = min(sum(inactive_FL==1));

ZHS = ZHS(1:end-cutoff_HS,:);
ZHL = ZHL(1:end-cutoff_HL,:);
ZFS = ZFS(1:end-cutoff_FS,:);
ZFL = ZFL(1:end-cutoff_FL,:);

[ZHS_length,~] = size(ZHS); 
[ZHL_length,~] = size(ZHL); 
[ZFS_length,~] = size(ZFS); 
[ZFL_length,~] = size(ZFL); 

inactive_HS = log(ZHS) < varphi_bar_H(1,small);
inactive_HL = log(ZHL) < varphi_bar_H(1,~small);
inactive_FS = log(ZFS) < varphi_bar_F(1,small);
inactive_FL = log(ZFL) < varphi_bar_F(1,~small);

VB_HS = repmat(varphi_bar_H(1,small),ZHS_length,1);
VB_HL = repmat(varphi_bar_H(1,~small),ZHL_length,1);
VB_FS = repmat(varphi_bar_F(1,small),ZFS_length,1);
VB_FL = repmat(varphi_bar_F(1,~small),ZFL_length,1);

% Length of dynamic simulation

% ZHS_Save = [];
rng(aseed);

% Years that are actually being recorded
RECORD = [1:11];
R_length = length(RECORD);
T = RECORD(end);

LAMBDAFVEC_t = zeros(R_length,S);
DVEC_t = zeros(R_length,S);
XVEC_t = zeros(R_length,S);
PHIFVEC_t = zeros(R_length,S);


% Determine productivity draws
su = sqrt(rho*nu^2);
sv = sqrt(nu^2*(1-rho));

counter = 1;

for t=1:T

    uz_HS = repmat(randn(1,S),ZHS_length,1);
    uz_HL = repmat(randn(1,S),ZHL_length,1);
    uz_FS = repmat(randn(1,S),ZFS_length,1);
    uz_FL = repmat(randn(1,S),ZFL_length,1);

    vz_HS = randn(ZHS_length,S);
    vz_HL = randn(ZHL_length,S);
    vz_FS = randn(ZFS_length,S);
    vz_FL = randn(ZFL_length,S);
    
    eps_HS = su*uz_HS+sv*vz_HS;
    eps_HL = su*uz_HL+sv*vz_HL;
    eps_FS = su*uz_FS+sv*vz_FS;
    eps_FL = su*uz_FL+sv*vz_FL;
    
    ZHS(~inactive_HS) = exp(VB_HS(~inactive_HS)+abs(log(ZHS(~inactive_HS))-VB_HS(~inactive_HS)+mu+eps_HS(~inactive_HS)));
%     ZHS_Save =[ZHS_Save,ZHS(:)];
    ZHL(~inactive_HL) = exp(VB_HL(~inactive_HL)+abs(log(ZHL(~inactive_HL))-VB_HL(~inactive_HL)+mu+eps_HL(~inactive_HL)));
    ZFS(~inactive_FS) = exp(VB_FS(~inactive_FS)+abs(log(ZFS(~inactive_FS))-VB_FS(~inactive_FS)+mu+eps_FS(~inactive_FS)));
    ZFL(~inactive_FL) = exp(VB_FL(~inactive_FL)+abs(log(ZFL(~inactive_FL))-VB_FL(~inactive_FL)+mu+eps_FL(~inactive_FL)));
   
%     ZHS(~inactive_HS) = exp(VB_HS(~inactive_HS)+abs(log(ZHS(~inactive_HS))-VB_HS(~inactive_HS)+mu+nu*randn(length(ZHS(~inactive_HS)),1)));
% %     ZHS_Save =[ZHS_Save,ZHS(:)];
%     ZHL(~inactive_HL) = exp(VB_HL(~inactive_HL)+abs(log(ZHL(~inactive_HL))-VB_HL(~inactive_HL)+mu+nu*randn(length(ZHL(~inactive_HL)),1)));
%     ZFS(~inactive_FS) = exp(VB_FS(~inactive_FS)+abs(log(ZFS(~inactive_FS))-VB_FS(~inactive_FS)+mu+nu*randn(length(ZFS(~inactive_FS)),1)));
%     ZFL(~inactive_FL) = exp(VB_FL(~inactive_FL)+abs(log(ZFL(~inactive_FL))-VB_FL(~inactive_FL)+mu+nu*randn(length(ZFL(~inactive_FL)),1)));
%     
 
    % If the year is part of RECORD, record PE results
    if any(RECORD==t)
        [~,~,~,~,~,~,~,~,~,~,~,~,~,~,LAMBDAFVEC,PHIFVEC,~,~,XVEC,DVEC,DSHM_small,DSHM_large,DSHMS_small,DSHMS_large]=PEreplication_vectorized(sigma,theta,F,tau,ALPHA,RTS,RTL,ZHS,ZFS,ZHL,ZFL,w,wF,Y,YF,small,vMU,BER,0,0);
        LAMBDAFVEC_t(t,:) = LAMBDAFVEC;
        PHIFVEC_t(t,:) = PHIFVEC;
        XVEC_t(t,:) = XVEC; 
        DVEC_t(t,:) = DVEC;
        save_share_small{t} = DSHM_small;
        save_share_large{t} = DSHM_large;
        disp(['Loop ',num2str(t),' is finished']);
    end
    
end

% %% MEAN REVERSION
% 
% % GAMMAF
% GAMMAFVEC_t = LAMBDAFVEC_t-PHIFVEC_t;
% 
% % Delta LAMBDAF
% DELTA_LAMBDAFVEC_t = LAMBDAFVEC_t-repmat(LAMBDAFVEC_t(1,:),R_length+1,1);
% 
% % Prepare plotting variables
% ptiles = 10;
% LPCT = prctile(GAMMAFVEC_t(1,:),[1/ptiles:1/ptiles:1-1/ptiles]*100)';
% 
% mean_20 = zeros(ptiles,1);         
% mean_50 = zeros(ptiles,1);
% 
% min_20 = zeros(ptiles,1);
% min_50 = zeros(ptiles,1);
% 
% max_20 = zeros(ptiles,1);
% max_50 = zeros(ptiles,1);
% 
% DELTA_20 = DELTA_LAMBDAFVEC_t(2,:);
% DELTA_50 = DELTA_LAMBDAFVEC_t(3,:);
% 
% mean_20(1) = mean(DELTA_20(GAMMAFVEC_t(1,:)<LPCT(1)));
% mean_50(1) = mean(DELTA_50(GAMMAFVEC_t(1,:)<LPCT(1)));
% min_20(1) = prctile(DELTA_20(GAMMAFVEC_t(1,:)<LPCT(1)),10);
% min_50(1) = prctile(DELTA_50(GAMMAFVEC_t(1,:)<LPCT(1)),10);
% max_20(1) = prctile(DELTA_20(GAMMAFVEC_t(1,:)<LPCT(1)),90);
% max_50(1) = prctile(DELTA_50(GAMMAFVEC_t(1,:)<LPCT(1)),90);
% for i=2:ptiles-1
%    mean_20(i) = mean(DELTA_20(GAMMAFVEC_t(1,:)<LPCT(i) & GAMMAFVEC_t(1,:)>=LPCT(i-1))); 
%    mean_50(i) = mean(DELTA_50(GAMMAFVEC_t(1,:)<LPCT(i) & GAMMAFVEC_t(1,:)>=LPCT(i-1))); 
%    min_20(i) = prctile(DELTA_20(GAMMAFVEC_t(1,:)<LPCT(i) & GAMMAFVEC_t(1,:)>=LPCT(i-1)),10); 
%    min_50(i) = prctile(DELTA_50(GAMMAFVEC_t(1,:)<LPCT(i) & GAMMAFVEC_t(1,:)>=LPCT(i-1)),10); 
%    max_20(i) = prctile(DELTA_20(GAMMAFVEC_t(1,:)<LPCT(i) & GAMMAFVEC_t(1,:)>=LPCT(i-1)),90); 
%    max_50(i) = prctile(DELTA_50(GAMMAFVEC_t(1,:)<LPCT(i) & GAMMAFVEC_t(1,:)>=LPCT(i-1)),90); 
% end
% mean_20(ptiles) = mean(DELTA_20(GAMMAFVEC_t(1,:)>=LPCT(ptiles-1)));
% mean_50(ptiles) = mean(DELTA_50(GAMMAFVEC_t(1,:)>=LPCT(ptiles-1)));
% min_20(ptiles) = prctile(DELTA_20(GAMMAFVEC_t(1,:)>=LPCT(ptiles-1)),10);
% min_50(ptiles) = prctile(DELTA_50(GAMMAFVEC_t(1,:)>=LPCT(ptiles-1)),10);
% max_20(ptiles) = prctile(DELTA_20(GAMMAFVEC_t(1,:)>=LPCT(ptiles-1)),90);
% max_50(ptiles) = prctile(DELTA_50(GAMMAFVEC_t(1,:)>=LPCT(ptiles-1)),90);
% 
% 
% % Plot mean reversion graph
% figure(5)
% bar(mean_50,'FaceColor','[0, 0.4470, 0.7410]','FaceAlpha',0.33)
% hold on
% bar(mean_20,'FaceColor','[0.900, 0.20, 0.05]')%,'FaceAlpha',0.65)
% lgd = legend('50 years','20 years');
% set(lgd,'box','off','location','northeast','interpreter','latex','fontsize',20)
% % ylim([0, .25])
% % set(gca,'ytick',[0.05:.05:.25])
% xlim([0.5, 10.55])
% ylim([-.08, .065])
% xlabel('Deciles of sectors, by granular $\Gamma_z^\ast$','interpreter','latex','fontsize',20)
% ylabel('Expected change in export share, $\Delta \Lambda_z^\ast$','interpreter','latex','fontsize',20)
% set(gca,'FontSize',16)
% set(gca,'xtick',[0.5:1:9.5 10.55])
% LPCT0 = LPCT;
% LPCT0(5)=0;
% set(gca,'xticklabel',[{ -1} round(LPCT0(1:4)'*100)/100 {  0} {  0} round(LPCT0(7:end)'*100)/100 1])
% saveas(gcf,'Results/mean_reversion.png')
% 
% figure(6)
% hold on
% bar(mean_20,'FaceColor','[0.900, 0.20, 0.05]')%,'FaceAlpha',0.65)
% errorbar(1:10,mean_20,min_20,max_20,'.','LineWidth',2,'Color',[0, 0.4470, 0.7410])
% lgd = legend('20 years');
% set(lgd,'box','off','location','northeast','interpreter','latex','fontsize',20)
% % ylim([0, .25])
% % set(gca,'ytick',[0.05:.05:.25])
% xlim([0.5, 10.55])
% ylim([-.2, .2])
% xlabel('Deciles of sectors, by granular $\Gamma_z^\ast$','interpreter','latex','fontsize',20)
% ylabel('Expected change in export share, $\Delta \Lambda_z^\ast$','interpreter','latex','fontsize',20)
% set(gca,'FontSize',16)
% set(gca,'xtick',[0.5:1:9.5 10.55])
% LPCT0 = LPCT;
% LPCT0(5)=0;
% set(gca,'xticklabel',[{ -1} round(LPCT0(1:4)'*100)/100 {  0} {  0} round(LPCT0(7:end)'*100)/100 1])
% saveas(gcf,'Results/mean_reversion_20.png')
% 
% figure(7)
% hold on
% bar(mean_50,'FaceColor','[0, 0.4470, 0.7410]','FaceAlpha',0.33)
% errorbar(1:10,mean_50,min_50,max_50,'.','LineWidth',2,'Color',[0.8500, 0.3250, 0.0980])
% lgd = legend('50 years');
% set(lgd,'box','off','location','northeast','interpreter','latex','fontsize',20)
% % ylim([0, .25])
% % set(gca,'ytick',[0.05:.05:.25])
% xlim([0.5, 10.55])
% ylim([-.3, .3])
% xlabel('Deciles of sectors, by granular $\Gamma_z^\ast$','interpreter','latex','fontsize',20)
% ylabel('Expected change in export share, $\Delta \Lambda_z^\ast$','interpreter','latex','fontsize',20)
% set(gca,'FontSize',16)
% set(gca,'xtick',[0.5:1:9.5 10.55])
% LPCT0 = LPCT;
% LPCT0(5)=0;
% set(gca,'xticklabel',[{ -1} round(LPCT0(1:4)'*100)/100 {  0} {  0} round(LPCT0(7:end)'*100)/100 1])
% saveas(gcf,'Results/mean_reversion_50.png')


%% Compute Firm Dynamics Moments (Table 4)

keepstd_1 = zeros(S,T-1)*NaN;
keepstd_2 = zeros(S,T-1)*NaN;
keepstd_3 = zeros(S,T-1)*NaN;
keepstd_4 = zeros(S,T-1)*NaN;
keepstd_5 = zeros(S,T-1)*NaN;
keepstd_6 = zeros(S,T-1)*NaN;
keepstd_7 = zeros(S,T-1)*NaN;
keepstd_8 = zeros(S,T-1)*NaN;
keepstd_9 = zeros(S,T-1)*NaN;
keepcorr_10 = zeros(S,1)*NaN;
keepcorr_11 = zeros(S,1)*NaN;
keepcorr_12 = zeros(S,1)*NaN;
keepcorr_13 = zeros(S,1)*NaN;
keepcorr_14 = zeros(S,1)*NaN;
keepcorr_15 = zeros(S,1)*NaN;

% avg_ms = zeros(S,1);
% avg_ms_50 = zeros(S,1);
% avg_ms_80 = zeros(S,1);

% avg_ms = nanmean(avg_ms_int);
% avg_ms_50 = nanmean(avg_ms_50_int);
% avg_ms_80 = nanmean(avg_ms_80_int);
    
% First the small sectors
for z=1:Nsmall
    
    idx_50_small = save_share_small{1}(:,z)>=quantile(save_share_small{1}(save_share_small{1}(:,z)>0,z),0.5);
    idx_80_small = save_share_small{1}(:,z)>=quantile(save_share_small{1}(save_share_small{1}(:,z)>0,z),0.8);
    avg_ms = mean(save_share_small{1}(save_share_small{1}(:,z)>0,z));
    avg_ms_50 = mean(save_share_small{1}(idx_50_small,z));
    avg_ms_80 = mean(save_share_small{1}(idx_80_small,z));
    
    % Standard deviation across firms within time-sector (Moment 1)
    for tt=1:T-1
        a2 = save_share_small{tt+1}(:,z);
        a1 = save_share_small{tt}(:,z);
        idx = (a2>0).*(a1>0);               % Only if firm exists in both periods
        keepstd_1(z,tt) = std(a2(idx>0)-a1(idx>0));
    end
    
    % Standard deviation across firms (top 20%) time-sector (Moment 2)
    for tt=1:T-1
       a2 = save_share_small{tt+1}(idx_80_small,z);
       a1 = save_share_small{t}(idx_80_small,z);
       idx = (a2>0).*(a1>0); 
       keepstd_2(z,tt) = std(a2(idx>0)-a1(idx>0));
    end
    
    % Standard deviation across firms (top 50%) time-sector (Moment 3)
    for tt=1:T-1
       a2 = save_share_small{tt+1}(idx_50_small,z);
       a1 = save_share_small{t}(idx_50_small,z);
       idx = (a2>0).*(a1>0); 
       keepstd_3(z,tt) = std(a2(idx>0)-a1(idx>0));
    end
    
     % Standard deviation of log across firms within time-sector (Moment 4)
    for tt=1:T-1
        a2 = save_share_small{tt+1}(:,z);
        a1 = save_share_small{tt}(:,z);
        idx = (a2>0).*(a1>0);               % Only if firm exists in both periods
        keepstd_4(z,tt) = std(log(a2(idx>0))-log(a1(idx>0)));
    end
    
    % Standard deviation of log across firms (top 20%) time-sector (Moment 5)
    for tt=1:T-1
       a2 = save_share_small{tt+1}(idx_80_small,z);
       a1 = save_share_small{t}(idx_80_small,z);
       idx = (a2>0).*(a1>0); 
       keepstd_5(z,tt) = std(log(a2(idx>0))-log(a1(idx>0)));
    end
    
    % Standard deviation of log across firms (top 50%) time-sector (Moment 6)
    for tt=1:T-1
       a2 = save_share_small{tt+1}(idx_50_small,z);
       a1 = save_share_small{t}(idx_50_small,z);
       idx = (a2>0).*(a1>0); 
       keepstd_6(z,tt) = std(log(a2(idx>0))-log(a1(idx>0)));
    end
 
    % Normalized standard deviation across firms within time-sector (Moment 7)
    for tt=1:T-1
        a2 = save_share_small{tt+1}(:,z);
        a1 = save_share_small{tt}(:,z);
        idx = (a2>0).*(a1>0);               % Only if firm exists in both periods
        keepstd_7(z,tt) = std((a2(idx>0)-a1(idx>0)))/avg_ms;
    end
    
    % Normalized standard deviation across firms (top 20%) time-sector (Moment 8)
    for tt=1:T-1
       a2 = save_share_small{tt+1}(idx_80_small,z);
       a1 = save_share_small{t}(idx_80_small,z);
       idx = (a2>0).*(a1>0); 
       keepstd_8(z,tt) = std(a2(idx>0)-a1(idx>0))/avg_ms_80;
    end
    
    % Normalized standard deviation across firms (top 50%) time-sector (Moment 9)
    for tt=1:T-1
       a2 = save_share_small{tt+1}(idx_50_small,z);
       a1 = save_share_small{t}(idx_50_small,z);
       idx = (a2>0).*(a1>0); 
       keepstd_9(z,tt) = std(a2(idx>0)-a1(idx>0))/avg_ms_50;
    end

    % Correlation over firms (Moment 10)
    a0 = save_share_small{1}(:,z);
    aT = save_share_small{T}(:,z);
    idx = (a0>0).*(aT>0);                   % Again, only if firm exists in both periods
    if sum(idx)>1
        keepcorr_10(z) = corr(aT(idx>0),a0(idx>0));
    end
    
    % Correlation over firms (top 20%) (Moment 11)
    a0 = save_share_small{1}(idx_80_small,z);
    aT = save_share_small{T}(idx_80_small,z);
    idx = (a0>0).*(aT>0);                   % Again, only if firm exists in both periods
    if sum(idx)>1
        keepcorr_11(z) = corr(aT(idx>0),a0(idx>0));
    end
    
    % Correlation over firms (top 50%) (Moment 12)
    a0 = save_share_small{1}(idx_50_small,z);
    aT = save_share_small{T}(idx_50_small,z);
    idx = (a0>0).*(aT>0);                   % Again, only if firm exists in both periods
    if sum(idx)>1
        keepcorr_12(z) = corr(aT(idx>0),a0(idx>0));
    end
    
    % Correlation over firms log (Moment 13)
    a0 = save_share_small{1}(:,z);
    aT = save_share_small{T}(:,z);
    idx = (a0>0).*(aT>0);                   % Again, only if firm exists in both periods
    if sum(idx)>1
        keepcorr_13(z) = corr(log(aT(idx>0)),log(a0(idx>0)));
    end
    
    % Correlation over firms log (top 20%) (Moment 14)
    a0 = save_share_small{1}(idx_80_small,z);
    aT = save_share_small{T}(idx_80_small,z);
    idx = (a0>0).*(aT>0);                   % Again, only if firm exists in both periods
    if sum(idx)>1
        keepcorr_14(z) = corr(log(aT(idx>0)),log(a0(idx>0)));
    end
    
    % Correlation over firms log (top 50%) (Moment 15)
    a0 = save_share_small{1}(idx_50_small,z);
    aT = save_share_small{T}(idx_50_small,z);
    idx = (a0>0).*(aT>0);                   % Again, only if firm exists in both periods
    if sum(idx)>1
        keepcorr_15(z) = corr(log(aT(idx>0)),log(a0(idx>0)));
    end
end

% Now do the same for large sectors

for z=1:Nlarge
    
    idx_50_large = save_share_large{1}(:,z)>=quantile(save_share_large{1}(save_share_large{1}(:,z)>0,z),0.5);
    idx_80_large = save_share_large{1}(:,z)>=quantile(save_share_large{1}(save_share_large{1}(:,z)>0,z),0.8);
    avg_ms = mean(save_share_large{1}(save_share_large{1}(:,z)>0,z));
    avg_ms_50 = mean(save_share_large{1}(idx_50_large,z));
    avg_ms_80 = mean(save_share_large{1}(idx_80_large,z));
    
    % Moment 1
    for tt=1:T-1
        a2 = save_share_large{tt+1}(:,z);
        a1 = save_share_large{tt}(:,z);
        idx = (a2>0).*(a1>0);
        keepstd_1(Nsmall+z,tt) = std(a2(idx>0)-a1(idx>0));
    end
    
    % Moment 2
    for tt=1:T-1
       a2 = save_share_large{tt+1}(idx_80_large,z);
       a1 = save_share_large{t}(idx_80_large,z);
       idx = (a2>0).*(a1>0); 
       keepstd_2(Nsmall+z,tt) = std(a2(idx>0)-a1(idx>0));
    end
    
    % Moment 3
    for tt=1:T-1
       a2 = save_share_large{tt+1}(idx_50_large,z);
       a1 = save_share_large{t}(idx_50_large,z);
       idx = (a2>0).*(a1>0); 
       keepstd_3(Nsmall+z,tt) = std(a2(idx>0)-a1(idx>0));
    end
    
    % Moment 4
    for tt=1:T-1
        a2 = save_share_large{tt+1}(:,z);
        a1 = save_share_large{tt}(:,z);
        idx = (a2>0).*(a1>0);
        keepstd_4(Nsmall+z,tt) = std(log(a2(idx>0))-log(a1(idx>0)));
    end
    
    % Moment 5
    for tt=1:T-1
       a2 = save_share_large{tt+1}(idx_80_large,z);
       a1 = save_share_large{t}(idx_80_large,z);
       idx = (a2>0).*(a1>0); 
       keepstd_5(Nsmall+z,tt) = std(log(a2(idx>0))-log(a1(idx>0)));
    end
    
    % Moment 6
    for tt=1:T-1
       a2 = save_share_large{tt+1}(idx_50_large,z);
       a1 = save_share_large{t}(idx_50_large,z);
       idx = (a2>0).*(a1>0); 
       keepstd_6(Nsmall+z,tt) = std(log(a2(idx>0))-log(a1(idx>0)));
    end
    
    % Moment 7
    for tt=1:T-1
        a2 = save_share_large{tt+1}(:,z);
        a1 = save_share_large{tt}(:,z);
        idx = (a2>0).*(a1>0);
        keepstd_7(Nsmall+z,tt) = std(a2(idx>0)-a1(idx>0))/avg_ms;
    end
    
    % Moment 8
    for tt=1:T-1
       a2 = save_share_large{tt+1}(idx_80_large,z);
       a1 = save_share_large{t}(idx_80_large,z);
       idx = (a2>0).*(a1>0); 
       keepstd_8(Nsmall+z,tt) = std(a2(idx>0)-a1(idx>0))/avg_ms_80;
    end
    
    % Moment 9
    for tt=1:T-1
       a2 = save_share_large{tt+1}(idx_50_large,z);
       a1 = save_share_large{t}(idx_50_large,z);
       idx = (a2>0).*(a1>0); 
       keepstd_9(Nsmall+z,tt) = std(a2(idx>0)-a1(idx>0))/avg_ms_50;
    end
    
    % Moment 10
    a0 = save_share_large{1}(:,z);
    aT = save_share_large{T}(:,z);
    idx = (a0>0).*(aT>0);
    if sum(idx)>1
        keepcorr_10(Nsmall+z) = corr(aT(idx>0),a0(idx>0));
    end
    
    % Moment 11
    a0 = save_share_large{1}(idx_80_large,z);
    aT = save_share_large{T}(idx_80_large,z);
    idx = (a0>0).*(aT>0);
    if sum(idx)>1
        keepcorr_11(Nsmall+z) = corr(aT(idx>0),a0(idx>0));
    end
    
    % Moment 12
    a0 = save_share_large{1}(idx_50_large,z);
    aT = save_share_large{T}(idx_50_large,z);
    idx = (a0>0).*(aT>0);
    if sum(idx)>1
        keepcorr_12(Nsmall+z) = corr(aT(idx>0),a0(idx>0));
    end
    
    % Moment 13
    a0 = save_share_large{1}(:,z);
    aT = save_share_large{T}(:,z);
    idx = (a0>0).*(aT>0);
    if sum(idx)>1
        keepcorr_13(Nsmall+z) = corr(log(aT(idx>0)),log(a0(idx>0)));
    end
    
    % Moment 14
    a0 = save_share_large{1}(idx_80_large,z);
    aT = save_share_large{T}(idx_80_large,z);
    idx = (a0>0).*(aT>0);
    if sum(idx)>1
        keepcorr_14(Nsmall+z) = corr(log(aT(idx>0)),log(a0(idx>0)));
    end
    
    % Moment 15
    a0 = save_share_large{1}(idx_50_large,z);
    aT = save_share_large{T}(idx_50_large,z);
    idx = (a0>0).*(aT>0);
    if sum(idx)>1
        keepcorr_15(Nsmall+z) = corr(log(aT(idx>0)),log(a0(idx>0)));
    end
         
end

mom1 = nanmedian(nanmedian(keepstd_1,2));
mom2 = nanmedian(nanmedian(keepstd_2,2));
mom3 = nanmedian(nanmedian(keepstd_3,2));
mom4 = nanmedian(nanmedian(keepstd_4,2));
mom5 = nanmedian(nanmedian(keepstd_5,2));
mom6 = nanmedian(nanmedian(keepstd_6,2));
mom7 = nanmedian(nanmedian(keepstd_7,2));
mom8 = nanmedian(nanmedian(keepstd_8,2));
mom9 = nanmedian(nanmedian(keepstd_9,2));
mom10 = nanmedian(keepcorr_10);
mom11 = nanmedian(keepcorr_11);
mom12 = nanmedian(keepcorr_12);
mom13 = nanmedian(keepcorr_13);
mom14 = nanmedian(keepcorr_14);
mom15 = nanmedian(keepcorr_15);

DynMom = [mom1,mom2,mom3,mom4,mom5,mom6,mom7,mom8,mom9,mom10,mom11,mom12,mom13,mom14,mom15];

%% Compute Sectoral Dynamics Moments

% Create Time Fixed Effects 
T = zeros((R_length-1)*S,R_length-1);
for t=1:R_length-1
    for z=1:S
        T((z-1)*(R_length-1)+t,t) = 1;
    end
end

% % Create Sectoral Fixed Effects
Z = zeros((R_length-1)*S,S);
for z=1:S
    Z((z-1)*(R_length-1)+1:z*(R_length-1),z) = ones(R_length-1,1);
end

% Get rid of overdue FE (since we have a constant)
T = T(:,1:end-1);
Z = Z(:,1:end-1);
 
% X_t
X_dep = XVEC_t(2:end,:);

%X_{t-1}
X_indep = XVEC_t(1:end-1,:);

% Correcting for zeros
ind = (X_dep>0)&(X_indep>0);
X_dep = log(X_dep(ind));
X_indep = log(X_indep(ind));

T1 = zeros(length(X_dep),R_length-2);
for t=1:R_length-2
    T1(:,t) = T(ind,t);
end
T=T1;

Z1 = zeros(length(X_dep),S-1);
for z=1:S-1
   Z1(:,z)=Z(ind,z); 
end
Z=Z1;

% Regression w/o fixed effects
stats = regstats(X_dep,X_indep,'linear',{'beta','covb'});
smom1 = stats.beta(2);
sd = diag(sqrt(stats.covb));
sd1 = sd(2);

% Regression w/ time fixed effects
stats = regstats(X_dep,[X_indep,T],'linear',{'beta','covb'});
smom2 = stats.beta(2);
sd = diag(sqrt(stats.covb));
sd2 = sd(2);

% Regression w/ sector fixed effects
stats = regstats(X_dep,[X_indep,Z],'linear',{'beta','covb'});
smom3 = stats.beta(2);
sd = diag(sqrt(stats.covb));
sd3 = sd(2);

% Regression w/ time and sector fixed effects
stats = regstats(X_dep,[X_indep,T,Z],'linear',{'beta','covb'});
smom4 = stats.beta(2);
sd = diag(sqrt(stats.covb));
sd4 = sd(2);

% DynMom2 = [smom1,smom2,smom3,smom4];
% DynMom2_sd = [sd1,sd2,sd3,sd4];
%% Generate data for Stata (verification of the above in Stata)
ID = (1:S);
ID = repmat(ID,R_length-1,1);
YEAR = repmat(RECORD(2:end)',S,1);
X_t = XVEC_t(2:end,:);
X_t1 = XVEC_t(1:end-1,:);
DATA = [ID(:),YEAR,X_t(:),X_t1(:)];
fname = sprintf('Results/Data/sectoral_regdata_%d',S);
fname2 = sprintf('Results/Data/sectoral_regdata_%d.csv',S);
save(fname,'DATA');
title1 = {'ID','Year','X_t','X_t1'};
TT = cell2table(num2cell(DATA),'VariableNames',title1);
writetable(TT,fname2);

% Sectoral Standard Deviation

X_t = XVEC_t(2:end,:);
X_t1 = XVEC_t(1:end-1,:);

sectoral_sd = zeros(R_length-1,1);
for t=1:R_length-1
    ind = (X_t(t,:)>0)&(X_t1(t,:)>0);
    sectoral_sd(t) = std(log(X_t(t,ind))-log(X_t1(t,ind)));
end
smom5 = nanmedian(sectoral_sd);

% Sectoral Correlation

ind = (XVEC_t(1,:)>0)&(XVEC_t(end,:)>0);
smom6 = corr(log(XVEC_t(1,ind))',log(XVEC_t(end,ind))');
DynMom2 = [smom1,smom2,smom3,smom4,smom5,smom6];
% DynMom2 = [smom5,smom6];
DynMom2_sd = [sd1,sd2,sd3,sd4];
DynMom2_t = [smom1/sd1,smom2/sd2,smom3/sd3,smom4/sd4];

%% 3rd Class of Moments

K=20;
DynMom3_beta = zeros(K,1);
DynMom3_20_beta = zeros(K,1);
TOP = zeros(K,S);
TOP20_share_small = zeros(1,Nsmall);
TOP20_share_large = zeros(1,Nlarge);
TOP20_share = zeros(1,S);
Conf_lower = zeros(K,1);
Conf_upper = zeros(K,1);
Conf_lower_20 = zeros(K,1);
Conf_upper_20 = zeros(K,1);
R2 = zeros(K,1);
R2_20 = zeros(K,1);

X = XVEC_t(end,:);
D = DVEC_t(end,:);
ind = (X>0)&(D>0);

for z = 1:Nsmall
    ind_20_small = DSHMS_small(:,z)>=quantile(DSHMS_small(DSHMS_small(:,z)>0,z),0.8);
    TOP20_share_small(z) = sum(DSHMS_small(ind_20_small,z),1);
end

for z = 1:Nlarge
    ind_20_large = DSHMS_large(:,z)>=quantile(DSHMS_large(DSHMS_large(:,z)>0,z),0.8);
    TOP20_share_large(z) = sum(DSHMS_large(ind_20_large,z),1);
end

TOP20_share(small) = TOP20_share_small;
TOP20_share(~small) = TOP20_share_large;

for k = 1:K

    TOP(k,small) = sum(DSHMS_small(1:k,:),1);
    TOP(k,~small) = sum(DSHMS_large(1:k,:),1);

    stats = regstats(log(X(ind))',[TOP(k,ind)',log(D(ind))'],'linear',{'beta','covb','rsquare'});
    DynMom3_beta(k) = stats.beta(2);
    se = diag(sqrt(stats.covb));
    se = se(2);
    Conf_lower(k) = stats.beta(2)-1.980272249272974619185*se;
    Conf_upper(k) = stats.beta(2)+1.980272249272974619185*se;
    R2(k) = stats.rsquare;
    
    stats = regstats(log(X(ind))',[TOP(k,ind)',log(D(ind))',TOP20_share(ind)'],'linear',{'beta','covb','rsquare'});
    DynMom3_20_beta(k) = stats.beta(2); 
    se = diag(sqrt(stats.covb));
    se = se(2);
    Conf_lower_20(k) = stats.beta(2)-1.980272249272974619185*se;
    Conf_upper_20(k) = stats.beta(2)+1.980272249272974619185*se;
    R2_20(k) = stats.rsquare; 
end

% Generate data for Stata
DDATA = [X',D',TOP(1,:)',TOP(3,:)',TOP20_share'];
fname = sprintf('Results/Data/CS_regdata_%d',S);
fname2 = sprintf('Results/Data/CS_regdata_%d.csv',S);
save(fname,'DATA');
title1 = {'X','D','TOP1','TOP3','T20share'};
TTT = cell2table(num2cell(DDATA),'VariableNames',title1);
writetable(TTT,fname2);

% Plot Coefficients and R2

figure(8)
x=[1:K];
hold on
plot(x,DynMom3_beta','-or','LineWidth',3.0)
plot(x,DynMom3_20_beta,'-og','LineWidth',3.0)
plot(x,Conf_lower,'r')
plot(x,Conf_upper,'r')
plot(x,Conf_lower_20,'g')
plot(x,Conf_upper_20,'g')
legend('Only log(D) as control','Also Top20 share as control');
name = sprintf('Regression of log-exports on top-k firm market share, $rho = %d$',rho);
title(name,'interpreter','latex');
xlabel('k');
ylabel('Regression Coefficient');
name = sprintf('Results/coefficients_rho%d.png',rho);
saveas(gcf,name)

figure(9)
hold on 
plot(x,R2','-or','LineWidth',3.0)
plot(x,R2_20,'-og','LineWidth',3.0)
legend('Only log(D) as control','Also Top20 share as control');
name = sprintf('R2 of Regression of log-exports on top-k firm market share, $rho = %d$',rho);
title(name,'interpreter','latex');
xlabel('k');
ylabel('R2');
name = sprintf('Results/R2_rho%d.png',rho);
saveas(gcf,name)

