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
kappa=bestParams(4);
f=bestParams(5);

sigma=5;
theta=kappa*(sigma-1);
f=f*4.93*.43e-5;
F=f/sigma;

% New parameters for dynamic model
nu = 0.05;
mu = -theta*nu^2/2;


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
scale = 50;

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

% % Compute GE given initial productivities
[~,~,~,~,~,~,~,~,~,~,~,~,~,~,LAMBDAFVEC_0,PHIFVEC_0,~,varphi_BB]=PEreplication_vectorized(sigma,theta,F,tau,ALPHA,RTS,RTL,ZHS,ZFS,ZHL,ZFL,w,wF,Y,YF,small,vMU,BER,0,0);

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

figure(1)
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

figure(2)
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

figure(3)
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

figure(4)
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
RECORD = 1+[20,50];
R_length = length(RECORD);
T = RECORD(end);

LAMBDAFVEC_t = zeros(R_length+1,S);
LAMBDAFVEC_t(1,:) = LAMBDAFVEC_0;

PHIFVEC_t = zeros(R_length+1,S);
PHIFVEC_t(1,:) = PHIFVEC_0;

counter = 1;

for t=1:T
    ZHS(~inactive_HS) = exp(VB_HS(~inactive_HS)+abs(log(ZHS(~inactive_HS))-VB_HS(~inactive_HS)+mu+nu*randn(length(ZHS(~inactive_HS)),1)));
%     ZHS_Save =[ZHS_Save,ZHS(:)];
    ZHL(~inactive_HL) = exp(VB_HL(~inactive_HL)+abs(log(ZHL(~inactive_HL))-VB_HL(~inactive_HL)+mu+nu*randn(length(ZHL(~inactive_HL)),1)));
    ZFS(~inactive_FS) = exp(VB_FS(~inactive_FS)+abs(log(ZFS(~inactive_FS))-VB_FS(~inactive_FS)+mu+nu*randn(length(ZFS(~inactive_FS)),1)));
    ZFL(~inactive_FL) = exp(VB_FL(~inactive_FL)+abs(log(ZFL(~inactive_FL))-VB_FL(~inactive_FL)+mu+nu*randn(length(ZFL(~inactive_FL)),1)));
    
    % If the year is part of RECORD, record PE results
    if any(RECORD==t)
        counter = counter+1;
        [~,~,~,~,~,~,~,~,~,~,~,~,~,~,LAMBDAFVEC,PHIFVEC,~,~]=PEreplication_vectorized(sigma,theta,F,tau,ALPHA,RTS,RTL,ZHS,ZFS,ZHL,ZFL,w,wF,Y,YF,small,vMU,BER,0,0);
        LAMBDAFVEC_t(counter,:) = LAMBDAFVEC;
        PHIFVEC_t(counter,:) = PHIFVEC;
    end
    
end

%% MEAN REVERSION

% GAMMAF
GAMMAFVEC_t = LAMBDAFVEC_t-PHIFVEC_t;

% Delta LAMBDAF
DELTA_LAMBDAFVEC_t = LAMBDAFVEC_t-repmat(LAMBDAFVEC_t(1,:),R_length+1,1);

% Prepare plotting variables
ptiles = 10;
LPCT = prctile(GAMMAFVEC_t(1,:),[1/ptiles:1/ptiles:1-1/ptiles]*100)';

mean_20 = zeros(ptiles,1);          % TO DO: Automatically find 20 and 50 in RECORD
mean_50 = zeros(ptiles,1);

DELTA_20 = DELTA_LAMBDAFVEC_t(2,:);
DELTA_50 = DELTA_LAMBDAFVEC_t(3,:);

mean_20(1) = mean(DELTA_20(GAMMAFVEC_t(1,:)<LPCT(1)));
mean_50(1) = mean(DELTA_50(GAMMAFVEC_t(1,:)<LPCT(1)));
for i=2:ptiles-1
   mean_20(i) = mean(DELTA_20(GAMMAFVEC_t(1,:)<LPCT(i) & GAMMAFVEC_t(1,:)>=LPCT(i-1))); 
   mean_50(i) = mean(DELTA_50(GAMMAFVEC_t(1,:)<LPCT(i) & GAMMAFVEC_t(1,:)>=LPCT(i-1))); 
end
mean_20(ptiles) = mean(DELTA_20(GAMMAFVEC_t(1,:)>LPCT(ptiles-1)));
mean_50(ptiles) = mean(DELTA_50(GAMMAFVEC_t(1,:)>LPCT(ptiles-1)));

% Plot mean reversion graph
figure(5)
bar(mean_50,'FaceColor','[0, 0.4470, 0.7410]','FaceAlpha',0.33)
hold on
bar(mean_20,'FaceColor','[0.900, 0.20, 0.05]')%,'FaceAlpha',0.65)
lgd = legend('50 years','20 years')
set(lgd,'box','off','location','northeast','interpreter','latex','fontsize',20)
% ylim([0, .25])
% set(gca,'ytick',[0.05:.05:.25])
xlim([0.5, 10.55])
ylim([-.08, .04])
xlabel('Deciles of sectors, by granular $\Gamma_z^\ast$','interpreter','latex','fontsize',20)
ylabel('Expected change in export share, $\Delta \Lambda_z^\ast$','interpreter','latex','fontsize',20)
set(gca,'FontSize',16)
set(gca,'xtick',[0.5:1:9.5 10.55])
LPCT0 = LPCT;
LPCT0(5)=0;
set(gca,'xticklabel',[{ -1} round(LPCT0(1:4)'*100)/100 {  0} round(LPCT0(6:end)'*100)/100 1])






