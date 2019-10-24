% This code compute 30 target moments on a (rho,nu) grid.

%% Housekeeping
clear all;
close all;
clc;

tstart0 = tic;

addpath('Data');
addpath('Auxiliary Functions');

%% Parameter Vector from Static Model
load('GE_Results')

muT=bestParams(1);
sigmaT=bestParams(2);
tau=bestParams(3);
theta=bestParams(4);
F=bestParams(5);

sigma=5;

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

% Variable Threshold
var_thresh = 1;

% Scale of model (i.e. number of sectors)
scale = 21; %1 %3/4 %21

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

% Solve for GE quantities (can be skipped if computed previously and
% results imported above)
% [~,Y,YF,LF,~,~,~,~,~,~,LAMBDAFVEC_0,varphi_BB]=GEreplication_vectorized(sigma,theta,F,tau,ALPHA,RTS,RTL,ZHS,ZFS,ZHL,ZFL,w,wF,L0,Y0,YF0,small,vMU,BER,0);
% save('Data/GE_Results','bestParams','Y','YF','LF','varphi_BB')

% % Compute PE given initial productivities
[~,~,~,~,~,~,~,~,~,~,~,~,~,~,LAMBDAFVEC_0,PHIFVEC_0,~,varphi_BB,~,~]=PEreplication_vectorized(sigma,theta,F,tau,ALPHA,RTS,RTL,ZHS,ZFS,ZHL,ZFL,w,wF,Y,YF,small,vMU,BER,0,0);
disp('Loop 0 is finished')


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
net_intercept_HS = gross_intercept + 0.5*value_HS;
varphi_bar_HS = (1/theta)*log(RT(small))+net_intercept_HS;

% Large Sectors Home
[value_HL,~] = min(abs(log(min(varphi(1:2,Nsmall+1:end)))-log(min(ZHL))));
[~,zero_index_HL] = min(abs(log(RT(~small))));
gross_intercept = log(min(ZHL(:,zero_index_HL)));
net_intercept_HL = gross_intercept + 0.5*value_HL;
varphi_bar_HL = (1/theta)*log(RT(~small))+net_intercept_HL;

% Small Sectors Foreign
[value_FS,~] = min(abs(log(min(varphi(3:4,1:Nsmall)))-log(min(ZFS))));
[~,zero_index_FS] = min(abs(log(RT(small))));
gross_intercept = log(min(ZFS(:,zero_index_FS)));
net_intercept_FS = gross_intercept + 0.5*value_FS;
varphi_bar_FS = net_intercept_FS*ones(1,Nsmall);

% Large Sectors Foreign
[value_FL,~] = min(abs(log(min(varphi(3:4,Nsmall+1:end)))-log(min(ZFL))));
[~,zero_index_FL] = min(abs(log(RT(~small))));
gross_intercept = log(min(ZFL(:,zero_index_FL)));
net_intercept_FL = gross_intercept + 0.5*value_FL;
varphi_bar_FL = net_intercept_FL*ones(1,Nlarge);

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

rng(aseed);

% Years that are actually being recorded
RECORD = [10,20,50];
R_length = length(RECORD);
T = RECORD(end);

% Set up grid for nu and rho
alpha_u_vec = [0.07,0.1,0.075];
alpha_v_vec = [0.045,0.025,0.04];
% [alpha_v_mat,alpha_u_mat] = meshgrid(alpha_v_vec,alpha_u_vec); 
% [row,col] = size(alpha_v_mat);
% num_param = row*col;
num_param = length(alpha_u_vec);


ZHS_start = ZHS;
ZHL_start = ZHL;
ZFS_start = ZFS;
ZFL_start = ZFL;

VB_HS_start = VB_HS;
VB_HL_start = VB_HL;
VB_FS_start = VB_FS;
VB_FL_start = VB_FL;

Pareto = zeros(T,num_param);

mom = zeros(4,1);

% Iterate over this parameter grid
for itparam = 1:num_param
    
    % Set the correct parameters for this iteration
    alpha_v = alpha_v_vec(itparam);
    alpha_u = alpha_u_vec(itparam);
    mu = -theta*alpha_u^2/2;
    
    % Initialize matrices of interest
    LAMBDAFVEC_t = zeros(R_length+1,S);
    PHIFVEC_t = zeros(R_length+1,S);
    LAMBDAFVEC_t(1,:) = LAMBDAFVEC_0;
    PFIFVEC_t(1,:) = PHIFVEC_0;
    
    % Re-initialize matrices, so every iteration starts at the same
    % productivity 
    
    ZHS = ZHS_start;
    ZHL = ZHL_start;
    ZFS = ZFS_start;
    ZFL = ZFL_start;
    
    if var_thresh == 1
        VB_HS = VB_HS_start;
        VB_HL = VB_HL_start;
        VB_FS = VB_FS_start;
        VB_FL = VB_FL_start;
    end
    
    aseed = 1;
    rng(aseed);

    counter = 1;
    
    for t=1:T

        % Determine productivity draws
        vz_HS = repmat(randn(1,Nsmall),ZHS_length,1);
        vz_HL = repmat(randn(1,Nlarge),ZHL_length,1);
        vz_FS = repmat(randn(1,Nsmall),ZFS_length,1);
        vz_FL = repmat(randn(1,Nlarge),ZFL_length,1);

        uz_HS = randn(ZHS_length,Nsmall);
        uz_HL = randn(ZHL_length,Nlarge);
        uz_FS = randn(ZFS_length,Nsmall);
        uz_FL = randn(ZFL_length,Nlarge);

        eps_HS = alpha_u*uz_HS+alpha_v*vz_HS;
        eps_HL = alpha_u*uz_HL+alpha_v*vz_HL;
        eps_FS = alpha_u*uz_FS+alpha_v*vz_FS;
        eps_FL = alpha_u*uz_FL+alpha_v*vz_FL;
        
        if var_thresh == 1
            VB_HS = VB_HS + alpha_v * vz_HS;
            VB_HL = VB_HL + alpha_v * vz_HL;
            VB_FS = VB_FS + alpha_v * vz_FS;
            VB_FL = VB_FL + alpha_v * vz_FL;
        end

        % Evolution of productivity draws
        ZHS(~inactive_HS) = exp(VB_HS(~inactive_HS)+abs(log(ZHS(~inactive_HS))-VB_HS(~inactive_HS)+mu+eps_HS(~inactive_HS)));
        ZHL(~inactive_HL) = exp(VB_HL(~inactive_HL)+abs(log(ZHL(~inactive_HL))-VB_HL(~inactive_HL)+mu+eps_HL(~inactive_HL)));
        ZFS(~inactive_FS) = exp(VB_FS(~inactive_FS)+abs(log(ZFS(~inactive_FS))-VB_FS(~inactive_FS)+mu+eps_FS(~inactive_FS)));
        ZFL(~inactive_FL) = exp(VB_FL(~inactive_FL)+abs(log(ZFL(~inactive_FL))-VB_FL(~inactive_FL)+mu+eps_FL(~inactive_FL)));

        % If the year is part of RECORD, record PE results
        if any(RECORD==t)
            counter = counter+1;
            [~,~,~,~,~,~,~,~,~,~,~,~,~,~,LAMBDAFVEC,PHIFVEC,~,~,XVEC,DVEC,DSHM_small,DSHM_large,DSHMS_small,DSHMS_large]=PEreplication_vectorized(sigma,theta,F,tau,ALPHA,RTS,RTL,ZHS,ZFS,ZHL,ZFL,w,wF,Y,YF,small,vMU,BER,0,0);
            LAMBDAFVEC_t(counter,:) = LAMBDAFVEC;
            PHIFVEC_t(counter,:) = PHIFVEC;
            disp(['Loop ',num2str(t),' is finished']);
        end
    
    end
    disp(['Iteration ',num2str(itparam),' is finished']);
    
    %% MEAN REVERSION
    
    % GAMMAF
    GAMMAFVEC_t = LAMBDAFVEC_t-PHIFVEC_t;
    
    % Delta LAMBDAF
    DELTA_LAMBDAFVEC_t = LAMBDAFVEC_t-repmat(LAMBDAFVEC_t(1,:),R_length+1,1);
    
    % Prepare plotting variables
    ptiles = 10;
    LPCT = prctile(GAMMAFVEC_t(1,:),[1/ptiles:1/ptiles:1-1/ptiles]*100)';
    
    mean_20 = zeros(ptiles,1);         
    mean_50 = zeros(ptiles,1);
    
    DELTA_20 = DELTA_LAMBDAFVEC_t(3,:);
    DELTA_50 = DELTA_LAMBDAFVEC_t(4,:);
    
    mean_20(1) = mean(DELTA_20(GAMMAFVEC_t(1,:)<LPCT(1)));
    mean_50(1) = mean(DELTA_50(GAMMAFVEC_t(1,:)<LPCT(1)));

    for i=2:ptiles-1
       mean_20(i) = mean(DELTA_20(GAMMAFVEC_t(1,:)<LPCT(i) & GAMMAFVEC_t(1,:)>=LPCT(i-1))); 
       mean_50(i) = mean(DELTA_50(GAMMAFVEC_t(1,:)<LPCT(i) & GAMMAFVEC_t(1,:)>=LPCT(i-1))); 
    end
    mean_20(ptiles) = mean(DELTA_20(GAMMAFVEC_t(1,:)>=LPCT(ptiles-1)));
    mean_50(ptiles) = mean(DELTA_50(GAMMAFVEC_t(1,:)>=LPCT(ptiles-1)));
    
    % Plot mean reversion graph
    figure(itparam)
    bar(mean_50,'FaceColor','[0, 0.4470, 0.7410]','FaceAlpha',0.33)
    hold on
    bar(mean_20,'FaceColor','[0.900, 0.20, 0.05]')%,'FaceAlpha',0.65)
    lgd = legend('50 years','20 years');
    set(lgd,'box','off','location','northeast','interpreter','latex','fontsize',20)
    % ylim([0, .25])
    % set(gca,'ytick',[0.05:.05:.25])
    xlim([0.5, 10.55])
    ylim([-.15, .2])
    xlabel('Deciles of sectors, by granular $\Gamma_z^\ast$','interpreter','latex','fontsize',20)
    ylabel('Expected change in export share, $\Delta \Lambda_z^\ast$','interpreter','latex','fontsize',20)
    set(gca,'FontSize',16)
    set(gca,'xtick',[0.5:1:9.5 10.55])
    LPCT0 = LPCT;
    LPCT0(5)=0;
    set(gca,'xticklabel',[{ -1} round(LPCT0(1:4)'*100)/100 {  0} {  0} round(LPCT0(7:end)'*100)/100 1])
    ffname = sprintf('Results/Calibrating_Graphs/mean_reversion%d.png',itparam);
    saveas(gcf,ffname)
    
    idx_0 = LAMBDAFVEC_t(1,:)>prctile(LAMBDAFVEC_t(1,:),95);
    idx_10 = LAMBDAFVEC_t(2,:)>prctile(LAMBDAFVEC_t(2,:),95);
    idx_20 = LAMBDAFVEC_t(3,:)>prctile(LAMBDAFVEC_t(3,:),95);
    
    mom(1:2,1) = [alpha_v;alpha_u];
    mom(3,1) = sum(idx_0 & idx_10)/sum(idx_0);
    mom(4,1) = sum(idx_0 & idx_20)/sum(idx_0);
    fname = sprintf('Results/Calibrating_Graphs/turnover%d.csv',itparam);
    TTT = cell2table(num2cell(mom));
    writetable(TTT,fname);
    
end

