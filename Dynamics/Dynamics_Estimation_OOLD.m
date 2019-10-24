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
scale = 3/4; %1 %3/4 %21

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
RECORD = 1:11;
R_length = length(RECORD);
T = RECORD(end);

% Set up grid for nu and rho
alpha_v_vec = 0.035:0.005:0.05;
alpha_u_vec = [0.07,0.075,0.08:0.01:0.12];
[alpha_v_mat,alpha_u_mat] = meshgrid(alpha_v_vec,alpha_u_vec); 
[row,col] = size(alpha_v_mat);
num_param = row*col;

DATA = zeros(38,num_param);

ZHS_start = ZHS;
ZHL_start = ZHL;
ZFS_start = ZFS;
ZFL_start = ZFL;

VB_HS_start = VB_HS;
VB_HL_start = VB_HL;
VB_FS_start = VB_FS;
VB_FL_start = VB_FL;

Pareto = zeros(T,num_param);

% Iterate over this parameter grid
for itparam = 1:num_param
    
    % Set the correct parameters for this iteration
    alpha_v = alpha_v_mat(itparam);
    alpha_u = alpha_u_mat(itparam);
    mu = -theta*alpha_u^2/2;
    
    % Initialize matrices of interest
    LAMBDAFVEC_t = zeros(R_length,S);
    DVEC_t = zeros(R_length,S);
    XVEC_t = zeros(R_length,S);
    TOP1_t = zeros(R_length,S);
    TOP3_t = zeros(R_length,S);
    PHIFVEC_t = zeros(R_length,S);
    
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
            [~,~,~,~,~,~,~,~,~,~,~,~,~,~,LAMBDAFVEC,PHIFVEC,~,~,XVEC,DVEC,DSHM_small,DSHM_large,DSHMS_small,DSHMS_large]=PEreplication_vectorized(sigma,theta,F,tau,ALPHA,RTS,RTL,ZHS,ZFS,ZHL,ZFL,w,wF,Y,YF,small,vMU,BER,0,0);
            LAMBDAFVEC_t(t,:) = LAMBDAFVEC;
            PHIFVEC_t(t,:) = PHIFVEC;
            XVEC_t(t,:) = XVEC; 
            DVEC_t(t,:) = DVEC;
            DSHM_small = DSHM_small.*repmat(sum(DSHM_small),size(DSHM_small,1),1);
            DSHM_large = DSHM_large.*repmat(sum(DSHM_large),size(DSHM_large,1),1);
            save_share_small{t} = DSHM_small;
            save_share_large{t} = DSHM_large;
            save_share_small_S{t} = DSHMS_small;
            save_share_large_S{t} = DSHMS_large;
            save_share_small_all{t} = DSHM_small;
            save_share_large_all{t} = DSHM_large; 
            disp(['Loop ',num2str(t),' is finished']);
        end
    
    end
    disp(['Iteration ',num2str(itparam),' is finished']);
    
    idx_active_small = true(size(DSHM_small));
    idx_active_large = true(size(DSHM_large));
    
    for t=1:T
       idx_active_small = idx_active_small & (save_share_small{t}>0);
       idx_active_large = idx_active_large & (save_share_large{t}>0);
    end
    
%     save_share_small_X = save_share_small;
%     save_share_large_X = save_share_large;
    for z=1:Nsmall
        for t=1:T
            int_save_share = zeros(length(save_share_small{t}(:,z)),1);
            int_save_share(idx_active_small(:,z)) = save_share_small{t}(idx_active_small(:,z),z);
            save_share_small{t}(:,z) = int_save_share;
        end
    end
    
    for z=1:Nlarge
        for t=1:T
            int_save_share = zeros(length(save_share_large{t}(:,z)),1);
            int_save_share(idx_active_large(:,z)) = save_share_large{t}(idx_active_large(:,z),z);
            save_share_large{t}(:,z) = int_save_share;
        end
    end
    
    for t=1:T
        save_share_small_all{t} = save_share_small{t}./repmat(sum(save_share_small_all{t}),size(save_share_small_all{t},1),1);
        save_share_large_all{t} = save_share_large{t}./repmat(sum(save_share_large_all{t}),size(save_share_large_all{t},1),1);
        save_share_small{t} = save_share_small{t}./repmat(sum(save_share_small{t}),size(save_share_small{t},1),1);
        save_share_large{t} = save_share_large{t}./repmat(sum(save_share_large{t}),size(save_share_large{t},1),1);
        TOP1_t(t,small) = save_share_small{t}(1,:);
        TOP3_t(t,small) = sum(save_share_small{t}(1:3,:));
        TOP1_t(t,~small) = save_share_large{t}(1,:);
        TOP3_t(t,~small) = sum(save_share_large{t}(1:3,:));
    end
    
%     %%% Save Regression data for Stata
%     ID = (1:S);
%     ID = repmat(ID,R_length,1);
%     YEAR = repmat(RECORD',S,1);
%     DATA1 = [ID(:),YEAR,XVEC_t(:),DVEC_t(:),TOP1_t(:),TOP3_t(:)];
%     fname = sprintf('Results/Data/calibrated_regdata_%d',S);
%     fname2 = sprintf('Results/Data/calibrated_regdata_%d.csv',S);
%     save(fname,'DATA');
%     title1 = {'ID','Year','X','D','TOP1','TOP3'};
%     TT = cell2table(num2cell(DATA1),'VariableNames',title1); 
%     writetable(TT,fname2);
    
    %%%% Compute 30 target moments
    
    
    %%% Moments 1 & 2
    
    Dep = XVEC_t(2:end,:);
    Indep = XVEC_t(1:end-1,:);
    
    % Correcting for zeros
    ind = (Dep>0) & (Indep>0);
    Dep = log(Dep(ind));
    Indep = log(Indep(ind));
    
    % Regression w/o fixed effects
    stats = regstats(Dep,Indep,'linear',{'beta','covb'});
    DATA(1,itparam) = stats.beta(2);
    DATA(2,itparam) = std(Dep-stats.beta(1)-stats.beta(2)*Indep);
    
    %%% Moments 3 & 4
    
    Dep = XVEC_t(2:end,:)./DVEC_t(2:end,:);
    Indep = XVEC_t(1:end-1,:)./DVEC_t(1:end-1,:);
    
    % Correcting for zeros
    ind = (XVEC_t(2:end,:)>0) & (DVEC_t(2:end,:)>0) & (XVEC_t(1:end-1,:)>0) & (DVEC_t(1:end-1,:)>0);
    Dep = Dep(ind);
    Indep = Indep(ind);
    
    % Regression w/o fixed effects
    stats = regstats(Dep,Indep,'linear',{'beta','covb'});
    DATA(3,itparam) = stats.beta(2);
    DATA(4,itparam) = std(Dep-stats.beta(1)-stats.beta(2)*Indep);    
    
    %%% Moments 5 & 6
    
    % Pick only those firms who are active in all years
    ind_small = true(size(DSHM_small));
    ind_large = true(size(DSHM_large));
    
    for t=1:T
        ind_small = ind_small & (save_share_small{t}>0);
        ind_large = ind_large & (save_share_large{t}>0);
    end
    
    clearvars small_shares large_shares
    for t=1:T
       small_shares(t,:) = save_share_small{t}(ind_small)';
       large_shares(t,:) = save_share_large{t}(ind_large)';
    end
    
    overall_shares = [small_shares,large_shares];
    
    Dep = overall_shares(2:end,:);
    Indep = overall_shares(1:end-1,:);
    
    % Regression w/o fixed effects
    stats = regstats(Dep(:),Indep(:),'linear',{'beta','covb'});
    DATA(5,itparam) = stats.beta(2);
    DATA(6,itparam) = std(Dep(:)-stats.beta(1)-stats.beta(2)*Indep(:));    
    
    %%% Moments 7 & 8
    
    small_shares = zeros(T,0);
    for z=1:Nsmall
        clearvars small_shares1
        % Indicator for being active in all periods
        dim = size(DSHM_small);
        idx_z_active = true(dim(1),1);
        for t=1:T
           idx_z_active = idx_z_active & (save_share_small{t}(:,z)>0);
        end
        % Indicator for belonging to the top 20% among all active firms in
        % year 1
        idx20 = save_share_small{1}(idx_z_active,z)>=quantile(save_share_small{1}(idx_z_active,z),0.8);
        for t=1:T
           int = save_share_small{t}(idx_z_active,z);
           int2 = int(idx20);
           small_shares1(t,:) = int2';
        end
        % This is the target matrix that contains a matrix with T rows and
        % where each column represents the development of one active firm
        % in any sector over time
        small_shares = [small_shares,small_shares1];
    end
    
    large_shares = zeros(T,0);
    for z=1:Nlarge
        clearvars large_shares1
        dim = size(DSHM_large);
        idx_z_active = true(dim(1),1);
        for t=1:T
           idx_z_active = idx_z_active & (save_share_large{t}(:,z)>0); 
        end
        idx20 = save_share_large{1}(idx_z_active,z)>=quantile(save_share_large{1}(idx_z_active,z),0.8);
        for t=1:T
           int = save_share_large{t}(idx_z_active,z);
           int2 = int(idx20);
           large_shares1(t,:) = int2';
        end
        large_shares = [large_shares,large_shares1];
    end
    
    overall_shares = [small_shares,large_shares];
    
    Dep = overall_shares(2:end,:);
    Indep = overall_shares(1:end-1,:);
    
    % Regression w/o fixed effects
    stats = regstats(Dep(:),Indep(:),'linear',{'beta','covb'});
    DATA(7,itparam) = stats.beta(2);
    DATA(8,itparam) = std(Dep(:)-stats.beta(1)-stats.beta(2)*Indep(:));    
    
    
    %%% Moments 9 & 10
    
    XX = XVEC_t(2:end,:);
    XX1 = XVEC_t(1:end-1,:);
    
    clearvars Dep Indep ind
    % Demeaning to control for fixed effects

    for z=1:S
        ind = (XVEC_t(2:end,z)>0) & (XVEC_t(1:end-1,z)>0);
        if sum(ind)>0
            Dep(ind,z) = log(XX(ind,z))-mean(log(XX(ind,z)));
            Indep(ind,z) = log(XX1(ind,z))-mean(log(XX1(ind,z)));
        end
    end

    % Regression 
    stats = regstats(Dep(:),Indep(:),'linear',{'beta','covb'});
    DATA(9,itparam) = stats.beta(2);
    DATA(10,itparam) = std(Dep(:)-stats.beta(2)*Indep(:));

    %%% Moments 11 & 12
    
    XX = XVEC_t(2:end,:)./DVEC_t(2:end,:);
    XX1 = XVEC_t(1:end-1,:)./DVEC_t(1:end-1,:);
    
    clearvars Dep Indep ind
    % Demeaning to control for fixed effects
    for z=1:S
        ind = (XVEC_t(2:end,z)>0)&(DVEC_t(2:end,z)>0)&(XVEC_t(1:end-1,z)>0)&(DVEC_t(1:end-1,z)>0);
        if sum(ind)>0
            Dep(ind,z) = XX(ind,z)-mean(XX(ind,z));
            Indep(ind,z) = XX1(ind,z)-mean(XX1(ind,z));
        end
    end

    % Regression
    stats = regstats(Dep(:),Indep(:),'linear',{'beta','covb'});
    DATA(11,itparam) = stats.beta(2);
    DATA(12,itparam) = std(Dep(:)-stats.beta(2)*Indep(:));
  
    
    %%% Moments 13 - 16
    
    small_shares_dep = zeros(T-1,0);
    small_shares_indep = zeros(T-1,0);
    small_shares_dep_f = zeros(T-1,0);
    small_shares_indep_f = zeros(T-1,0);
    for z=1:Nsmall
        clearvars small_shares1 small_shares2
        % Indicator for being active in all periods
        dim = size(DSHM_small);
        idx_z_active = true(dim(1),1);
        for t=1:T
           idx_z_active = idx_z_active & (save_share_small{t}(:,z)>0);
        end
        % Indicator for belonging to the top 20% among all active firms in
        % year 1
        idx20 = save_share_small{1}(idx_z_active,z)>=quantile(save_share_small{1}(idx_z_active,z),0.8);
        for t=1:T
           int = save_share_small{t}(idx_z_active,z);
           int2 = int(idx20);
           small_shares1(t,:) = int';
           small_shares2(t,:) = int2';
        end
        small_shares1_dep = small_shares1(2:end,:) - mean(mean(small_shares1(2:end,:)));
        small_shares1_indep = small_shares1(1:end-1,:) - mean(mean(small_shares1(1:end-1,:)));
        small_shares2_dep = small_shares2(2:end,:) - mean(mean(small_shares2(2:end,:)));
        small_shares2_indep = small_shares2(1:end-1,:) - mean(mean(small_shares2(1:end-1,:)));
        % This is the target matrix that contains a matrix with T rows and
        % where each column represents the development of one active firm
        % in any sector over time
        small_shares_dep = [small_shares_dep,small_shares1_dep];
        small_shares_indep = [small_shares_indep,small_shares1_indep];
        small_shares_dep_f = [small_shares_dep_f,small_shares2_dep];
        small_shares_indep_f = [small_shares_indep_f,small_shares2_indep];
    end
    
    large_shares_dep = zeros(T-1,0);
    large_shares_indep = zeros(T-1,0);
    large_shares_dep_f = zeros(T-1,0);
    large_shares_indep_f = zeros(T-1,0);
    for z=1:Nlarge
        clearvars large_shares1 large_shares2
        dim = size(DSHM_large);
        idx_z_active = true(dim(1),1);
        for t=1:T
           idx_z_active = idx_z_active & (save_share_large{t}(:,z)>0); 
        end
        idx20 = save_share_large{1}(idx_z_active,z)>=quantile(save_share_large{1}(idx_z_active,z),0.8);
        for t=1:T
           int = save_share_large{t}(idx_z_active,z);
           int2 = int(idx20);
           large_shares1(t,:) = int';
           large_shares2(t,:) = int2';
        end
        large_shares1_dep = large_shares1(2:end,:) - mean(mean(large_shares1(2:end,:)));
        large_shares1_indep = large_shares1(1:end-1,:) - mean(mean(large_shares1(1:end-1,:)));
        large_shares2_dep = large_shares2(2:end,:) - mean(mean(large_shares2(2:end,:)));
        large_shares2_indep = large_shares2(1:end-1,:) - mean(mean(large_shares2(1:end-1,:)));
        large_shares_dep = [large_shares_dep,large_shares1_dep];
        large_shares_indep = [large_shares_indep,large_shares1_indep];
        large_shares_dep_f = [large_shares_dep_f,large_shares2_dep];
        large_shares_indep_f = [large_shares_indep_f,large_shares2_indep];
    end
    
    Dep = [small_shares_dep,large_shares_dep];
    Indep = [small_shares_indep,large_shares_indep];
    
    Dep_f = [small_shares_dep_f,large_shares_dep_f];
    Indep_f = [small_shares_indep_f,large_shares_indep_f];
    
    % Regression
    stats = regstats(Dep(:),Indep(:),'linear',{'beta','covb'});
    DATA(13,itparam) = stats.beta(2);
    DATA(14,itparam) = std(Dep(:)-stats.beta(2)*Indep(:));
    
    stats = regstats(Dep_f(:),Indep_f(:),'linear',{'beta','covb'});
    DATA(15,itparam) = stats.beta(2);
    DATA(16,itparam) = std(Dep_f(:)-stats.beta(2)*Indep_f(:));
    
    % Moments 17 - 22
    k_vec = [5,10];
    mom = zeros(length(k_vec),3);
    counter = 0;
    for k = k_vec
        Dep = zeros(0);
        Indep = zeros(0);
        Dep_f = zeros(0);
        Indep_f = zeros(0);
        Control = zeros(0);
        counter = counter + 1;
        for t=1:T-k
            Dep1 = log(XVEC_t(t+k,:))-log(XVEC_t(t,:));
            Dep2 = (XVEC_t(t+k,:)./DVEC_t(t+k,:))-(XVEC_t(t,:)./DVEC_t(t,:));
            Indep1 = log(XVEC_t(t,:));
            Indep2 = XVEC_t(t,:)./DVEC_t(t,:);
            Control1 = log(DVEC_t(t,:));
            ind = (XVEC_t(t+k,:)>0) & (XVEC_t(t,:)>0);
            ind2 = (XVEC_t(t+k,:)>0) & (DVEC_t(t+k,:)>0) & (XVEC_t(t,:)>0) & (DVEC_t(t,:)>0);
            Dep1 = Dep1(ind);
            Dep2 = Dep2(ind2);
            Indep1 = Indep1(ind);
            Indep2 = Indep2(ind2);
            Control1 = Control1(ind);
            Dep = [Dep,Dep1];
            Indep = [Indep,Indep1]; 
            Dep_f = [Dep_f,Dep2];
            Indep_f = [Indep_f,Indep2];
            Control = [Control,Control1];
        end
        % Regression w/o controls
        stats = regstats(Dep',Indep','linear',{'beta','covb'});
        mom(counter,1) = stats.beta(2);

        % Regression w/ controls
        IndCont = [Indep', Control'];
        stats = regstats(Dep',[Indep', Control'],'linear',{'beta','covb'});
        mom(counter,2) = stats.beta(2);
        
        % Regression for moments 21 & 22
        stats = regstats(Dep_f',Indep_f','linear',{'beta','covb'});
        mom(counter,3) = stats.beta(2);
        
    end
    DATA(17,itparam) = mom(1,1);
    DATA(18,itparam) = mom(2,1);
    DATA(19,itparam) = mom(1,2);
    DATA(20,itparam) = mom(2,2);
    DATA(21,itparam) = mom(1,3);
    DATA(22,itparam) = mom(2,3);
    
    
    % Moments 23 - 26
    mom = zeros(length(k_vec),2);
    counter = 0;
    for k = k_vec
        
        small_shares_dep = zeros(T-k,0);
        small_shares_indep = zeros(T-k,0);
        small_shares_dep_f = zeros(T-k,0);
        small_shares_indep_f = zeros(T-k,0);
        counter = counter + 1;
        for z=1:Nsmall
            clearvars small_shares1_dep small_shares2_dep small_shares1_indep small_shares2_indep
            % Indicator for being active in all periods
            dim = size(DSHM_small);
            idx_z_active = true(dim(1),1);
            for t=1:T
               idx_z_active = idx_z_active & (save_share_small{t}(:,z)>0);
            end
            % Indicator for belonging to the top 20% among all active firms in
            % year 1
            idx20 = save_share_small{1}(idx_z_active,z)>=quantile(save_share_small{1}(idx_z_active,z),0.8);
            for t=1:T-k
               int_dep = save_share_small{t+k}(idx_z_active,z) - save_share_small{t}(idx_z_active,z);
               int2_dep = int_dep(idx20);
               int_indep = save_share_small{t}(idx_z_active,z);
               int2_indep = int_indep(idx20);
               small_shares1_dep(t,:) = int_dep';
               small_shares1_indep(t,:) = int_indep';
               small_shares2_dep(t,:) = int2_dep';
               small_shares2_indep(t,:) = int2_indep';
            end
            % This is the target matrix that contains a matrix with T rows and
            % where each column represents the development of one active firm
            % in any sector over time
            small_shares_dep = [small_shares_dep,small_shares1_dep];
            small_shares_indep = [small_shares_indep,small_shares1_indep];
            small_shares_dep_f = [small_shares_dep_f,small_shares2_dep];
            small_shares_indep_f = [small_shares_indep_f,small_shares2_indep];
        end

        large_shares_dep = zeros(T-k,0);
        large_shares_indep = zeros(T-k,0);
        large_shares_dep_f = zeros(T-k,0);
        large_shares_indep_f = zeros(T-k,0);
        for z=1:Nlarge
            clearvars large_shares1_dep large_shares2_dep large_shares1_indep large_shares2_indep 
            dim = size(DSHM_large);
            idx_z_active = true(dim(1),1);
            for t=1:T
               idx_z_active = idx_z_active & (save_share_large{t}(:,z)>0); 
            end
            idx20 = save_share_large{1}(idx_z_active,z)>=quantile(save_share_large{1}(idx_z_active,z),0.8);
            for t=1:T-k
                int_dep = save_share_large{t+k}(idx_z_active,z) - save_share_large{t}(idx_z_active,z);
                int2_dep = int_dep(idx20);
                int_indep = save_share_large{t}(idx_z_active,z);
                int2_indep = int_indep(idx20);
                large_shares1_dep(t,:) = int_dep';
                large_shares1_indep(t,:) = int_indep';
                large_shares2_dep(t,:) = int2_dep';
                large_shares2_indep(t,:) = int2_indep';
            end
            large_shares_dep = [large_shares_dep,large_shares1_dep];
            large_shares_indep = [large_shares_indep,large_shares1_indep];
            large_shares_dep_f = [large_shares_dep_f,large_shares2_dep];
            large_shares_indep_f = [large_shares_indep_f,large_shares2_indep];
        end

        Dep = [small_shares_dep,large_shares_dep];
        Indep = [small_shares_indep,large_shares_indep];

        Dep_f = [small_shares_dep_f,large_shares_dep_f];
        Indep_f = [small_shares_indep_f,large_shares_indep_f];

        % Regression
        stats = regstats(Dep(:),Indep(:),'linear',{'beta','covb'});
        mom(counter,1) = stats.beta(2);

        stats = regstats(Dep_f(:),Indep_f(:),'linear',{'beta','covb'});
        mom(counter,2) = stats.beta(2);
    end
    
    DATA(23,itparam) = mom(1,1);
    DATA(24,itparam) = mom(2,1);
    DATA(25,itparam) = mom(1,2);
    DATA(26,itparam) = mom(2,2);
    
        
    % Moments 23' - 26'
    mom = zeros(length(k_vec),2);
    counter = 0;
    for k = k_vec
        
        small_shares_dep = zeros(T-k,0);
        small_shares_indep = zeros(T-k,0);
        small_shares_dep_f = zeros(T-k,0);
        small_shares_indep_f = zeros(T-k,0);
        counter = counter + 1;
        for z=1:Nsmall
            clearvars small_shares1_dep small_shares2_dep small_shares1_indep small_shares2_indep
            % Indicator for being active in all periods
            dim = size(DSHM_small);
            idx_z_active = true(dim(1),1);
            for t=1:T
               idx_z_active = idx_z_active & (save_share_small_all{t}(:,z)>0);
            end
            % Indicator for belonging to the top 20% among all active firms in
            % year 1
            idx20 = save_share_small_all{1}(idx_z_active,z)>=quantile(save_share_small_all{1}(idx_z_active,z),0.8);
            for t=1:T-k
               int_dep = save_share_small_all{t+k}(idx_z_active,z) - save_share_small_all{t}(idx_z_active,z);
               int2_dep = int_dep(idx20);
               int_indep = save_share_small_all{t}(idx_z_active,z);
               int2_indep = int_indep(idx20);
               small_shares1_dep(t,:) = int_dep';
               small_shares1_indep(t,:) = int_indep';
               small_shares2_dep(t,:) = int2_dep';
               small_shares2_indep(t,:) = int2_indep';
            end
            % This is the target matrix that contains a matrix with T rows and
            % where each column represents the development of one active firm
            % in any sector over time
            small_shares_dep = [small_shares_dep,small_shares1_dep];
            small_shares_indep = [small_shares_indep,small_shares1_indep];
            small_shares_dep_f = [small_shares_dep_f,small_shares2_dep];
            small_shares_indep_f = [small_shares_indep_f,small_shares2_indep];
        end

        large_shares_dep = zeros(T-k,0);
        large_shares_indep = zeros(T-k,0);
        large_shares_dep_f = zeros(T-k,0);
        large_shares_indep_f = zeros(T-k,0);
        for z=1:Nlarge
            clearvars large_shares1_dep large_shares2_dep large_shares1_indep large_shares2_indep 
            dim = size(DSHM_large);
            idx_z_active = true(dim(1),1);
            for t=1:T
               idx_z_active = idx_z_active & (save_share_large_all{t}(:,z)>0); 
            end
            idx20 = save_share_large_all{1}(idx_z_active,z)>=quantile(save_share_large_all{1}(idx_z_active,z),0.8);
            for t=1:T-k
                int_dep = save_share_large_all{t+k}(idx_z_active,z) - save_share_large_all{t}(idx_z_active,z);
                int2_dep = int_dep(idx20);
                int_indep = save_share_large_all{t}(idx_z_active,z);
                int2_indep = int_indep(idx20);
                large_shares1_dep(t,:) = int_dep';
                large_shares1_indep(t,:) = int_indep';
                large_shares2_dep(t,:) = int2_dep';
                large_shares2_indep(t,:) = int2_indep';
            end
            large_shares_dep = [large_shares_dep,large_shares1_dep];
            large_shares_indep = [large_shares_indep,large_shares1_indep];
            large_shares_dep_f = [large_shares_dep_f,large_shares2_dep];
            large_shares_indep_f = [large_shares_indep_f,large_shares2_indep];
        end

        Dep = [small_shares_dep,large_shares_dep];
        Indep = [small_shares_indep,large_shares_indep];

        Dep_f = [small_shares_dep_f,large_shares_dep_f];
        Indep_f = [small_shares_indep_f,large_shares_indep_f];

        % Regression
        stats = regstats(Dep(:),Indep(:),'linear',{'beta','covb'});
        mom(counter,1) = stats.beta(2);

        stats = regstats(Dep_f(:),Indep_f(:),'linear',{'beta','covb'});
        mom(counter,2) = stats.beta(2);
    end
    
    DATA(27,itparam) = mom(1,1);
    DATA(28,itparam) = mom(2,1);
    DATA(29,itparam) = mom(1,2);
    DATA(30,itparam) = mom(2,2);
    
    % Moment 27
    
    Delta = log(XVEC_t(2:end,:))-log(XVEC_t(1:end-1,:));
    ind = (XVEC_t(2:end,:)>0) & (XVEC_t(1:end-1,:)>0);
    Delta = Delta(ind);
    DATA(31,itparam) = std(Delta);
    
    % Moment 28
    
    Delta = (XVEC_t(2:end,:)./DVEC_t(2:end,:))-(XVEC_t(1:end-1,:)./DVEC_t(1:end-1,:));
    ind = (XVEC_t(2:end,:)>0) & (DVEC_t(2:end,:)>0) & (XVEC_t(1:end-1,:)>0) & (DVEC_t(1:end-1,:)>0);
    Delta = Delta(ind);
    DATA(32,itparam) = std(Delta);
    
    % Moments 29 & 30
    
    small_shares = zeros(T-1,0);
    small_shares_f = zeros(T-1,0);
    for z=1:Nsmall
        clearvars small_shares1 small_shares2
        % Indicator for being active in all periods
        dim = size(DSHM_small);
        idx_z_active = true(dim(1),1);
        for t=1:T
           idx_z_active = idx_z_active & (save_share_small{t}(:,z)>0);
        end
        % Indicator for belonging to the top 20% among all active firms in
        % year 1
        idx20 = save_share_small{1}(idx_z_active,z)>=quantile(save_share_small{1}(idx_z_active,z),0.8);
        for t=1:T-1
           int = save_share_small{t+1}(idx_z_active,z)-save_share_small{t}(idx_z_active,z);
           int2 = int(idx20);
           small_shares1(t,:) = int';
           small_shares2(t,:) = int2';
        end
        % This is the target matrix that contains a matrix with T rows and
        % where each column represents the development of one active firm
        % in any sector over time
        small_shares = [small_shares,small_shares1];
        small_shares_f = [small_shares_f,small_shares2];
    end
    
    large_shares = zeros(T-1,0);
    large_shares_f = zeros(T-1,0);
    for z=1:Nlarge
        clearvars large_shares1 large_shares2
        dim = size(DSHM_large);
        idx_z_active = true(dim(1),1);
        for t=1:T
           idx_z_active = idx_z_active & (save_share_large{t}(:,z)>0); 
        end
        idx20 = save_share_large{1}(idx_z_active,z)>=quantile(save_share_large{1}(idx_z_active,z),0.8);
        for t=1:T-1
           int = save_share_large{t+1}(idx_z_active,z)-save_share_large{t}(idx_z_active,z);
           int2 = int(idx20);
           large_shares1(t,:) = int';
           large_shares2(t,:) = int2';
        end
        large_shares = [large_shares,large_shares1];
        large_shares_f = [large_shares_f,large_shares2]; 
    end
    
    Delta = [small_shares,large_shares];
    Delta_f = [small_shares_f,large_shares_f];
    
    DATA(33,itparam) = std(Delta(:));
    DATA(34,itparam) = std(Delta_f(:));
    
    % Moments 29' & 30'
    
    small_shares = zeros(T-1,0);
    small_shares_f = zeros(T-1,0);
    for z=1:Nsmall
        clearvars small_shares1 small_shares2
        % Indicator for being active in all periods
        dim = size(DSHM_small);
        idx_z_active = true(dim(1),1);
        for t=1:T
           idx_z_active = idx_z_active & (save_share_small_all{t}(:,z)>0);
        end
        % Indicator for belonging to the top 20% among all active firms in
        % year 1
        idx20 = save_share_small_all{1}(idx_z_active,z)>=quantile(save_share_small_all{1}(idx_z_active,z),0.8);
        for t=1:T-1
           int = save_share_small_all{t+1}(idx_z_active,z)-save_share_small_all{t}(idx_z_active,z);
           int2 = int(idx20);
           small_shares1(t,:) = int';
           small_shares2(t,:) = int2';
        end
        % This is the target matrix that contains a matrix with T rows and
        % where each column represents the development of one active firm
        % in any sector over time
        small_shares = [small_shares,small_shares1];
        small_shares_f = [small_shares_f,small_shares2];
    end
    
    large_shares = zeros(T-1,0);
    large_shares_f = zeros(T-1,0);
    for z=1:Nlarge
        clearvars large_shares1 large_shares2
        dim = size(DSHM_large);
        idx_z_active = true(dim(1),1);
        for t=1:T
           idx_z_active = idx_z_active & (save_share_large_all{t}(:,z)>0); 
        end
        idx20 = save_share_large_all{1}(idx_z_active,z)>=quantile(save_share_large_all{1}(idx_z_active,z),0.8);
        for t=1:T-1
           int = save_share_large_all{t+1}(idx_z_active,z)-save_share_large_all{t}(idx_z_active,z);
           int2 = int(idx20);
           large_shares1(t,:) = int';
           large_shares2(t,:) = int2';
        end
        large_shares = [large_shares,large_shares1];
        large_shares_f = [large_shares_f,large_shares2]; 
    end
    
    Delta = [small_shares,large_shares];
    Delta_f = [small_shares_f,large_shares_f];
    
    DATA(35,itparam) = std(Delta(:));
    DATA(36,itparam) = std(Delta_f(:));
    
    DATA(3:38,itparam) = DATA(1:36,itparam);
    DATA(1:2,itparam) = [alpha_v;alpha_u];
    
%     % Compute Pareto-Shares
%     
%     for t=1:T
%         [row_s,col_s] = size(save_share_small_S{t});
%         [row_l,col_l] = size(save_share_large_S{t});
%         Rank_small = repmat(log((1:row_s)-1/2)',1,col_s);
%         Rank_large = repmat(log((1:row_l)-1/2)',1,col_l);
%         index_small = save_share_small_S{t}>0; 
%         index_large = save_share_large_S{t}>0;
%         thresh_small = prctile(save_share_small_S{t}(index_small),25);
%         thresh_large = prctile(save_share_large_S{t}(index_large),25);
%         index_small = save_share_small_S{t}>0 & sum(index_small)>10 & save_share_small_S{t} > thresh_small;
%         index_large = save_share_large_S{t}>0 & sum(index_large)>10 & save_share_large_S{t} > thresh_large;
%         Share_small = log(save_share_small_S{t}(index_small));
%         Share_large = log(save_share_large_S{t}(index_large));
%         Rank = [Rank_small(index_small);Rank_large(index_large)];
%         Share = [Share_small;Share_large];
%         stats = regstats(Rank,Share,'linear',{'beta','covb'});
%         Pareto(t,itparam) = stats.beta(2);       
%     end
%     DATA(33:end,itparam) = Pareto(:,itparam);
%     hold on
%     scatter(Rank,Share)
%     plot(Rank,-theta/4*Rank)
end

fname = sprintf('Results/Moments/30DynMom_rel_stand%d.csv',S);
TTT = cell2table(num2cell(DATA));
writetable(TTT,fname);

