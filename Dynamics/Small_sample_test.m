% This code compute 30 target moments on a (rho,nu) grid.

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
scale = 3/4; %1 %3/4 %21
rep_size = 10;

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
rho_vec = 0:0.05:0.1;
nu_vec = 0.04:0.005:0.06;
[rho_mat,nu_mat] = meshgrid(rho_vec,nu_vec); 
[row,col] = size(rho_mat);
num_param = row*col;

DATA = zeros(6,num_param);

ZHS_start = ZHS;
ZHL_start = ZHL;
ZFS_start = ZFS;
ZFL_start = ZFL;

% Iterate over this parameter grid
for itparam = 1:num_param
    
    % Set the correct parameters for this iteration
    rho = rho_mat(itparam);
    nu = nu_mat(itparam);
    mu = -theta*nu^2/2;
    
    mom = zeros(2,rep_size);
    
    aseed = 1;
    rng(aseed);
    
    for rep = 1:rep_size
        % Initialize matrices of interest
        LAMBDAFVEC_t = zeros(R_length,S);
        DVEC_t = zeros(R_length,S);
        XVEC_t = zeros(R_length,S);
        PHIFVEC_t = zeros(R_length,S);

        su = sqrt(rho*nu^2);
        sv = sqrt(nu^2*(1-rho));

        % Re-initialize matrices, so every iteration starts at the same
        % productivity

        ZHS = ZHS_start;
        ZHL = ZHL_start;
        ZFS = ZFS_start;
        ZFL = ZFL_start;

        counter = 1;

        for t=1:T

            % Determine productivity draws
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
                save_share_small{t} = DSHM_small;
                save_share_large{t} = DSHM_large;
                disp(['Loop ',num2str(t),' is finished']);
            end

        end
        disp(['Iteration ',num2str(itparam),' is finished']);

        % Moment 22
        k = 10;
        Dep_f = zeros(0);
        Indep_f = zeros(0);
        for t=1:T-k
            Dep2 = (XVEC_t(t+k,:)./DVEC_t(t+k,:))-(XVEC_t(t,:)./DVEC_t(t,:));
            Indep2 = XVEC_t(t,:)./DVEC_t(t,:);
            ind2 = (XVEC_t(t+k,:)>0) & (DVEC_t(t+k,:)>0) & (XVEC_t(t,:)>0) & (DVEC_t(t,:)>0);
            Dep2 = Dep2(ind2);
            Indep2 = Indep2(ind2);
            Dep_f = [Dep_f,Dep2];
            Indep_f = [Indep_f,Indep2];
        end
        
        % Regression for moments 22
        mom(1,rep) = Indep_f*Indep_f'\(Indep_f*Dep_f');
    
    
        % Moment 24
        small_shares_dep = zeros(T-k,0);
        small_shares_indep = zeros(T-k,0);
        counter = counter + 1;
        for z=1:Nsmall
            clearvars small_shares1_dep small_shares1_indep
            % Indicator for being active in all periods
            dim = size(DSHM_small);
            idx_z_active = true(dim(1),1);
            for t=1:T
               idx_z_active = idx_z_active & (save_share_small{t}(:,z)>0);
            end
            % Indicator for belonging to the top 20% among all active firms in
            % year 1
            for t=1:T-k
               int_dep = save_share_small{t+k}(idx_z_active,z) - save_share_small{t}(idx_z_active,z);
               int_indep = save_share_small{t}(idx_z_active,z);
               small_shares1_dep(t,:) = int_dep';
               small_shares1_indep(t,:) = int_indep';
            end
            % This is the target matrix that contains a matrix with T rows and
            % where each column represents the development of one active firm
            % in any sector over time
            small_shares_dep = [small_shares_dep,small_shares1_dep];
            small_shares_indep = [small_shares_indep,small_shares1_indep];
        end

        large_shares_dep = zeros(T-k,0);
        large_shares_indep = zeros(T-k,0);
        for z=1:Nlarge
            clearvars large_shares1_dep large_shares2_dep large_shares1_indep large_shares2_indep 
            dim = size(DSHM_large);
            idx_z_active = true(dim(1),1);
            for t=1:T
               idx_z_active = idx_z_active & (save_share_large{t}(:,z)>0); 
            end
            for t=1:T-k
                int_dep = save_share_large{t+k}(idx_z_active,z) - save_share_large{t}(idx_z_active,z);
                int_indep = save_share_large{t}(idx_z_active,z);
                large_shares1_dep(t,:) = int_dep';
                large_shares1_indep(t,:) = int_indep';
            end
            large_shares_dep = [large_shares_dep,large_shares1_dep];
            large_shares_indep = [large_shares_indep,large_shares1_indep];
        end

        Dep = [small_shares_dep,large_shares_dep];
        Indep = [small_shares_indep,large_shares_indep];


        % Regression
        mom(2,rep) = Indep(:)'*Indep(:)\(Indep(:)'*Dep(:));
    end
    
    DATA(1:2,itparam) = [rho;nu];
    DATA(3,itparam) = mean(mom(1,:));
    DATA(4,itparam) = std(mom(1,:));
    DATA(5,itparam) = mean(mom(2,:));
    DATA(6,itparam) = std(mom(2,:));
    
end

fname = sprintf('Results/Moments/MomTest_%d.csv',S);
TTT = cell2table(num2cell(DATA));
writetable(TTT,fname);

