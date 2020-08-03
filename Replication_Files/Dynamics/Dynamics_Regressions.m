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

%% Setup of model

% Number of run
nRuns = 20;

% Results
table1_const = zeros(nRuns,4);
table1_top3 = zeros(nRuns,4);
table1_D = zeros(nRuns,2);

table2_const = zeros(nRuns,4);
table2_X = zeros(nRuns,3);
table2_top3 = zeros(nRuns,2);
table2_D = zeros(nRuns,4);

for n=1:nRuns
    
% Initialize random number generator and set seed
aseed = n;
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
% var_thresh = 1;

% Scale of model (i.e. number of sectors)
scale = 1; %1 %3/4 %21

% Compute number of sectors 
cdshares_init = csvread('cdshare.csv');             % Cobb-Douglas shares from external data source.

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
% alpha_u_vec = [0.05:0.0001:0.051];
% 	_vec = [0.034:0.0001:0.035];
alpha_u_vec = [0.0502];
alpha_v_vec = [0.0342];

[alpha_v_mat,alpha_u_mat] = meshgrid(alpha_v_vec,alpha_u_vec); 
% [row,col] = size(alpha_v_mat);
% num_param = row*col;
num_param = length(alpha_u_vec); 

DATA = zeros(32,num_param);

ZHS_start = ZHS;
ZHL_start = ZHL;
ZFS_start = ZFS;
ZFL_start = ZFL;

% Need to keep track of initial threshold only for home firms (stays the
% same for foreign firms)
VB_HS_start = VB_HS;
VB_HL_start = VB_HL;

RTS_start = RTS;
RTL_start = RTL;

% Need to keep track of initial AR(1) component
VZ_small_start = (log(RTS)-muT)/theta;
VZ_large_start = (log(RTL)-muT)/theta;

Pareto = zeros(T,num_param);

% Iterate over this parameter grid
for itparam = 1:num_param
    
    % Set the correct parameters for this iteration
%     alpha_v = alpha_v_mat(itparam);
%     alpha_u = alpha_u_mat(itparam);
    alpha_v = alpha_v_vec(itparam);
    alpha_u = alpha_u_vec(itparam);
    mu = -theta*alpha_u^2/2;
    rho_v = sqrt(1-(theta*alpha_v/sigmaT)^2);
    
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
    
    % Re-initialize thresholds (only for home firms)
    VB_HS = VB_HS_start;
    VB_HL = VB_HL_start;
    
    % Re-initialize AR(1)-components (only for home firms)
    VZ_HS = repmat(VZ_small_start,ZHS_length,1);
    VZ_HL = repmat(VZ_large_start,ZHL_length,1);
    
    % Re-initialize comparative advantage
    RTS = RTS_start;
    RTL = RTL_start;

    aseed = 1;
    rng(aseed);

    counter = 1;
    
    for t=1:T

        % Determine productivity draws
        
        % Sectoral Draws (only for home firms)
        epsv_HS = repmat(randn(1,Nsmall),ZHS_length,1);
        epsv_HL = repmat(randn(1,Nlarge),ZHL_length,1);

        % Firm-level Draws
        epsu_HS = randn(ZHS_length,Nsmall);
        epsu_HL = randn(ZHL_length,Nlarge);
        epsu_FS = randn(ZFS_length,Nsmall);
        epsu_FL = randn(ZFL_length,Nlarge);
        
        % Save previous threshold (needed for evolution)
        VB_HS_old = VB_HS;
        VB_HL_old = VB_HL;
        
        % Evolution of log productivity thresholds (only for home firms)
        VB_HS = VB_HS - (1-rho_v)*VZ_HS + alpha_v*epsv_HS;
        VB_HL = VB_HL - (1-rho_v)*VZ_HL + alpha_v*epsv_HL;
        
        % Keep track of AR(1) component
        VZ_HS = rho_v*VZ_HS + alpha_v*epsv_HS;
        VZ_HL = rho_v*VZ_HL + alpha_v*epsv_HL;      
        
        % Evolution of comparative advantage
        RTS = exp(log(RTS)+theta*(VB_HS(1,:)-VB_HS_old(1,:)));
        RTL = exp(log(RTL)+theta*(VB_HL(1,:)-VB_HL_old(1,:)));

        % Evolution of productivity draws
        ZHS(~inactive_HS) = exp(VB_HS(~inactive_HS)+abs(log(ZHS(~inactive_HS))-VB_HS_old(~inactive_HS)+mu+alpha_u*epsu_HS(~inactive_HS)));
        ZHL(~inactive_HL) = exp(VB_HL(~inactive_HL)+abs(log(ZHL(~inactive_HL))-VB_HL_old(~inactive_HL)+mu+alpha_u*epsu_HL(~inactive_HL)));
        ZFS(~inactive_FS) = exp(VB_FS(~inactive_FS)+abs(log(ZFS(~inactive_FS))-VB_FS(~inactive_FS)+mu+alpha_u*epsu_FS(~inactive_FS)));
        ZFL(~inactive_FL) = exp(VB_FL(~inactive_FL)+abs(log(ZFL(~inactive_FL))-VB_FL(~inactive_FL)+mu+alpha_u*epsu_FL(~inactive_FL)));

        % If the year is part of RECORD, record PE results
        if any(RECORD==t)
            [~,~,~,~,~,~,~,~,~,~,~,~,~,~,LAMBDAFVEC,PHIFVEC,~,~,XVEC,DVEC,DSHM_small,DSHM_large,DSHMS_small,DSHMS_large]=PEreplication_vectorized(sigma,theta,F,tau,ALPHA,RTS,RTL,ZHS,ZFS,ZHL,ZFL,w,wF,Y,YF,small,vMU,BER,0,0);
            LAMBDAFVEC_t(t,:) = LAMBDAFVEC;
            PHIFVEC_t(t,:) = PHIFVEC;
            XVEC_t(t,:) = XVEC; 
            DVEC_t(t,:) = DVEC;
            save_share_small{t} = DSHM_small;
            save_share_large{t} = DSHM_large;
            save_share_small_S{t} = DSHMS_small;
            save_share_large_S{t} = DSHMS_large;
            disp(['Loop ',num2str(t),' is finished']);
        end
    
    end
    disp(['Iteration ',num2str(itparam),' is finished']);
    

    for t=1:T
        TOP1_t(t,small) = save_share_small{t}(1,:);
        TOP3_t(t,small) = sum(save_share_small{t}(1:3,:));
        TOP1_t(t,~small) = save_share_large{t}(1,:);
        TOP3_t(t,~small) = sum(save_share_large{t}(1:3,:));
    end
    
    
end

%% Table 1
% Regression 1

XVEC_t=XVEC_t';
TOP3_t=TOP3_t';
DVEC_t=DVEC_t';

ind = (XVEC_t(:,8)>0) & (DVEC_t(:,8)>0);

Y_reg = log(XVEC_t(ind,8));
X_reg = [ones(sum(ind),1),TOP3_t(ind,8),log(DVEC_t(ind,8))];

res = (X_reg'*X_reg)\X_reg'*Y_reg;
table1_const(n,1) = res(1);
table1_top3(n,1) = res(2);
table1_D(n,1) = res(3);

% Regression 2
ind = (XVEC_t>0) & (DVEC_t>0);
Y_reg = log(XVEC_t(ind));
X_reg = [ones(sum(sum(ind)),1),TOP3_t(ind),log(DVEC_t(ind))];

res = (X_reg'*X_reg)\X_reg'*Y_reg;
table1_const(n,2) = res(1);
table1_top3(n,2) = res(2);
table1_D(n,2) = res(3);

% Regression 3
Y_reg = log(XVEC_t)-mean(log(XVEC_t),2);
ind = isfinite(Y_reg);
Y_reg = Y_reg(ind);
TOP_3_dm = TOP3_t-mean(TOP3_t,2);
X_reg = [ones(sum(sum(ind)),1),TOP_3_dm(ind)];

res = (X_reg'*X_reg)\X_reg'*Y_reg;
table1_const(n,3) = res(1);
table1_top3(n,3) = res(2);

% Regression 4
Y_reg = log(XVEC_t(:,2:end))-log(XVEC_t(:,1:end-1));
ind = isfinite(Y_reg);
Y_reg = Y_reg(ind);
TOP3_d = TOP3_t(:,2:end)-TOP3_t(:,1:end-1);
X_reg = [ones(sum(sum(ind)),1),TOP3_d(ind)];

res = (X_reg'*X_reg)\X_reg'*Y_reg;
table1_const(n,4) = res(1);
table1_top3(n,4) = res(2);

%% Table 2

% Regression 1
Y_reg = log(XVEC_t(:,11))-log(XVEC_t(:,1));

ind = isfinite(Y_reg) & isfinite(log(XVEC_t(:,1))) & isfinite(log(DVEC_t(:,1)));
X_reg = [ones(sum(ind),1),log(XVEC_t(ind,1)),log(DVEC_t(ind,1))];
Y_reg = Y_reg(ind);
res = (X_reg'*X_reg)\X_reg'*Y_reg;
table2_const(n,1) = res(1);
table2_X(n,1) = res(2);
table2_D(n,1) = res(3);

% Regression 2
Y_reg = log(XVEC_t(:,11))-log(XVEC_t(:,1));
ind = isfinite(Y_reg) & isfinite(log(XVEC_t(:,1))) & isfinite(log(DVEC_t(:,1)));
X_reg = [ones(sum(ind),1),TOP3_t(ind,1),log(DVEC_t(ind,1))];
Y_reg = Y_reg(ind);

res = (X_reg'*X_reg)\X_reg'*Y_reg;
table2_const(n,2) = res(1);
table2_top3(n,1) = res(2);
table2_D(n,2) = res(3);

% Regression 3
Y_reg = log(XVEC_t(:,11))-log(XVEC_t(:,1));
ind = isfinite(Y_reg) & isfinite(log(XVEC_t(:,1))) & isfinite(log(DVEC_t(:,1)));
X_reg = [ones(sum(ind),1),log(XVEC_t(ind,1)),TOP3_t(ind,1),log(DVEC_t(ind,1))];
Y_reg = Y_reg(ind);

res = (X_reg'*X_reg)\X_reg'*Y_reg;
table2_const(n,3) = res(1);
table2_X(n,2) = res(2);
table2_top3(n,2) = res(3);
table2_D(n,3) = res(4);

% Regression 4
Y_reg = log(XVEC_t(:,11))-log(XVEC_t(:,1));

ind = isfinite(Y_reg) & isfinite(log(XVEC_t(:,1))) & isfinite(log(DVEC_t(:,1)));
X_reg = [ones(sum(ind),1),TOP3_t(ind,1),log(DVEC_t(ind,1))];

endo = log(XVEC_t(ind,1));
Y_reg = Y_reg(ind);

beta = (X_reg'*X_reg)\X_reg' * endo;

endo_hat = X_reg * beta;


X_reg = [ones(sum(ind),1),endo_hat,log(DVEC_t(ind,1))];

res = (X_reg'*X_reg)\X_reg'*Y_reg;
table2_const(n,4) = res(1);
table2_X(n,3) = res(2);
table2_D(n,4) = res(3);

disp(['Run ',num2str(n),' is finished'])


end


table1_top3_sort = sort(table1_top3);
table1_D_sort = sort(table1_D);
table2_D_sort = sort(table2_D);
table2_X_sort = sort(table2_X);
table2_top3_sort = sort(table2_top3);

