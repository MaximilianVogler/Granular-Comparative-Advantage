%% This code performs local minimization around the 20 best grid points as in step 6 of the estimation procedure.

%% Housekeeping

clear all;
close all;
clc;

addpath('Data');
addpath('Auxiliary Functions');
addpath('Results');

%% Load Results from Grid_Optimization and data moments

load('GridOptimization_seed1_grid6');
datamoments = csvread('Data_Moments.csv');

%% Set up Estimation as in Estimation.m

tic

% Initialize random number generator and set seed
aseed = 1;
rng(aseed);

% Relative size of economies
k = 1.75;

% Parameter governing the number of firms in each sector
M = 350;   

% Wages
wR = 1.13;      
w = 1;          
wF = 1/wR; 

% Markups
vMU = 1;

% Competition Structure
BER = 1;

% Fix sigma
SIGMA = 5;

% Compute number of sectors 
cdshares_init = csvread('cdshares_v3.csv');             % Cobb-Douglas shares from external data source.

S_multiple = 4;                                         
S_init = length(cdshares_init);                         % Number of sectors in data
S = S_init*S_multiple;                                  % Number of sectors used (see footnote 56)

% Assign CD-shares across sectors 
ALPHA = cdshares_init;      

for iloop = 1:S_multiple-1;
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

% Get number of draws
MH = round(M*S*ALPHA);          
MF = round(k*M*S*ALPHA); 

% Draw normal mean productivity of home firms
mu_H = muT+sigmaT*rtdraws;

% Draw log-normal productivities of home and foreign firms
phiH = ones(max(MH),S)*epsilon;
phiF = ones(max(MF),S)*epsilon;

indexH = false(max(MH),S);
indexF = false(max(MF),S);

rng(aseed);
for j = 1:S
    phiH(1:MH(j),j) = exp(mu_H(j)+theta*randn(MH(j),1));
    indexH(1:MH(j),j) = true(MH(j),1);
    phiF(1:MF(j),j) = exp(theta*randn(MF(j),1));
    indexF(1:MF(j),j) = true(MF(j),1);
end

ZHS = phiH(1:max(sum(indexH(:,small))),small);
ZFS = phiF(1:max(sum(indexF(:,small))),small);
ZHL = phiH(:,~small);
ZFL = phiF(:,~small);

RTS = exp(mu_H(small));
RTL = exp(mu_H(~small));
%% Set up weight matrix W

ref_weight = datamoments;
%ref_weight(1:2) = 1;                                        % CAREFUL, do we really need this?
ref_weight(14:15) = datamoments(12:13);
std_indices = [1,3,5,7,9];
ref_weight(std_indices) = datamoments(std_indices)/3;
W = diag((1./(ref_weight.^2)));

lossfun=@(x) Loss_Function(x(1),x(2),x(3),x(4),x(5),SIGMA,ZHS,ZHL,ZFS,ZFL,RTS,RTL,rtdraws,small,MH_small,MH_large,ALPHA,w,wF,vMU,BER,datamoments,W);

options=optimset('TolFun',10^(-5),'TolX',10^(-3));

NN=size(paramsP,1);

paramsP_init=zeros(NN,6);
paramsP_localmin=zeros(NN,5);
loss_value=zeros(NN,1);

parfor i=1:NN
    
    init_point=paramsP(i,1:5);
    [Xbest,FVAL]=fminsearch(lossfun,init_point,options);
    paramsP_localmin(i,:)=Xbest;
    loss_value(i)=FVAL;
    paramsP_init(i,:)=paramsP(i,1:6);
    
end

time=toc;

save('local_min_seed1_grid6','paramsP_localmin','paramsP_init','loss_value','time')
