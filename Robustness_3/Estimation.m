% This code performs the estimation routine as outlined in the appendix,
% steps 1-3.

%% Housekeeping
clear all;
close all;
clc;

tstart0 = tic

addpath('Data');
addpath('Auxiliary Functions');

%% Setup of model

% Initialize random number generator and set seed
aseed = 1;
rng(aseed);

% Number of shadow firms
MF_large=10000;
MF_small=1400;
MH_large=5000;
MH_small=700;

% Wages
wR = 1.13;      
w = 1;          
wF = 1/wR; 

% Markups
vMU = 1;

% Competition Structure
BER = 1;

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

UH0S = exprnd(1,MH_small,Nsmall);                         % Draw U of most productive small home shadow firm and spacings in each sector from exponential with mean 1                        
UHS = cumsum(UH0S);                                       % Cumulate to get U of each home firm in each sector (see footnote 57)

UH0L = exprnd(1,MH_large,Nlarge);
UHL = cumsum(UH0L);

% Fix sigma
sigma = 5;

% Labor Normalization
L0 = 100;

load('GE_Results');
muT=bestParams(1);
sigmaT=bestParams(2);
tau=bestParams(3);
theta=bestParams(4);
F=bestParams(5);
                              
%% Estimation

% Initialize moment and parameter array (the only things we want from the estimation so that we know which parameter combination corresponds to which set of moments)
Momarray = zeros(1,15); 

% Given mu_T and sigma_T draw sectoral productivity T_z for each sector z (step 1 of estimation procedure)
RT=exp(muT+sigmaT*rtdraws);
RTS=RT(small);
RTL=RT(~small);

% Draw productivities phi (step 2 of estimation procedure)
ZHS=(UHS./(repmat(RTS,MH_small,1))).^(-1/theta);
ZHL=(UHL./(repmat(RTL,MH_large,1))).^(-1/theta);   

underlineVarphiSmall = (1/MF_small)^(1/theta);
underlineVarphiLarge = (1/MF_large)^(1/theta);

k = 1/theta;
sigmaMatlabSmall = underlineVarphiSmall/theta;
sigmaMatlabLarge = underlineVarphiLarge/theta;
ParetoSmall = makedist('GeneralizedPareto','k',k,'sigma',sigmaMatlabSmall,'theta',underlineVarphiSmall);
ParetoLarge = makedist('GeneralizedPareto','k',k,'sigma',sigmaMatlabLarge,'theta',underlineVarphiLarge);

percentilesSmall = [1-(1/(2*MF_small)) : -(1/MF_small) : 1/(2*MF_small)];
percentilesLarge = [1-(1/(2*MF_large)) : -(1/MF_large) : 1/(2*MF_large)];

prodFSmall = icdf(ParetoSmall,percentilesSmall)';
prodFLarge = icdf(ParetoLarge,percentilesLarge)';

ZFS = repmat(prodFSmall,1,Nsmall);
ZFL = repmat(prodFLarge,1,Nlarge);


tstart = tic

% Run loops in PE
% [K,KF,~,~,LAMBDA,LAMBDAF,MU,MUF,KHH,TOP1,TOP3,XS,YXS,LAMBDAHVEC,LAMBDAFVEC,~,~] = PEreplication_vectorized(sigma,theta,F,tau,ALPHA,RTS,RTL,ZHS,ZFS,ZHL,ZFL,w,wF,Y,YF,small,vMU,BER,75,0);

% Run loops in GE
% Guess home and foreign output Y
 Y0=123;
 YF0=2*Y0;
[iter,Y,YF,LF,KHH,TOP1,TOP3,XS,YXS,LAMBDAHVEC,LAMBDAFVEC]=GEreplication_vectorized(sigma,theta,F,tau,ALPHA,RTS,RTL,ZHS,ZFS,ZHL,ZFL,w,wF,L0,Y0,YF0,small,vMU,BER,0);
    


% Compute 15 target moments
[Momarray(:)]=Moments(KHH,TOP1,TOP3,LAMBDAHVEC,LAMBDAFVEC,XS,YXS,ALPHA,Y,YF);   

time = toc(tstart)


fprintf('Code has finished running')

% save('Results/PEMoments','Momarray','bestParams')  
save('Results/GEMoments','Momarray','bestParams','LF','Y','YF','wF')  