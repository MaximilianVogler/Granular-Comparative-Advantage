% This code performs the estimation routine as outlined in the appendix,
% steps 1-3.

%% Housekeeping
clear all;
close all;
clc;

tstart0 = tic;

addpath('Data');
addpath('Auxiliary Functions');

%% Setup of model

% Initialize random number generator and set seed
aseed = 1;
rng(aseed);

% Relative size of economies
k = 1.75;

% Parameter governing the number of firms in each sector
M = 350;    

% Productivities for added firms
epsilon = 1e-10;
    
% Wages
wR = 1.13;      
w = 1;          
wF = 1/wR; 

% Markups
vMU = 1;

% Competition Structure
BER = 1;

% Compute number of sectors 
cdshares_init = csvread('cdshare.csv');             % Cobb-Douglas shares from external data source.

S_multiple = 4;                                         
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

% Get number of draws
MH = round(M*S*ALPHA);          
MF = round(k*M*S*ALPHA);        

% Fix sigma
SIGMA = 5;

% Labor Normalization
L0 = 100;

%% Grids
size_grid = 100;

%Sparse grid 1

% muT_vec = linspace(-0.5,0.5,size_grid);
% sigmaT_vec = linspace(1.0,1.8,size_grid);
% tau_vec = linspace(1.1,2.2,size_grid);
% kappa_vec = linspace(0.05,0.3,size_grid);
% f_vec = linspace(.5,6,size_grid)*4.93*.43*10^(-5);
 

%Sparse grid 2
% muT_vec = linspace(-0.2,0.8,size_grid);
% sigmaT_vec = linspace(0.7,1.7,size_grid);
% tau_vec = linspace(1.1,3.0,size_grid);
% kappa_vec = linspace(0.05,0.2,size_grid);
% f_vec = linspace(0.45,3.53,size_grid)*4.93*.43*10^(-5);
 

%Fine grid 3
muT_vec = linspace(-0.05,0.5,size_grid);
sigmaT_vec = linspace(0.5,1.4,size_grid);
tau_vec = linspace(1.3,3.2,size_grid);
kappa_vec = linspace(0.08,0.175,size_grid);
f_vec = linspace(0.35,4.2,size_grid)*4.93*.43*10^(-5);


%Fine grid 4
% muT_vec = linspace(-0.03,0.17,size_grid);
% sigmaT_vec = linspace(1.36,1.51,size_grid);
% tau_vec = linspace(1.31,1.39,size_grid);
% kappa_vec = linspace(1.05,1.10,size_grid);
% f_vec = linspace(1.99,2.81,size_grid)*4.93*.43*10^(-5);


%Fine grid 5
% muT_vec = linspace(0.06,0.16,size_grid);
% sigmaT_vec = linspace(1.37,1.45,size_grid);
% tau_vec = linspace(1.32,1.37,size_grid);
% kappa_vec = linspace(1.06,1.09,size_grid);
% f_vec = linspace(2.21,2.79,size_grid)*4.93*.43*10^(-5);

%Fine grid 6
% muT_vec=linspace(0.06,0.18,size_grid);
% sigmaT_vec=linspace(1.3025,1.4581,size_grid);
% tau_vec=linspace(1.3172,1.3883,size_grid);
% kappa_vec=linspace(1.0511,1.1211,size_grid);
% f_vec=linspace(1.95,3.82,size_grid)*4.93*.43*10^(-5);

% Set up Halton grid
Nbparam = 5;                  
NumGrid = 20000;                                        % Number of estimation points

p = haltonset(Nbparam);                                 % Set up Halton sequence that is choose which points you actually pick for estimation
p0 = net(p,NumGrid);
index = floor(p0*size_grid)+1;                          % This index tells you which points you actually choose. It is a NumGrid x Nbparam matrix, which tells you for each estimation point which 
                                                        % values of the parameters you choose.
                              
%% Estimation

% Initialize moment and parameter array (the only things we want from the estimation so that we know which parameter combination corresponds to which set of moments)
Momarray = zeros(NumGrid,15); 
Paramarray = zeros(NumGrid,5);

% Get parameter vectors with the parameter values chosen by the Halton sequence
muT_grid = muT_vec(index(:,1));
sigmaT_grid = sigmaT_vec(index(:,2));
tau_grid = tau_vec(index(:,3));
kappa_grid = kappa_vec(index(:,4));
f_grid = f_vec(index(:,5));


parfor i = 1:NumGrid
    
    muT=muT_grid(i);
    sigmaT=sigmaT_grid(i);
    tau=tau_grid(i);
    kappa=kappa_grid(i);
    f=f_grid(i);
    
    theta=(SIGMA-1)*kappa;
    sigma=SIGMA;
    
    Paramarray(i,:)=[muT sigmaT tau kappa f];
    
    F=f/sigma;

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
     
%     PHI = [phiH;phiF];
  
    % Split into large and small sectors (we can actually keep the old
    % splitting rule)
%     split_param = 500;
%     split_param2 = 10000;
%     PHI1 = PHI;
%     PHI1(PHI~=epsilon) = 1;
%     small = (sum(PHI1)<split_param);                         
%     Nsmall = sum(small);                           
%     Nlarge = S-Nsmall;
    
%     ZHS = PHI(1:max(sum(indexH(:,small))),small);
%     ZFS = PHI(max(MH)+1:max(MH)+max(sum(indexF(:,small))),small);
%     ZHL = PHI(1:max(MH),~small);
%     ZFL = PHI(max(MH)+1:end,~small);
    
    ZHS = phiH(1:max(sum(indexH(:,small))),small);
    ZFS = phiF(1:max(sum(indexF(:,small))),small);
    ZHL = phiH(:,~small);
    ZFL = phiF(:,~small);

    RTS = exp(mu_H(small));
    RTL = exp(mu_H(~small));
    
    % Guess home and foreign output Y
    Y0=123;
    YF0=2*Y0;
    
    tstart=tic
    
    % Run loops to solve model (step 3 of estimation procedure)
    [iter,Y,YF,LF,KHH,TOP1,TOP3,XS,YXS,LAMBDAHVEC,LAMBDAFVEC]=GEreplication_vectorized(sigma,theta,F,tau,ALPHA,RTS,RTL,ZHS,ZFS,ZHL,ZFL,w,wF,L0,Y0,YF0,small,vMU,BER,0);
    
    % Compute 15 target moments
    [Momarray(i,:)]=Moments(KHH,TOP1,TOP3,LAMBDAHVEC,LAMBDAFVEC,XS,YXS,ALPHA,Y,YF);   
    
    time=toc(tstart)

end

toc(tstart0)

Paramarray(:,5)=Paramarray(:,5)/(4.93*.43*10^(-5));

fprintf('Code has finished running')

save('Results/estimation_seed1_grid3','Momarray','Paramarray')                                                 