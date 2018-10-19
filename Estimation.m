%% Housekeeping
clear all;
close all;
clc;

tstart0 = tic

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
vMU=1;

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
rtdraws = randn(1,S);

UH0S = exprnd(1,MH_small,Nsmall);                         % Draw U of most productive small home shadow firm and spacings in each sector from exponential with mean 1
UF0S = exprnd(1,MF_small,Nsmall);                        
UHS = cumsum(UH0S);                                       % Cumulate to get U of each home firm in each sector (see footnote 57)
UFS = cumsum(UF0S);

UH0L = exprnd(1,MH_large,Nlarge);
UF0L = exprnd(1,MF_large,Nlarge);
UHL = cumsum(UH0L);
UFL = cumsum(UF0L);

% Fix sigma
SIGMA = 5;

% Labor Normalization
L0 = 100;

%% Grids
size_grid = 100;

%Sparse grid 1

% muT_vec = linspace(-0.5,0.5,size_grid);
% sigmaT_vec = linspace(1,1.8,size_grid);
% tau_vec = linspace(1.1,2.2,size_grid);
% kappa_vec = linspace(.95,1.2,size_grid);
% f_vec = linspace(.5,6,size_grid)*4.93*.43*10^(-5);
%  

%Sparse grid 2
% muT_vec = linspace(-0.16,0.3,size_grid);
% sigmaT_vec = linspace(1.2,1.7,size_grid);
% tau_vec = linspace(1.2,1.6,size_grid);
% kappa_vec = linspace(1.01,1.13,size_grid);
% f_vec = linspace(1.1,4.8,size_grid)*4.93*.43*10^(-5);
 

%Fine grid 3
% muT_vec = linspace(-0.06,0.22,size_grid);
% sigmaT_vec = linspace(1.29,1.55,size_grid);
% tau_vec = linspace(1.3,1.4,size_grid);
% kappa_vec = linspace(1.05,1.11,size_grid);
% f_vec = linspace(1.55,3.95,size_grid)*4.93*.43*10^(-5);


%Fine grid 4
% muT_vec = linspace(-0.03,0.17,size_grid);
% sigmaT_vec = linspace(1.36,1.51,size_grid);
% tau_vec = linspace(1.31,1.39,size_grid);
% kappa_vec = linspace(1.05,1.10,size_grid);
% f_vec = linspace(1.99,2.81,size_grid)*4.93*.43*10^(-5);


%Fine grid 5
muT_vec = linspace(0.06,0.16,size_grid);
sigmaT_vec = linspace(1.37,1.45,size_grid);
tau_vec = linspace(1.32,1.37,size_grid);
kappa_vec = linspace(1.06,1.09,size_grid);
f_vec = linspace(2.21,2.79,size_grid)*4.93*.43*10^(-5);

% Set up Halton grid
Nbparam = 5;                  
NumGrid = 20000;                                        % Number of estimation points

p = haltonset(Nbparam);                                 % Set up Halton sequence that is choose which points you actually pick for estimation
p0 = net(p,NumGrid);
index = floor(p0*size_grid)+1;   % This index tells you which points you actually choose. It is a NumGrid x Nbparam matrix, which tells you for each estimation point which 
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
    
    % What is going on here?
    F=f/sigma;
    
    % Given mu_T and sigma_T draw sectoral productivity T_z for each sector z (step 1 of estimation procedure)
    RT=exp(muT+sigmaT*rtdraws);
    RTS=RT(small);
    RTL=RT(~small);
    
    % Draw productivities phi (step 2 of estimation procedure)
    ZHS=(UHS./(ones(MH_small,1)*RTS)).^(-1/theta);
    ZFS=UFS.^(-1/theta);
    ZHL=(UHL./(ones(MH_large,1)*RTL)).^(-1/theta);
    ZFL=UFL.^(-1/theta);
    
    % Guess home and foreign output Y
    Y0=123;
    YF0=2*Y0;
    
    tstart=tic
    
    % Run loops to solve model (step 3 of estimation procedure)
    [Y,YF,~]=GEreplication_vectorized(sigma,theta,F,tau,ALPHA,RTS,RTL,ZHS,ZFS,ZHL,ZFL,w,wF,L0,Y0,YF0,vMU,BER);
    
    % Compute 15 target moments
    [Momarray(i,:)]=PEmoments_vectorized_new_extra_v2(sigma,theta,F,tau,ALPHA,RTS,RTL,ZHS,ZFS,ZHL,ZFL,w,wF,Y,YF,L0,vMU,BER);
    
    time=toc(tstart)
end

toc(tstart0)

% Need to understand why we need this
Paramarray(:,5)=Paramarray(:,5)/(4.93*.43*10^(-5));

save('estimation_seed1_grid5','Momarray','Paramarray')                                                 