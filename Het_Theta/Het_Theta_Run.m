% This code performs the estimation routine as outlined in the appendix,
% steps 1-3.

%% Housekeeping
clear;
close all;
clc;

tstart0 = tic;

addpath('Data');
addpath('Auxiliary Functions');

%% Load results with homogeneous theta
load('estimated_vmu_nopareto','sigma','theta','F','tau','muT','sigmaT','Y','YF','LF') % XXX Need to input correct file, make sure sigma=5

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

% Labor
L0 = 100;

% Markups
vMU = 1;

% Competition Structure
BER = 1;

% Pareto
paretonb = 75;

% Get CD-shares and Pareto shape from data 
cdtheta_init=csvread('cdshares_and_theta_v3.csv');            % Cobb-Douglas shares and Pareto Shape from external data source.

S_multiple = 4;                                         
S_init = length(cdtheta_init);                                % Number of sectors in data
S = S_init*S_multiple;                                   % Number of sectors used (see footnote 56)

% Get rid of zeros in Pareto shape
kappa_data=cdtheta_init(:,2);
kappa_data(kappa_data==0)=mean(kappa_data(kappa_data>0));
cdtheta_init(:,2)=kappa_data;

% Assign CD-shares (alpha) and Pareto shapes (kappa) across sectors 
doublevec = cdtheta_init;                      

for iloop = 1:S_multiple-1;
    doublevec = [doublevec;[cdtheta_init(randperm(S_init),1),cdtheta_init(randperm(S_init),2)]]; %#ok<AGROW>
end
doublevec(:,1) = doublevec(:,1)/S_multiple;

load('datamoments_56.csv')
datamoments=datamoments_56;
datamoments(55)=-datamoments_56(55);

% Split sectors into large and small based on CD-share.
ALPHA = doublevec(:,1);
split_param = 1.25;
small = (ALPHA<split_param/S);                            
Nsmall = sum(small);                           
Nlarge = S-Nsmall;

% Make random draws
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

% Given mu_T and sigma_T draw sectoral productivity T_z for each sector z 
RT=exp(muT+sigmaT*rtdraws);
RTS=RT(small);
RTL=RT(~small);

                              
%% Introduce Theta Heterogeneity

kappa_frenchdata=1.015;  
kappa_model_sim=1.096;


% CHOICE 1: Read theta from the datas

choice = 1;
correc = 1;

THETA = doublevec(:,2)*(sigma-1)*correc;         %% Translate kappa into theta

% Prepare inputs for loops
VEC_CD_THETA(:,1) = ALPHA;
VEC_CD_THETA(:,2) = THETA;

THETAS = THETA(small)';
THETAL = THETA(~small)';

% Draw productivities phi 
ZHS=(UHS./(ones(MH_small,1)*RTS)).^(-1./(ones(MH_small,1)*THETAS));
ZFS=UFS.^(-1./(ones(MF_small,1)*THETAS));

ZHL=(UHL./(ones(MH_large,1)*RTL)).^(-1./(ones(MH_large,1)*THETAL));
ZFL=UFL.^(-1./(ones(MF_large,1)*THETAL));
    
% Run loops (with GE of homogeneous theta as starting point)
[Y_het,YF_het,~,moments_het,~,~]=GEreplication_thetas(sigma,F,tau,RTS,RTL,ZHS,ZFS,ZHL,ZFL,w,wF,L0,Y,YF,small,vMU,BER,paretonb,VEC_CD_THETA); %% NEED TO WRITE THIS FUNCTION


%[~,moments_het]=PEmoments_thetas_v2(sigma,F,tau,RT,ZHS,ZFS,ZHL,ZFL,w,wF,Y_het,YF_het,L0,vMU,75,VEC_CD_THETA);  %% WHAT MOMENTS DO I ACTUALLY NEED?
    
name=strcat('Results/heterogenoustheta_',num2str(choice));
save(name,'moments_het','Y_het','YF_het','ALPHA','THETA','correc');


% CHOICE 2: Scaling theta down to get correct mean

pareto_shape(1)=moments_het(4);

choice = 2;
correc = 1/1.066*(datamoments(55)-0.009); % Need to understand and need to check datamoments(55)

THETA = doublevec(:,2)*(sigma-1)*correc; 

% Prepare inputs for loops
VEC_CD_THETA(:,1) = ALPHA;
VEC_CD_THETA(:,2) = THETA;

THETAS = THETA(small)';
THETAL = THETA(~small)';

% Draw productivities phi 
ZHS=(UHS./(ones(MH_small,1)*RTS)).^(-1./(ones(MH_small,1)*THETAS));
ZFS=UFS.^(-1./(ones(MF_small,1)*THETAS));

ZHL=(UHL./(ones(MH_large,1)*RTL)).^(-1./(ones(MH_large,1)*THETAL));
ZFL=UFL.^(-1./(ones(MF_large,1)*THETAL));

% Run loops (with GE of homogeneous theta as starting point)
[Y_het,YF_het,~,moments_het,~,~]=GEreplication_thetas(sigma,F,tau,RTS,RTL,ZHS,ZFS,ZHL,ZFL,w,wF,L0,Y,YF,small,vMU,BER,paretonb,VEC_CD_THETA); %% NEED TO WRITE THIS FUNCTION

%[~,moments_het]=PEmoments_thetas_v2(sigma,F,tau,RT,ZHS,ZFS,ZHL,ZFL,w,wF,Y_het,YF_het,L0,vMU,75,VEC_CD_THETA); %#ok<ASGLU> %% WHAT MOMENTS DO I ACTUALLY NEED?

name=strcat('Results/heterogenoustheta_',num2str(choice));
save(name,'moments_het','Y_het','YF_het','ALPHA','THETA','correc')


% CHOICE 3: Get right top market share

choice = 3;
correc = 1/pareto_shape(1)*1.22; % Need to understand and need to check pareto_shape(1)

THETA = doublevec(:,2)*(sigma-1)*correc; 

% Prepare inputs for loops
VEC_CD_THETA(:,1) = ALPHA;
VEC_CD_THETA(:,2) = THETA;

THETAS = THETA(small)';
THETAL = THETA(~small)';

% Draw productivities phi 
ZHS=(UHS./(ones(MH_small,1)*RTS)).^(-1./(ones(MH_small,1)*THETAS));
ZFS=UFS.^(-1./(ones(MF_small,1)*THETAS));

ZHL=(UHL./(ones(MH_large,1)*RTL)).^(-1./(ones(MH_large,1)*THETAL));
ZFL=UFL.^(-1./(ones(MF_large,1)*THETAL));

% Run loops (with GE of homogeneous theta as starting point)
[Y_het,YF_het,~,moments_het,~,~]=GEreplication_thetas(sigma,F,tau,RTS,RTL,ZHS,ZFS,ZHL,ZFL,w,wF,L0,Y,YF,small,vMU,BER,paretonb,VEC_CD_THETA); %% NEED TO WRITE THIS FUNCTION

%[~,moments_het]=PEmoments_thetas_v2(sigma,F,tau,RT,ZHS,ZFS,ZHL,ZFL,w,wF,Y_het,YF_het,L0,vMU,75,VEC_CD_THETA); %% WHAT MOMENTS DO I ACTUALLY NEED?

name=strcat('Results/heterogenoustheta_',num2str(choice));
save(name,'moments_het','Y_het','YF_het','ALPHA','THETA','correc')

%% NO Theta Heterogeneity

% CHOICE 1: Baseline

choice = 1;
correc = 1;

THETA = ones(S,1)*theta*correc;

VEC_CD_THETA(:,1) = ALPHA;
VEC_CD_THETA(:,2) = THETA;

THETAS = THETA(small)';
THETAL = THETA(~small)';

% Draw productivities phi 
ZHS=(UHS./(ones(MH_small,1)*RTS)).^(-1./(ones(MH_small,1)*THETAS));
ZFS=UFS.^(-1./(ones(MF_small,1)*THETAS));

ZHL=(UHL./(ones(MH_large,1)*RTL)).^(-1./(ones(MH_large,1)*THETAL));
ZFL=UFL.^(-1./(ones(MF_large,1)*THETAL));

% Run loops (with GE of homogeneous theta as starting point)
[Y_hom,YF_hom,~,moments_hom,~,~]=GEreplication_thetas(sigma,F,tau,RTS,RTL,ZHS,ZFS,ZHL,ZFL,w,wF,L0,Y,YF,small,vMU,BER,paretonb,VEC_CD_THETA); %% NEED TO WRITE THIS FUNCTION

%[~,moments_hom]=PEmoments_thetas_v2(sigma,F,tau,RT,ZHS,ZFS,ZHL,ZFL,w,wF,Y_hom,YF_hom,L0,vMU,75,VEC_CD_THETA); %% WHAT MOMENTS DO I ACTUALLY NEED?

name=strcat('Results/homogeneoustheta_',num2str(choice));
save(name,'moments_hom','Y_hom','YF_hom','ALPHA','THETA','correc')




% CHOICE 2: Get correct mean

choice = 2;
correc=1/kappa_model_sim*(kappa_frenchdata-0.0003);

THETA=ones(S,1)*theta*correc;

VEC_CD_THETA(:,1) = ALPHA;
VEC_CD_THETA(:,2) = THETA;

THETAS = THETA(small)';
THETAL = THETA(~small)';

% Draw productivities phi 
ZHS=(UHS./(ones(MH_small,1)*RTS)).^(-1./(ones(MH_small,1)*THETAS));
ZFS=UFS.^(-1./(ones(MF_small,1)*THETAS));

ZHL=(UHL./(ones(MH_large,1)*RTL)).^(-1./(ones(MH_large,1)*THETAL));
ZFL=UFL.^(-1./(ones(MF_large,1)*THETAL));

% Run loops (with GE of homogeneous theta as starting point)
[Y_hom,YF_hom,~,moments_hom,~,~]=GEreplication_thetas(sigma,F,tau,RTS,RTL,ZHS,ZFS,ZHL,ZFL,w,wF,L0,Y,YF,small,vMU,BER,paretonb,VEC_CD_THETA); %% NEED TO WRITE THIS FUNCTION

%[~,moments_hom]=PEmoments_thetas_v2(sigma,F,tau,RT,ZHS,ZFS,ZHL,ZFL,w,wF,Y_hom,YF_hom,L0,vMU,75,VEC_CD_THETA); %% WHAT MOMENTS DO I ACTUALLY NEED?

name=strcat('Results/homogeneoustheta_',num2str(choice));
save(name,'moments_hom','Y_hom','YF_hom','ALPHA','THETA','correc')

                                                 