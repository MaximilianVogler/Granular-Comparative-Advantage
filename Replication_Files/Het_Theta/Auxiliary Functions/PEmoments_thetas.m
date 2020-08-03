function [mom] = PEmoments_thetas(sigma,F,tau,theta,muT,sigmaT,w,wF,Y,YF,vMU,BER,paretonb,correc,het)
% This codes computes all relevant moments for table 3 in PE with more
% sectors.

%% Setup of model

% Initialize random number generator and set seed
aseed = 1;
rng(aseed);

% Number of shadow firms
MF_large=10000;
MF_small=1400;
MH_large=5000;
MH_small=700;

% Get CD-shares and Pareto shape from data 
cdtheta_init=csvread('cdshares_and_theta.csv');            % Cobb-Douglas shares and Pareto Shape from external data source.

S_multiple = 4;                                         
S_init = length(cdtheta_init);                                % Number of sectors in data
S = S_init*S_multiple;                                        % Number of sectors used (see footnote 56)

% Get rid of zeros in Pareto shape
kappa_data=cdtheta_init(:,2);
kappa_data(kappa_data==0)=mean(kappa_data(kappa_data>0));
cdtheta_init(:,2)=kappa_data;

% Assign CD-shares (alpha) and Pareto shapes (kappa) across sectors 
doublevec = cdtheta_init;                      

for iloop = 1:S_multiple-1;
    doublevec = [doublevec;[cdtheta_init(randperm(S_init),1),cdtheta_init(randperm(S_init),2)]]; 
end
doublevec(:,1) = doublevec(:,1)/S_multiple;

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

% Determine correct theta vector
if het == 0
    THETA = ones(S,1)*theta*correc;    
elseif het == 1
    THETA = doublevec(:,2)*(sigma-1)*correc;
else 
    error('Het has to be either 1 or 2.')
end
    

VEC_CD_THETA(:,1) = ALPHA;
VEC_CD_THETA(:,2) = THETA;

THETAS = THETA(small)';
THETAL = THETA(~small)';

% Draw productivities phi 
ZHS=(UHS./(ones(MH_small,1)*RTS)).^(-1./(ones(MH_small,1)*THETAS));
ZFS=UFS.^(-1./(ones(MF_small,1)*THETAS));

ZHL=(UHL./(ones(MH_large,1)*RTL)).^(-1./(ones(MH_large,1)*THETAL));
ZFL=UFL.^(-1./(ones(MF_large,1)*THETAL));

%% Run PE model

[~,~,~,~,~,~,~,~,~,mom] = PEreplication_thetas(sigma,F,tau,RTS,RTL,ZHS,ZFS,ZHL,ZFL,w,wF,Y,YF,small,vMU,BER,paretonb,VEC_CD_THETA);


end