function[ALPHA,PHIFVEC,LAMBDAFVEC,KVEC,KFVEC,KHH,KFH,MH,RT,AdditionalMom] = PEmoments(sigma,F,tau,theta,muT,sigmaT,w,wF,Y,YF,vMU,BER,paretonb,scale)

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

% Random draws for loop
rng(aseed);                                               % Reset random number generator for consistency with old code
rtdraws = randn(1,S);

% Get number of draws
MH = round(M*S*ALPHA);          
MF = round(k*M*S*ALPHA);        

% Labor Normalization
L0 = 100;

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
RT = exp(mu_H);

%% Run PE model

[~,~,KVEC,KFVEC,~,~,~,~,~,~,KHH,KFH,~,~,~,~,~,LAMBDAFVEC,PHIFVEC,AdditionalMom] = PEreplication_vectorized(sigma,theta,F,tau,ALPHA,RTS,RTL,ZHS,ZFS,ZHL,ZFL,w,wF,Y,YF,small,vMU,BER,paretonb,1);

end