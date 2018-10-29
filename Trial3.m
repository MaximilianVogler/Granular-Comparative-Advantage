% Test of GE loop

clear all;
close all;
clc;

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

% Fix sigma
SIGMA = 5;

% Labor Normalization
L0 = 100;

% Input parameters
theta = 4.306666666666667;
F = 0.945561052525253*10^(-5);
tau = 1.340707070707071;

muT = 0.136767676767677;
sigmaT = 1.421717171717172;

RT=exp(muT+sigmaT*rtdraws);
RTS=RT(small);
RTL=RT(~small);

ZHS=(UHS./(repmat(RTS,MH_small,1))).^(-1/theta);
ZFS=UFS.^(-1/theta);
ZHL=(UHL./(repmat(RTL,MH_large,1))).^(-1/theta);                        
ZFL=UFL.^(-1/theta);

tstart=tic

BER = 1;

Y0 = 123;
YF0 = 2*Y0;



L = L0;
diff=1;
tol=1e-3;
iter=0;

tic
%while diff>tol && iter<51
    % Solve PE
    [K,KF,~,~,LAMBDA,LAMBDAF,MU,MUF,KHH,TOP1,TOP3,XS,YXS,LAMBDAHVEC,LAMBDAFVEC]=PEreplication_vectorized(SIGMA,theta,F,tau,ALPHA,RTS,RTL,ZHS,ZFS,ZHL,ZFL,w,wF,Y0,YF0,small,vMU,BER);
    
    % Set up linear system composed of (A10) and (A14)
    MAT=[(1-LAMBDA)*(MU-1)/MU-1 LAMBDAF*(MUF-1)/MUF;LAMBDA -LAMBDAF];
    VEC=[-w*L+LAMBDAF*wF*F*KF+(1-LAMBDA)*w*F*K;LAMBDA*w*F*K-LAMBDAF*wF*F*KF];
     SOL=MAT\VEC;
%     %SOL=MAT^(-1)*VEC;
% 
%     Y=SOL(1);
%     YF=SOL(2);
% 
%     diff=abs(Y-Y0)+abs(YF-YF0);
%     iter=iter+1;
%     if iter<=5
%         step=1;
%     else
%         step=0.5;
%     end
% 
%     Y0=Y0+step*(Y-Y0);
%     YF0=YF0+step*(YF-YF0);
% end
% 
Y=Y0;
YF=YF0;
 [Momarray]=Moments(KHH,TOP1,TOP3,LAMBDAHVEC,LAMBDAFVEC,XS,YXS,ALPHA,Y,YF); 
% 
% time=toc(tstart)