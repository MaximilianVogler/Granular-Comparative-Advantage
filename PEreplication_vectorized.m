function [K,KF,LAMBDA,LAMBDAF,MU,MUF] = PEreplication_vectorized(sigma,theta,F,tau,ALPHAS,ALPHAL,RTS,RTL,ZHS,ZFS,ZHL,ZFL,w,wF,Y0,YF0,small,vMU,BER)
% Takes parameters, random draws, and initial guess as inputs. Outputs are
% K, Lambda, and aggregate markup Mu in both countries. 

S=size(ZHL,2)+size(ZHS,2); % fix number of sectors

% Initialize vectors that contain PE data
PHIHVEC=zeros(1,S);
MUHVEC=zeros(1,S);
PHIFVEC=zeros(1,S);
MUFVEC=zeros(1,S);
KVEC=zeros(1,S);
KFVEC=zeros(1,S);
LAMBDAHVEC=zeros(1,S);
LAMBDAFVEC=zeros(1,S); 

% Solve inner loops separately for small and large sectors (for speed)
[KVEC(small),KFVEC(small),PHIHVEC(small),PHIFVEC(small),MUHVEC(small),MUFVEC(small),LAMBDAHVEC(small),LAMBDAFVEC(small)]=...
    Inner_Loops(sigma,theta,F,tau,ALPHAS',RTS,ZHS,ZFS,w,wF,Y0,YF0,vMU,BER);
[KVEC(~small),KFVEC(~small),PHIHVEC(~small),PHIFVEC(~small),MUHVEC(~small),MUFVEC(~small),LAMBDAHVEC(~small),LAMBDAHVEC(~small)]=...
    Inner_Loops(sigma,theta,F,tau,ALPHAL',RTL,ZHL,ZFL,w,wF,Z0,YF0,vMU,BER);

% Take averages of vectors to obtain aggregates       CAREFUL!!! CORRECT???
K=mean(KVEC); % number of firms in home
KF=mean(KFVEC);
PHI=mean(PHIHVEC); % expected import share in home
PHIF=mean(PHIFVEC);
LAMBDA=mean(LAMBDAHVEC); % realized import share in home
LAMBDAF=mean(LAMBDAFVEC);
MU=1./sum(ALPHA.*MUHVEC); % aggregate markup in home
MUF=1./sum(ALPHA.*MUFVEC);

end

