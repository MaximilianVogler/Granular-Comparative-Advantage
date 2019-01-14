function [K,KF,PHI,PHIF,LAMBDA,LAMBDAF,MU,MUF,LAMBDAFVEC,mom] = PEreplication_thetas(sigma,F,tau,RTS,RTL,ZHS,ZFS,ZHL,ZFL,w,wF,Y0,YF0,small,vMU,BER,paretonb,VEC_CD_THETA)
% Takes parameters, random draws, and initial guess as inputs. Outputs are
% K, Lambda, and aggregate markup Mu in both countries. 

S=size(ZHL,2)+size(ZHS,2); % fix number of sectors

ALPHA = VEC_CD_THETA(:,1);
THETA = VEC_CD_THETA(:,2);

ALPHAS=ALPHA(small);
ALPHAL=ALPHA(~small);

THETAS=THETA(small);
THETAL=THETA(~small);


% Initialize vectors that contain PE data
PHIHVEC=zeros(1,S);
MUHVEC=zeros(1,S);
PHIFVEC=zeros(1,S);
MUFVEC=zeros(1,S);
KVEC=zeros(1,S);
KFVEC=zeros(1,S);
LAMBDAHVEC=zeros(1,S);
LAMBDAFVEC=zeros(1,S); 
TOP1=zeros(1,S);
PARETO=zeros(1,S);

% Solve inner loops separately for small and large sectors (for speed)
[KVEC(small),KFVEC(small),PHIHVEC(small),PHIFVEC(small),MUHVEC(small),MUFVEC(small),LAMBDAHVEC(small),LAMBDAFVEC(small),TOP1(small),PARETO(small)]=...
    Inner_Loops_thetas(sigma,F,tau,ALPHAS',THETAS',RTS,ZHS,ZFS,w,wF,Y0,YF0,vMU,BER,paretonb);
[KVEC(~small),KFVEC(~small),PHIHVEC(~small),PHIFVEC(~small),MUHVEC(~small),MUFVEC(~small),LAMBDAHVEC(~small),LAMBDAFVEC(~small),TOP1(~small),PARETO(~small)]=...
    Inner_Loops_thetas(sigma,F,tau,ALPHAL',THETAL',RTL,ZHL,ZFL,w,wF,Y0,YF0,vMU,BER,paretonb);

%Take averages of vectors to obtain aggregates       
% K=mean(KVEC); % number of firms in home
% KF=mean(KFVEC);
% PHI=mean(PHIHVEC); % expected import share in home
% PHIF=mean(PHIFVEC);
% LAMBDA=mean(LAMBDAHVEC); % realized import share in home
% LAMBDAF=mean(LAMBDAFVEC);
% MU=1./sum(ALPHA'.*MUHVEC); % aggregate markup in home
% MUF=1./sum(ALPHA'.*MUFVEC);

GAMMAFVEC=LAMBDAFVEC-PHIFVEC;
K=sum(KVEC);            % number of firms in home
KF=sum(KFVEC);
PHI=sum(ALPHA'.*PHIHVEC);      % expected import share in home
PHIF=sum(ALPHA'.*PHIFVEC);
LAMBDA=sum(ALPHA'.*LAMBDAHVEC);
LAMBDAF=sum(ALPHA'.*LAMBDAFVEC);
MU=(1-LAMBDA)/sum(ALPHA'.*MUHVEC);
MUF=(1-LAMBDAF)/sum(ALPHA'.*MUFVEC);
DEC1=var(GAMMAFVEC(GAMMAFVEC>0))/var(LAMBDAFVEC(LAMBDAFVEC>0));
X = LAMBDAFVEC.*VEC_CD_THETA(:,1)'.*YF0;
DEC2 = var(log(LAMBDAFVEC(LAMBDAFVEC>0)))/var(log(X(X>0)));
TOP1=mean(TOP1(TOP1>0));
PARETO=mean(PARETO(PARETO>0));
KAPPA=mean(VEC_CD_THETA(:,2)/4);
mom=[DEC1,DEC2,KAPPA,PARETO,TOP1];

end

