function [K,KF,PHI,PHIF,LAMBDA,LAMBDAF,MU,MUF,KHH,TOP1,TOP3,XS,YXS,LAMBDAHVEC,LAMBDAFVEC,PHIFVEC,mom] = PEreplication_vectorized(sigma,theta,F,tau,ALPHA,RTS,RTL,ZHS,ZFS,ZHL,ZFL,w,wF,Y0,YF0,small,vMU,BER,paretonb,AddMom)
% Takes parameters, random draws, and initial guess as inputs. Outputs are
% K, Lambda, and aggregate markup Mu in both countries. 

S=size(ZHL,2)+size(ZHS,2); % fix number of sectors

ALPHAS=ALPHA(small);
ALPHAL=ALPHA(~small);

% Initialize vectors that contain PE data
PHIHVEC=zeros(1,S);
MUHVEC=zeros(1,S);
PHIFVEC=zeros(1,S);
MUFVEC=zeros(1,S);
KVEC=zeros(1,S);
KFVEC=zeros(1,S);
LAMBDAHVEC=zeros(1,S);
LAMBDAFVEC=zeros(1,S); 
KHH=zeros(1,S);
TOP1=zeros(1,S);
TOP3=zeros(1,S);
XS=zeros(1,S);
YXS=zeros(1,S);
PARETO=zeros(1,S);
mom=zeros(1,S);
if AddMom == 0
    % Solve inner loops separately for small and large sectors (for speed)
    [KVEC(small),KFVEC(small),PHIHVEC(small),PHIFVEC(small),MUHVEC(small),MUFVEC(small),LAMBDAHVEC(small),LAMBDAFVEC(small),KHH(small),TOP1(small),TOP3(small),XS(small),YXS(small),~]=...
        Inner_Loops(sigma,theta,F,tau,ALPHAS',RTS,ZHS,ZFS,w,wF,Y0,YF0,vMU,BER,paretonb,AddMom,S);
    [KVEC(~small),KFVEC(~small),PHIHVEC(~small),PHIFVEC(~small),MUHVEC(~small),MUFVEC(~small),LAMBDAHVEC(~small),LAMBDAFVEC(~small),KHH(~small),TOP1(~small),TOP3(~small),XS(~small),YXS(~small),~]=...
        Inner_Loops(sigma,theta,F,tau,ALPHAL',RTL,ZHL,ZFL,w,wF,Y0,YF0,vMU,BER,paretonb,AddMom,S);
elseif AddMom == 1
    [KVEC(small),KFVEC(small),PHIHVEC(small),PHIFVEC(small),MUHVEC(small),MUFVEC(small),LAMBDAHVEC(small),LAMBDAFVEC(small),KHH(small),TOP1(small),TOP3(small),XS(small),YXS(small),PARETO(small)]=...
        Inner_Loops(sigma,theta,F,tau,ALPHAS',RTS,ZHS,ZFS,w,wF,Y0,YF0,vMU,BER,paretonb,AddMom,S);
    [KVEC(~small),KFVEC(~small),PHIHVEC(~small),PHIFVEC(~small),MUHVEC(~small),MUFVEC(~small),LAMBDAHVEC(~small),LAMBDAFVEC(~small),KHH(~small),TOP1(~small),TOP3(~small),XS(~small),YXS(~small),PARETO(~small)]=...
        Inner_Loops(sigma,theta,F,tau,ALPHAL',RTL,ZHL,ZFL,w,wF,Y0,YF0,vMU,BER,paretonb,AddMom,S);
    PARETO=mean(PARETO(PARETO>0));
    DEC1=(1-var(PHIFVEC)/var(LAMBDAFVEC))*100;
    X = LAMBDAFVEC.*ALPHA';
    DEC2 = var(log(LAMBDAFVEC(LAMBDAFVEC>0)))/var(log(X(LAMBDAFVEC>0)))*100;
    KAPPA = theta/(sigma-1);
    TOP1=mean(TOP1(TOP1>0));
    mom=[DEC1,DEC2,KAPPA,PARETO,TOP1];    
else 
    error('Pareto needs to be either 0 or 1')
end
% %Take averages of vectors to obtain aggregates       
% K=mean(KVEC); % number of firms in home
% KF=mean(KFVEC);
% PHI=mean(PHIHVEC); % expected import share in home
% PHIF=mean(PHIFVEC);
% LAMBDA=mean(LAMBDAHVEC); % realized import share in home
% LAMBDAF=mean(LAMBDAFVEC);
% MU=1./sum(ALPHA'.*MUHVEC); % aggregate markup in home
% MUF=1./sum(ALPHA'.*MUFVEC);

K=sum(KVEC);            % number of firms in home
KF=sum(KFVEC);
PHI=sum(ALPHA'.*PHIHVEC);      % expected import share in home
PHIF=sum(ALPHA'.*PHIFVEC);
LAMBDA=sum(ALPHA'.*LAMBDAHVEC);
LAMBDAF=sum(ALPHA'.*LAMBDAFVEC);
MU=(1-LAMBDA)/sum(ALPHA'.*MUHVEC);
MUF=(1-LAMBDAF)/sum(ALPHA'.*MUFVEC);



end

