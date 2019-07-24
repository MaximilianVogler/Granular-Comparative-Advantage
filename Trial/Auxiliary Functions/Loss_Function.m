%% This is the loss function to be minimized
function [lfun] = Loss_Function(muT,sigmaT,tau,kappa,f,sigma,UHS,UHL,UFS,UFL,rtdraws,small,MH_small,MH_large,ALPHA,w,wF,vMU,BER,datamoments,W)

theta=(sigma-1)*kappa;
f = f * (4.93*.43*10^(-5));
F=f/sigma;

% Given mu_T and sigma_T draw sectoral productivity T_z for each sector z 
  RT=exp(muT+sigmaT*rtdraws);
  RTS=RT(small);
  RTL=RT(~small);
    
% Draw productivities phi 
  ZHS=(UHS./(repmat(RTS,MH_small,1))).^(-1/theta);
  ZFS=UFS.^(-1/theta);
  ZHL=(UHL./(repmat(RTL,MH_large,1))).^(-1/theta);                        
  ZFL=UFL.^(-1/theta);
    
% Guess home and foreign output Y and fix labor
  Y0=126;
  YF0=225;
  L0=100;
  
[~,Y,YF,~,KHH,TOP1,TOP3,XS,YXS,LAMBDAHVEC,LAMBDAFVEC] = GEreplication_vectorized(sigma,theta,F,tau,ALPHA,RTS,RTL,ZHS,ZFS,ZHL,ZFL,w,wF,L0,Y0,YF0,small,vMU,BER,0);
    
moments = Moments(KHH,TOP1,TOP3,LAMBDAHVEC,LAMBDAFVEC,XS,YXS,ALPHA,Y,YF);  
               
lfun = (moments'-datamoments)'*W*(moments'-datamoments);

end