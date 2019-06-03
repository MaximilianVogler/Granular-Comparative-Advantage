function [iter,Y,YF,LF,KHH,TOP1,TOP3,XS,YXS,LAMBDAHVEC,LAMBDAFVEC,varphi_bar] = GEreplication_vectorized(sigma,theta,F,tau,ALPHA,RTS,RTL,ZHS,ZFS,ZHL,ZFL,w,wF,L,Y0,YF0,small,vMU,BER,paretonb)
% Matrix code to compute GE with w/w* fixed. Takes parameters, random draws, 
% and initial guess as inputs. Outputs GE values of Y, Y*, L*.

diff=1;
tol=1e-3;
iter=0;

tic
while diff>tol && iter<51
    % Solve PE
    [K,KF,~,~,LAMBDA,LAMBDAF,MU,MUF,KHH,TOP1,TOP3,XS,YXS,LAMBDAHVEC,LAMBDAFVEC,~,~,varphi_bar]=PEreplication_vectorized(sigma,theta,F,tau,ALPHA,RTS,RTL,ZHS,ZFS,ZHL,ZFL,w,wF,Y0,YF0,small,vMU,BER,paretonb,0);
    
    % Set up linear system composed of (A10) and (A14)
    MAT=[(1-LAMBDA)*(MU-1)/MU-1 LAMBDAF*(MUF-1)/MUF;LAMBDA -LAMBDAF];
    VEC=[-w*L+LAMBDAF*wF*F*KF+(1-LAMBDA)*w*F*K;LAMBDA*w*F*K-LAMBDAF*wF*F*KF];
    SOL=MAT\VEC;

    Y=SOL(1);
    YF=SOL(2);

    diff=abs(Y-Y0)+abs(YF-YF0);
    iter=iter+1;
    if iter<=5
        step=1;
    else
        step=0.5;
    end

    Y0=Y0+step*(Y-Y0);
    YF0=YF0+step*(YF-YF0);
end

% Use (A13) to solve for L* 
LF=1/wF*(wF*F*KF+YF*(1-LAMBDAF)/MUF+Y*LAMBDA/MU);

time=toc;

end