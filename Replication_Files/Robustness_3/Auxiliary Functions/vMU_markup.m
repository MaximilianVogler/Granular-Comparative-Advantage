function [EPSX] = vMU_markup(sigma,SHX,BER)
% This function computes the variable markup depending on the kind of
% competition. BER==1 corresponds to Bertrand and BER==0 to Cournot
% competition.

    if BER==1
        EPSX=sigma*(1-SHX)+SHX;
    elseif BER==0
        EPSX=1./((1/sigma)*(1-SHX)+SHX); 
    else
        error('BER needs to equal either 0 or 1')
    end
end