function [moms] = moment_stats(vec)
% This function computes mean and standard deviation of input vector.

moms=zeros(2,1);
moms(1)=mean(vec);
moms(2)=std(vec);

end

