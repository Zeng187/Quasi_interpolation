function [value] = rbf(r, approxArg, s)
%RBF rbf距离场函数定义
sqr = r^2;
value = (1-r)^4*(4*r+1)*((approxArg+sqr)^s );
end

