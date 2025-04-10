function value = type2error(x,vector)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
prediction = (vector>x)
value = sum(prediction([1:33,35:79],1)==1)
end