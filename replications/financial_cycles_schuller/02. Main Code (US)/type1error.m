function value = type1error(x,vector)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
prediction = (vector>x)
indicator = (prediction(34,1)~=1)
value = indicator
end