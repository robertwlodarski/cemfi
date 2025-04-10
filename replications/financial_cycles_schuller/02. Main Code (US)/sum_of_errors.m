function output_arg = sum_of_errors(x,vector)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here
prediction = (vector>x)
value1 = (prediction(34,1)~=1)
value2 = sum(prediction([1:33,35:79],1)==1)
output_arg = (value1 + value2)
end