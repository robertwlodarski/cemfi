function result = sum_of_types(vector)
    crises_predicted = vector>0.5
    type_1 = (crises_predicted(34)~=1)/length(vector)
    big_sum = sum(crises_predicted)
    correct_prediction =(crises_predicted(34)>0)
    type_2 = (big_sum-correct_prediction)/length(vector)
    result(i,1) = type_1+type_2
end