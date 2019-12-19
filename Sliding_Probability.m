
function [Prob_Vector] = Sliding_Probability(vector, maxval, increments)
%This computes a sliding probability of whether the vector values are
%greater than a criterion value which moves from zero to maxval in
%increments of the input variable, increments
%   Detailed explanation goes here

critvalues = 0:maxval/increments:maxval;

if maxval == 0
    
    Prob_Vector = zeros(1, length(critvalues));
    return
end

for i = 1:length(critvalues)
    
    Prob_Vector(i) = sum(vector > critvalues(i))/length(vector);
    
end
end