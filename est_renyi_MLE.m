function est = est_renyi_MLE(samp,alpha)
%est_renyi_MLE  ML estimation of Renyi entropy (in bits) of the input data            
%
% This function returns a scalar MLE of the Renyi entropy of samp when samp 
% is a vector, or returns a (row-) vector containing the clolumn-wise MLEs
% of the Renyi entropy when samp is a matrix.
%
% Input:
% ----- samp: a vector or matrix which can only contain integers. The input
%             data type can be any interger classes such as uint8/int8/
%             uint16/int16/uint32/int32/uint64/int64, or floating-point 
%             such as single/double. 
% ----- alpha: the order of Renyi entropy (must be non-negative). 
% 
% Output:
% ----- est: Renyi entropy (in bits) of the input vector or that of each 
%            column of the input matrix. The output data type is double. 

if ~isequal(samp, fix(samp))
    error('Input sample must only contain integers.');
end

if isrow(samp)
    samp = samp.';
end
[n, wid] = size(samp);

f = find([diff(sort(samp))>0; true(1,wid)]);   
f = accumarray({filter([1;-1],1,f),ceil(f/n)},1);  

prob = (1:size(f,1))/n;
if alpha < 0
    error('The order of Renyi entropy must be non-negative.');
elseif alpha == 1
    prob_mat = -prob.*log2(prob);
    est = prob_mat * f;
else
    prob_mat = prob.^alpha;
    est = log2(prob_mat*f)/(1-alpha);
end