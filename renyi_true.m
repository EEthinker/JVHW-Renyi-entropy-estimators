function Ha = renyi_true(p,alpha)
% renyi_true  computes base-2 Renyi entropy of the input discrete distribution.
%
% This function returns the scalar Renyi entropy Ha when the input p is a 
% probability vector, or returns a row vector containing the Renyi entropy 
% of each column of the input probability matrix p. 
%
% Remarks: 
%   --- Renyi entropy is a non-decreasing function in alpha: 
%          log2(S) = H_0 >= H_1 >= H_2 >= ... >= H_infty
%   --- Speicial cases of Renyi entropy include: 
%              H_0 = log2(S) 
%              H_1 = H(p) --- Shannon entropy 
%              H_2 = -log2(sum(p.^2))
%              H_infty = -log2(max(p))

% Error-check of the input distribution
p0 = p(:);
if any(imag(p0)) || any(isinf(p0)) || any(isnan(p0)) || any(p0<0) || any(p0>1)
    error('The probability elements must be real numbers between 0 and 1.');
elseif any(abs(sum(p)-1) > sqrt(eps))
    error('Sum of the probability elements must equal 1.');
end

if alpha < 0
    error('The order of Renyi entropy must be non-negative.');
elseif alpha == 1  
    Ha = -sum(log2(p.^p));  
else
    Ha = log2(sum(p.^alpha))/(1-alpha);
end