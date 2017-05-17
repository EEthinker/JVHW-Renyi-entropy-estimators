function est = est_renyi_JVHW(samp,alpha)
% This function returns a scalar estimate of the Renyi entropy of samp when
% samp is a vector, or returns a row-vector containing the clolumn-wise
% estimates of the Renyi entropy when samp is a matrix.
%
% Input:
% ----- samp: a vector or matrix which can only contain integers. The input
%             data type can be any interger classes such as uint8/int8/
%             uint16/int16/uint32/int32/uint64/int64, or floating-point
%             such as single/double.
% ----- alpha: the order of Renyi entropy (must be non-negative).
%
% Output:
% ----- est:  base-2 Renyi entropy of the input vector or that of each
%             column of the input matrix. The output data type is double.

if alpha < 0
    error('The order of Renyi entropy must be non-negative.');
end
if ~isequal(samp, fix(samp))
    error('Input sample must only contain integers.');
end
if isrow(samp)
    samp = samp.';
end
[n, wid] = size(samp);

f = find([diff(sort(samp))>0; true(1,wid)]);
f = accumarray({filter([1;-1],1,f),ceil(f/n)},1);

if alpha == 1
    est = est_entro_JVHW(n,f); 
else
    chooseIntmethod = fix(alpha) == alpha || ...
                      2.00 < alpha && alpha <= 2.15 && n < 150   || ...
                      2.15 < alpha && alpha <= 2.25 && n < 65    || ...
                      2.25 < alpha && alpha <= 2.35 && n < 65    || ...
                      2.35 < alpha && alpha <= 2.45 && n < 50    || ...
                      2.45 < alpha && alpha <= 2.55 && n < 120   || ...
                      2.55 < alpha && alpha <= 2.65 && n < 200   || ...
                      2.65 < alpha && alpha <= 2.75 && n < 200   || ...
                      2.75 < alpha && alpha <= 2.85 && n < 450   || ...
                      2.85 < alpha && alpha <= 2.95 && n < 1000  || ...
                      2.95 < alpha && alpha <= 3.05              || ...
                      3.05 < alpha && alpha <= 3.15 && n < 3000  || ...
                      3.15 < alpha && alpha <= 3.25 && n < 3000  || ...
                      3.25 < alpha && alpha <= 3.35 && n < 1800  || ...
                      3.35 < alpha && alpha <= 3.45 && n < 1800  || ...
                      3.45 < alpha && alpha <= 3.55 && n < 3500  || ...
                      3.55 < alpha && alpha <= 3.65 && n < 6000  || ...
                      3.65 < alpha && alpha <= 3.75 && n < 6500  || ...
                      3.75 < alpha && alpha <= 3.85 && n < 1e4   || ...
                      3.85 < alpha;
    if chooseIntmethod
        powerSumEst = integerMethod(alpha, n, f);
    else
        powerSumEst = polyAppMethod(alpha, n, f);
    end
    est = log2(powerSumEst)/(1-alpha);
end 


function powerSumEst = integerMethod(alpha, n, f)
% The integerMethod function estimates the powerSum = sum(pi^alpha,i=1..S)
% by pluging in various modified forms of the original unbiased estimator
% of pi^alpha when alpha is integer. The integer-order unbiased estimator
% is given by: pi^alpha = prod((i-r)/(n-r),i=0..alpha-1)

n2f1LogRegime = {[]; [];
    [0  0.361  0.362  0.500  0.790  1.000  Inf];     
    [0  0.220  0.318  0.550  0.750  1.000  Inf];     
    [0  0.220  0.294  0.479  0.790  1.120  Inf];    
    [0  0.200  0.283  0.460  0.747  1.200  Inf];    
    [0  0.200  0.280  0.436  0.700  1.250  Inf];     
    [0  0.190  0.278  0.428  0.684  1.300  Inf];    
    [0  0.180  0.275  0.397  0.655  1.310  Inf];    
    [0  0.175  0.275  0.375  0.630  1.290  Inf];    
    [0  0.160  0.275  0.365  0.620  1.320  Inf];    
    [0  0.153  0.270  0.362  0.620  1.350  Inf];  
    [0  0.148  0.270  0.350  0.610  1.350  Inf];    
    [0  0.147  0.270  0.350  0.610  1.350  Inf];    
    [0  0.145  0.270  0.350  0.610  1.360  Inf];    
    [0  0.140  0.270  0.342  0.583  1.365  Inf];    
    [0  0.138  0.270  0.340  0.560  1.370  Inf];    
    [0  0.137  0.270  0.330  0.560  1.375  Inf];    
    [0  0.131  0.270  0.330  0.560  1.378  Inf];    
    [0  0.129  0.268  0.328  0.560  1.378  Inf];    
    };

a = round(alpha);
orderCorrection = alpha/a;
r = (0:a-1).';
X = 1:size(f,1);
powerEst = prod(bsxfun(@rdivide, bsxfun(@minus, X, r), n-r),1).^orderCorrection;
powerSumEst = powerEst*f;
isZero = powerSumEst == 0;  
f1 = f(:,isZero);

method1 = powerSumEst;    
method2 = powerSumEst;   
method3 = powerSumEst;  
method4 = powerSumEst;  
method5 = powerSumEst;   
method6 = powerSumEst;  
 
powerEst1 = bsxfun(@rdivide, bsxfun(@minus, X, r), n-r);
idn = powerEst1 <= 0;
powerEst1(idn) = repelem(powerEst1(1,:),sum(idn));
powerEst1 = prod(powerEst1,1).^orderCorrection;
powerEst2 = prod(bsxfun(@max, bsxfun(@rdivide,bsxfun(@minus, X, r),n-r), 1./(n-X+1)),1).^orderCorrection; 
powerEst3 = prod(max(bsxfun(@minus,X, r),1-1/a+bsxfun(@minus,X, r)/a)./bsxfun(@max,n-r,bsxfun(@minus,n-r/a,(1-1/a)*(X-1))),1).^orderCorrection;
powerEst4 = prod(bsxfun(@max,bsxfun(@minus,X, r),1-r/a)./bsxfun(@max,n-r,bsxfun(@minus,n-X+1,r/a)),1).^orderCorrection;  
powerEst5 = prod(bsxfun(@max, bsxfun(@minus, X, r), 1./(r+1))./bsxfun(@max, n-r, bsxfun(@plus, n-X, 1./(r+1)))).^orderCorrection;

method1(isZero) = powerEst1*f1;
method2(isZero) = powerEst2*f1;
method3(isZero) = powerEst3*f1;
method4(isZero) = powerEst4*f1;

temp = bsxfun(@times, powerEst4.', f1);
temp(temp==0) = Inf;
method5(isZero) = min(temp);

temp = bsxfun(@times, powerEst5.', f1);
temp(temp==0) = Inf;
method6(isZero) = min(temp);

allMethods = [method6; method5; method4; method3; method2; method1];
pos = discretize(log(n./f1(1,:)), n2f1LogRegime{min(max(a,3),20)}) + (0:sum(isZero)-1)*size(allMethods,1);
powerSumEst(isZero) = allMethods(pos);

isNaN = isnan(powerSumEst);  
if any(isNaN)
    powerSumEst(isNaN) = powerEst4*f(:,isNaN);   
end


function powerSumEst = polyAppMethod(alpha, n, f)
persistent poly_renyi current_poly_renyi alpha_samp alpha_old
if isempty(poly_renyi)
    load poly_coeff_renyi.mat poly_renyi alpha_samp;
end
if ~isequal(alpha,alpha_old)  
    alpha_old = alpha;
    alpha_find = alpha == alpha_samp;
    if any(alpha_find)
        current_poly_renyi = poly_renyi{alpha_find};
    else
        current_poly_renyi = renyi_poly(alpha);
    end
end 

K = min([4+ceil(1.2*log(n)), 22, n]); 
% Piecewise linear/quadratic fit
V1 = [0.3303 0.4679];   
V2 = [-0.530556484842359,1.09787328176926,0.184831781602259];    
f1nonzero = f(1,:) > 0;
c_1 = zeros(1, size(f,2));
if K >= 3 && any(f1nonzero)
    if alpha < 6
        if n < 200
            c_1(f1nonzero) = polyval(V1, log(n./f(1,f1nonzero)));   
        else
            n2f1_small = f1nonzero & log(n./f(1,:)) <= 2;
            n2f1_large = f1nonzero & log(n./f(1,:)) > 2;
            c_1(n2f1_small) = polyval(V2, log(n./f(1,n2f1_small))); 
            c_1(n2f1_large) = polyval(V1, log(n./f(1,n2f1_large))); 
        end
        c_1(f1nonzero) = max(c_1(f1nonzero), 1/(1.9*log(n)));  
    else
        c_1(f1nonzero) = polyval(V2, log(n./f(1,f1nonzero)));
    end
end
g_coeff = current_poly_renyi{K};   
cnt = (1:size(f,1)).';
thres = 4*c_1*log(n)/n;
[T, P] = meshgrid(thres,cnt/n);
ratio = min(max(2*P./T-1,0),1);  
MLE = P.^alpha;  
r = 0:K-1;
X_cum = cumprod([ones(length(cnt),1) bsxfun(@rdivide, bsxfun(@minus, cnt, r), n-r)],2);
G = diag(g_coeff)*bsxfun(@power, thres, (alpha:-1:alpha-K).');
polyApp = X_cum*G;

polyfail = isnan(polyApp) | isinf(polyApp);
polyApp(polyfail) = MLE(polyfail);
prob_mat = ratio.*MLE + (1-ratio).*polyApp;
powerSumEst = sum(f.*prob_mat, 1);


function est = est_entro_JVHW(n,f)
% This function returns a scalar JVHW estimate of Shannon entropy (in bits)
% when f is a column vector, or returns a row vector containing the JVHW
% estimate of each column of f when f is a matrix.

K = min([4+ceil(1.2*log(n)), 22, n]);  
persistent poly_entro;
if isempty(poly_entro)
    load poly_coeff_entro.mat poly_entro;
end
g_coeff = poly_entro{K};  
V1 = [0.3303 0.4679];
V2 = [-0.530556484842359,1.09787328176926,0.184831781602259];
f1nonzero = f(1,:) > 0;
c_1 = zeros(1, size(f,2));
if K >= 3 && any(f1nonzero)
    if n < 200
        c_1(f1nonzero) = polyval(V1, log(n./f(1,f1nonzero)));
    else
        n2f1_small = f1nonzero & log(n./f(1,:)) <= 1.5;
        n2f1_large = f1nonzero & log(n./f(1,:)) > 1.5;
        c_1(n2f1_small) = polyval(V2, log(n./f(1,n2f1_small)));
        c_1(n2f1_large) = polyval(V1, log(n./f(1,n2f1_large)));
    end
    c_1(f1nonzero) = max(c_1(f1nonzero), 1/(1.9*log(n)));  
end
cnt = (1:size(f,1)).';
thres = 4*c_1*log(n)/n;
[T, P] = meshgrid(thres,cnt/n);
ratio = min(max(2*P./T-1,0),1);  
MLE = -P.*log(P) + 1/(2*n);
r = 0:K-1;
X_cum = cumprod([ones(length(cnt),1) bsxfun(@rdivide, bsxfun(@minus, cnt, r), n-r)],2);
G = diag(g_coeff)*bsxfun(@power, thres, (1:-1:1-K).');
polyApp = X_cum*G - P.*log(T);
polyfail = ~isfinite(polyApp); 
polyApp(polyfail) = MLE(polyfail);
prob_mat = ratio.*MLE + (1-ratio).*polyApp;
prob_mat = max(prob_mat,0);
est = sum(f.*prob_mat, 1)/log(2);