clear, close all, clc

rng(0)

alpha = 1.5;    % The order of Renyi entropy
mc_times = 50;  % The number of Monte-Carlo trials for each alphabet size  
num = 15;
C = 1;  
record_S = ceil(logspace(2, 7, num));
record_n = ceil(C*record_S./log(record_S));
 
% The Chebfun package provides the functionality for best polynomial approximation
addpath(genpath(fullfile(pwd,'Chebfun v5.3.0')));  
twonum = rand(2,1); 

for iter = num:-1:1  
    S = record_S(iter)
    n = record_n(iter)
    dist = betarnd(twonum(1),twonum(2),S,1);
    dist = dist/sum(dist);                        
    true_S(iter) = renyi_true(dist,alpha);      
    samp = randsmpl(dist, n, mc_times, 'int32');           
    record_JVHW = est_renyi_JVHW(samp,alpha);    
    record_MLE = est_renyi_MLE(samp,alpha);
    JVHW_err(iter) = mean(abs(record_JVHW - true_S(iter)));  
    MLE_err(iter) = mean(abs(record_MLE - true_S(iter)));      
end
 
figure(1) 
p(2) = plot(record_S./record_n, JVHW_err,'r-','LineWidth',2,'MarkerSize',8);  hold on;
p(1) = plot(record_S./record_n, MLE_err,'k-s','LineWidth',2,'MarkerSize',10);
legend(p,{'MLE','JVHW'},'Location','northwest','Interpreter','latex','FontSize',14);
xlabel('S/n')
ylabel('Mean Absolute Error')
title(sprintf('Renyi entropy $H_\\alpha(P)$ estimation for $\\alpha=%g$ and $n=%g\\frac{S}{\\log S}$',alpha,C),'Interpreter','latex')
xlim([4, 16.5])
grid on 
