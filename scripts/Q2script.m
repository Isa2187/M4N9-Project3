load('PortfolioData.mat'); 

%Solve the optimisation problem using LDL' factorisation
gammaValues = (0.1:0.1:0.9)';
n = length(gammaValues);
returnValuesQ21 = zeros(n,1);
riskValuesQ21 = zeros(n,1);
iterCountsQ21 = zeros(n,1);
timesQ21 = zeros(n,1);
largest2x1 = zeros(n,2);

for i = 1:n
    gamma = gammaValues(i);
    
    tic;
    [xn1,ret,risk,iterCount] = portOptimiseLDLT(gamma,sigma,pbar);
    timesQ21(i) = toc;
    
    returnValuesQ21(i) = ret;
    riskValuesQ21(i) = risk;
    iterCountsQ21(i) = iterCount;
    
    xnSort1 = sort(xn1);
    largest2x1(i,1) = xnSort1(end);
    largest2x1(i,2) = xnSort1(end-1);
end

%Solve the optimisation problem using the CG method
returnValuesQ23 = zeros(n,1);
riskValuesQ23 = zeros(n,1);
iterCountsQ23 = zeros(n,1);
timesQ23 = zeros(n,1);
largest2x2 = zeros(n,2);

for i = 1:n
    gamma = gammaValues(i);
    
    tic;
    [xn2,ret,risk,iterCount] = portOptimiseCG(gamma,sigma,pbar);
    timesQ23(i) = toc;
    
    returnValuesQ23(i) = ret;
    riskValuesQ23(i) = risk;
    iterCountsQ23(i) = iterCount;
    
    xnSort2 = sort(xn2);
    largest2x2(i,1) = xnSort2(end);
    largest2x2(i,2) = xnSort2(end-1);
end

%Solve the optimisation problem using the preconditioned CG method
returnValuesQ24 = zeros(n,1);
riskValuesQ24 = zeros(n,1);
iterCountsQ24 = zeros(n,1);
timesQ24 = zeros(n,1);
largest2x3 = zeros(n,2);

for i = 1:n
    gamma = gammaValues(i);
    
    tic;
    [xn2,ret,risk,iterCount] = portOptimisePCCG(gamma,sigma,pbar);
    timesQ24(i) = toc;
    
    returnValuesQ24(i) = ret;
    riskValuesQ24(i) = risk;
    iterCountsQ24(i) = iterCount;
    
    xnSort3 = sort(xn2);
    largest2x3(i,1) = xnSort3(end);
    largest2x3(i,2) = xnSort3(end-1);
end

%Plot the run times of the LDL' and CG methods
figure()
plot(gammaValues,timesQ21)
hold on
plot(gammaValues,timesQ23)
title('Graph Of Run Time Against \gamma')
xlabel('Gamma, \gamma')
ylabel('Run Time, t')
legend('LDL^T Method','CG Method')
hold off

%Plot the run times of the preconditioned CG and normal CG methods
figure()
plot(gammaValues,timesQ24)
hold on
plot(gammaValues,timesQ23)
title('Graph Of Run Time Against \gamma')
xlabel('Gamma, \gamma')
ylabel('Run Time, t')
legend('Preconditioned CG Method','CG Method')
hold off