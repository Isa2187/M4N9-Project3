%Define a function that takes as its inputs a value of gamma, the spd
%covariance matrix sigma, and the mean prices pbar.
%The function returns the expected return, risk and optimal solution, where
%this is found using the CG method.
%It also outputs iterCount, which will be the total number of times CG was
%used to solve a linear system
function [xn,ret,risk,iterCount] = portOptimiseCG(gamma,sigma,pbar)

    clear val;
    %Set m to be the number of possible investments
    m = size(pbar,1);

    %Set the tolerance for Newton convergence
    tol = 1e-8;

    %Set the step size alhpa, and the initial guess of the minimizer
    alpha = 0.1;
    xn = ones(m,1)./m;

    i = 1;
    
    %Initialise iterCount, which will be the number of times the CG 
    %method is used to solve a linear system
    iterCount = 0;
    
    %Set up the matrix Z as described in my report
    Z = zeros(m,m);
    for j = 1:m-1
        Z(j,j) = 1;
        Z(j+1,j) = -1;
    end    
    Z(1,m) = -1;
    Z(m,m) = 1;

    for n = 4:11
    %Fix the parameter t 
    t = 10^-n;
    res = 1;
    
        while (res > tol)
            %Set up the right hand side, gbar
            g = t./(xn) + (1 - gamma)*pbar - 2.0*gamma*sigma*xn;
            gbar = Z' * g;
            
            %Set up for the matrix, Mbar
            Phi = diag(t./(xn.^2));
            M = 2.0*gamma*sigma + Phi;            
            Mbar = Z' * M * Z;

            %Solve the linear system using the CG method
            xbar = zeros(m,1);
            rkm2 = gbar;
            res1 = 1;
            k = 1;

            while (res1 > tol)
                if k == 1
                    p = rkm2;
                    rkm1 = rkm2;
                else
                    beta = ( norm(rkm1)^2 )/( norm(rkm2)^2 );
                    p = rkm1 + beta * p;
                end

                a = ( norm(rkm1)^2 )/ (p' * Mbar * p);
                xbar = xbar + a * p;
                rkm2 = rkm1;
                rkm1 = rkm1 - a * Mbar * p;
                k = k + 1;
                res1 = norm(Mbar*xbar - gbar)./norm(xbar);    
                
                %Update the iteration count
                iterCount = iterCount + 1;
            end
            
            %Obtain dx from the solution to the other above linear system
            dx = Z*xbar;
            
            % Update the minimizer
            xnp1 = xn + alpha*dx(1:m);

            %Compute change in the minimizer to check for convergence
            res = norm(xnp1 - xn)./norm(xnp1);

            %Compute the value of the total cost function
            val(i) = gamma*xn'*sigma*xn - (1 - gamma)*pbar'*xn - t.*sum(log(xn)); 
            xn = xnp1;
            i = i + 1;
        end
    end

    %Compute the portfolio return
    ret = pbar'*xn;

    %Compute the risk
    risk = xn'*sigma*xn;
end