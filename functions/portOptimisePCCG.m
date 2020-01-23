%Define a function that takes as its inputs a value of gamma, the spd
%covariance matrix sigma, and the mean prices pbar.
%The function returns the expected return, risk and optimal solution, where
%this is found using the preconditioned CG method.
%It also outputs iterCount, which will be the total number of times
%preconditioned CG was used to solve a linear system
function [xn,ret,risk,iterCount] = portOptimisePCCG(gamma,sigma,pbar)

    clear val;
    m = size(pbar,1);

    %Set the tolerance for Newton convergence
    tol = 1e-8;

    %Set the step size alhpa, and the initial guess of the minimizer
    alpha = 0.1;
    xn = ones(m,1)./m;

    i = 1;    
    
    %Initialise iterCount, which will be the number of times the 
    %preconditioned CG method is used to solve a linear system
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

            % Set up for the matrix, Mbar
            Phi = diag(t./(xn.^2));
            M = 2.0*gamma*sigma + Phi;
            Mbar = Z' * M * Z;

            %Solve the linear system using the preconditioned CG method,
            %where the preconditioner matrix, C, is a diagonal matrix whose
            %entries are the diagonal of Mbar
            k = 0;
            rkm2 = gbar;
            diagMbar = diag(Mbar);
            preconMmat = diag(diagMbar);
            zkm2 = rkm2./diagMbar;
            xbar = zeros(m,1);
            res1 = 1;

            while (res1 > tol)
                k = k + 1;
                if k == 1
                    p = zkm2;
                    zkm1 = zkm2;
                    rkm1 = rkm2;
                else
                    beta = ( rkm1' * zkm1 )/( rkm2' * zkm2 );
                    p = zkm1 + beta * p;
                end
                
                a = ( rkm1' * zkm1 )/ (p' * Mbar * p);
                xbar = xbar + a * p;
                rkm2 = rkm1;
                rkm1 = rkm1 - a * Mbar * p;
                zkm2 = zkm1;            
                zkm1 = rkm1./diagMbar;
                res1 = norm(rkm1);
                
                %Update the iteration count
                iterCount = iterCount + 1;
            end
            
            %Obtain dx from the solution to the other above linear system
            dx = Z*xbar;

            %Update the minimizer
            xnp1 = xn + alpha*dx(1:m);

            %Compute change in the minimizer to check for convergence. 
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