%Define a function that takes as its inputs a value of gamma, the spd
%covariance matrix sigma, and the mean prices pbar.
%The function returns the expected return, risk and optimal solution, where
%this is found using the LDL' factorisation method.
%It also outputs iterCount, which will be the total number of LDL'
%factorisations computed as we find the optimal solution
function [xn,ret,risk,iterCount] = portOptimiseLDLT(gamma,sigma,pbar)

    clear val;
    %Set m to be the number of possible investments
    m = size(pbar,1);

    %Set the tolerance for Newton convergence
    tol = 1e-8;

    % Set the step size alpha, the initial guess of the minimizer, and the
    %vector 1 appearing in A
    alpha = 0.1;
    xn = ones(m,1)./m;
    OneVec = ones(m,1);

    %Initialise iteration counter for the Newton method, iterCount, which
    %will also be the number of LDL' factorisations which are computed
    iterCount = 0;

    for n = 4:11
    %Fix the parameter t 
    t = 10^-n;
    res = 1;
    
        while (res > tol)
            %Set up the right hand side, r
            g = t./(xn) + (1 - gamma)*pbar - 2.0*gamma*sigma*xn;
            r = [g; 0];

            %Set up for the matrix, A
            A = ones(m+1);
            Phi = diag(t./(xn.^2));
            M = 2.0*gamma*sigma + Phi;
            A(1:m, 1:m) = M;
            A(m+1, m+1) = 0;

            %Solve the linear system A dx = r, using the algorithm provided
            %in Chapter 4 of Matrix Computations by Golub and Van Loan
            L = A;
            v = zeros(m+1,1);
            for j = 1:m+1
                for k = 1:m
                    v(k) = L(j,k)*L(k,k);
                end
                v(j) = L(j,j) - L(j,1:j-1)*v(1:j-1);
                L(j,j) = v(j);
                L(j+1:m+1,j) = ( L(j+1:m+1,j) - L(j+1:m+1,1:j-1)*v(1:j-1) )/v(j);
            end
            
            %Obtain the diagonal matrix, D
            D = diag(L);

            for j = 1:m+1
                L(j,j) = 1;
            end

            %Form the lower triangular matrix, L
            for j = 1:m+1
                for k = j+1:m+1
                    L(j,k) = 0;
                end
            end

            %Define LT to be L' in the LDL' of A
            LT = L';

            %Perform a forward substitution
            y = zeros(m+1,1);
            y(1) = r(1)/L(1,1);
            for k=2:m+1
                y(k) = ( r(k) - L(k,1:k-1)*y(1:k-1) )/L(k,k);
            end

            %Solve the diagonal system
            z = zeros(m+1,1);
            z = y ./D;

            %Perform a backward substitution
            dx = zeros(m+1,1);
            dx(m+1) = z(m+1)/LT(m+1,m+1);
            for k=m:-1:1
                dx(k) = ( z(k) - LT(k,k+1:m+1)*dx(k+1:m+1) )/LT(k,k);
            end

            %Update the minimizer
            xnp1 = xn + alpha*dx(1:m);

            %Compute change in the minimizer to check for convergence. 
            res = norm(xnp1 - xn)./norm(xnp1);
            
            %Update the iteration count since we have completed one more
            %LDL'factorisation
            iterCount = iterCount + 1;

            %Compute the value of the total cost function
            val(iterCount) = gamma*xn'*sigma*xn - (1 - gamma)*pbar'*xn - t.*sum(log(xn)); 
            xn = xnp1;
            
        end
    end

    %Compute the portfolio return
    ret = pbar'*xn;

    %Compute the risk
    risk = xn'*sigma*xn;
end