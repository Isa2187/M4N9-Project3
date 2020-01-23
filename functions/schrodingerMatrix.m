%Using the template provided, define a function which computes the
%pentadiagonal matrix L for any given value of N.
%The function outputs the pentadiagonal matrix, L, and the permutation
%matrix, P, such that L = P'AP.
function [L,P] = schrodingerMatrix(N)

    %Establish grid
    dx = 2*pi/N; 
    x = (0:dx:2*pi-dx)';

    %Sets potential parameters
    V_0 = 200;
    gamma = 10;

    %Set up D matrix
    D = -diag(2*ones(N,1)) + diag(ones(N-1,1),1) + diag(ones(N-1,1),-1);
    D(N,1) = 1;
    D(1,N) = 1;
    D = D./dx^2;

    %Evaluate potential and form discrete Schrodinger operator
    V = -V_0*exp(-gamma*(x - pi).^2);
    A = -D + diag(V);

    %Set up permutation matrix P
    P = zeros(N);
    for i = 1:2:N
        k = (i-1)/2 + 1;
        j = N - k + 1;
        P(k,i) = 1;
        P(j,i+1) = 1;
    end

    %Obtain pentadiagonal matrix
    L = P'*A*P;
end