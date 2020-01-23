%Define a function which solves a tridiagonal system Tx = v. 
%The inputs a,b and c represent the subdiagonal, main diagonal and
%superdiagonal of T, and v is the right-hand side.
%The output x is the solution to Tx = v.
%Note that this algorithm is not my own - it is known as the Thomas
%algorithm and is used to solve tridiagonal linear systems (see link below)
%https://www.quantstart.com/articles/Tridiagonal-Matrix-Solver-via-Thomas-Algorithm
function x = tridiagSolver(a,b,c,v)
    
    n = size(b,1);    
    c1 = zeros(n-1,1);
    d1 = zeros(n,1);
    x = zeros(n,1);
    
    for i = 1:n-1
        
        if i == 1
            c1(i) = c(i)/b(i);
        else
            c1(i) = c(i)/( b(i) - c1(i-1)*a(i-1) );
        end
    end
    
    for i = 1:n-1
        if i == 1
            d1(i) = v(i)/b(i);
        else
            d1(i) = ( v(i) - d1(i-1)*a(i-1) )/(b(i) - c1(i-1)*a(i-1));
        end
    end
    
    d1(n) = ( v(n) - d1(n-1)*a(n-1) )/( b(n) - c1(n-1)*a(n-1) );
    
    for i = n:-1:1
        if i == n
            x(i) = d1(i);
        else
            x(i) = d1(i) - c1(i)*x(i+1);
        end
    end
end