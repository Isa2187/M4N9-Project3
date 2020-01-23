%Define a function which computes the sine and cosine values for a Givens
%rotation in which we want to set one value to zero.
%Note that this algorithm is not my own, and is from Section 5.1.8 in the
%book by Golub and Van Loan.
function[c,s] = givensMatrix(a,b)
    if b == 0
        c = 1;
        s = 0;
    else
        if abs(b) > abs(a)
            r = -a/b;
            s = 1/sqrt(1+r^2);
            c = s*r;
        else
            r = -b/a;
            c = 1/sqrt(1+r^2);
            s = c*r;
        end
    end
end