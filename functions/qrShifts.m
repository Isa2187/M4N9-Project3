%Define a function qrShifts, which computes one eigenvalue of the input
%matrix T using the QR algorithm with shifts.
%The input tol is used as a measure of when to stop iterating.
%Note that the function leaves an eigenvalue in the bottom-right most 
%entry of T when T is outputted.
function T = qrShifts(T,tol)   
    
    %Set m to be the size of T
    m = size(T,1);
    
    %Assign mu to be the bottom-right most entry of T
    mu = T(m,m);
    
    %Initialise tZero to be the entry directly to the left of the 
    %bottom-right most entry of T, so that when this gets less than tol we
    %stop iterating and can readily obtain the desired eigenvalue
    tmm1 = abs(T(m,m-1));

    %Repeat the QR algorithm with shifts until we find the eigenvalue
    %within some measure of accuracy dependent on tol
    while(tmm1 > tol)

        %Shift the matrix T i.e. T - mu*I
        for j=1:m
            T(j,j) = T(j,j) - mu;
        end
        
        %Otain the QR factorisation of the shifted matrix i.e. T-mu*I =QR
        [Q,R] = tridiagQR(T);

        %Update T such that T = RQ + mu*I
        T = R*Q;
        for j=1:m
            T(j,j) = T(j,j) + mu;
        end

        %Apply the Wilkinson Shift to determine the next value of mu - this
        %algorithm is from lectures and I do not claim it to be my own
        am1 = T(m-1,m-1);
        am = T(m,m);
        bm1 = T(m,m-1);
        delta = (am1 - am)/2;
        temp = abs(delta) + sqrt(delta^2 + bm1^2);

        if(delta == 0)
            mu = am - bm1^2/temp;
        else
            mu = am - sign(delta)*bm1^2/temp;
        end

        tmm1 = abs(T(m,m-1));
    end
end