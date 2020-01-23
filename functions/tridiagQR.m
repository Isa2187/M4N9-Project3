%Define a function tridiagQR, which computes the QR factorisation of a
%tridiagonal matrix, T, such that T = QR
function [Q,R] = tridiagQR(T)

    %Set m to be the size of T
    m = size(T,1);
    
    %Initialise R to be T, and Q the n-by-n identity matrix
    Q = eye(m);
    R = T;

    %Since T is tridiagonal, we only have to update one entry in the first
    %m-1 columns
    for j = 1:m-1   

        %Find the entries we are performing a Givens rotation on
        a = R(j,j);
        b = R(j+1,j);
        
        %If b is sufficiently small, assume it is equal to zero and skip
        %this iteration, as there is no need to perform a Givens rotation
        %on such an entry
        if abs(b) < 1e-10
            continue
        end
        
        %Store the sine and cosine of this Givens rotation, and form the
        %matrix G with these values in their respective place
        [c,s] = givensMatrix(a,b);        
        G = [c s; -s c];       

        %Update the j-th and (j+1)-th rows of R and the and j-th and 
        %(j+1)-th columns Q by multiplying them by G' and G respectively
        R(j:j+1,:) = G' * R(j:j+1,:);
        Q(:,j:j+1) = Q(:,j:j+1) * G;
    end
end