%Define a function tridiagMatrix, which transforms the symmetric matrix L
%into a tridiagonal matrix T using Givens rotations. 
%The function also outputs the matrix Q, the matrix containing the combined
%Givens rotations such that T = Q'LQ
function [T,Q] = tridiagMatrix(L)
    
    %Set n to be the size of the matrix L
    n = size(L,1);
    
    %Initialise T to be L, and Q the n-by-n identity matrix
    T = L;
    Q = eye(n);
    
    %Iterate over these specific rows and columns, since we only require
    %that T is tridiagonal
    for i = 1:n-2 
        for j = i+2:n     
            
            %Find the entries we are performing a Givens rotation on
            a = T(i+1,i);
            b = T(j,i);

            %If b is sufficiently small, assume it is equal to zero and 
            %skip this iteration, as there is no need to perform a Givens 
            %rotation on such an entry
            if abs(b) < 1e-10
                continue
            end

            %Store the sine and cosine of this Givens rotation, and form
            %the matrix G with these values in their respective place
            [c,s] = givensMatrix(a,b);        
            G = [c s; -s c];        
            
            %Update the (i+1)-th and j-th rows and columns of T by
            %multiplying them by G' and G respectively
            T([i+1,j],:) = G' * T([i+1,j],:);
            T(:,[i+1,j]) = T(:,[i+1,j]) * G;
            
            %Update the (i+1)-th and j-th columns of Q by multiplying by G
            Q(:,[i+1,j]) = Q(:,[i+1,j]) * G;
        end
    end
end