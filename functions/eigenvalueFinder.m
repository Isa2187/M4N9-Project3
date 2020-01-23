%Define a function that takes a matrix T as an input, and returns a vector
%orderedEValues, which contains the eigenvalues of T found using the QR
%algorithm with shifts in ascending order.
%The input tol is a tolerance to be used when calling the function qrShifts
function [eValues,orderedEValues] = eigenvalueFinder(T,tol)
 
    %Set n to be the size of T, and eValues which is an be n-vector that
    %will store the n eigenvalues of T
    n = size(T,1);
    eValues = zeros(n,1);
    
    %Find the first n-1 eigenvalues
    for i = 1:n-1
        
        %Update T after applying the function qrShifts, which results in an
        %eigenvalue being in the bottom-right most entry of T
        T = qrShifts(T,tol);
        
        %Store this eigenvalue of T in eValues
        eValues(n-i+1) = T(n-i+1,n-i+1);
        
        %Delete the last row and column of T, so that the next eigenvalue
        %of T can be found on this slightly smaller matrix
        T(end,:) = [];
        T(:,end) = [];
    end
    
    %Eventually T will be a 1-by-1 matrix, so store this as the final
    %eigenvalue in orderedEValues
    eValues(1) = T;
    
    %Sort the eigenvalues of T in order of size, so they appear in
    %ascending order, in the n-vector orderedEValues
    orderedEValues = sort(eValues);
end