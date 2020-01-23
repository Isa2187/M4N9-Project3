%Define a function that takes as inputs a tridiagonal matrix, T, the
%eigenvalues of T, eValues, and a small number, diff.
%It applies inverse iteration to T to obtain its eigenvectors, which are
%outputted in the matrix eVectors
function eVectors = eigenvectorFinder(T,eValues,diff)
    
    %Set n to be the size of T
    n = size(T,1);
    
    %Initialise the matrix eVectors, whose columns will be the eigenvectors
    %of T
    eVectors = zeros(n,n);
    
    %If T is a real number, set the eigenvector to be 1
    if n == 1
        eVectors = 1;
    else
        
        for i = 1:n

            %Initialise Ti and lambda to be the current eigenvalue
            Ti = T;
            lambda = eValues(i);
            
            %Set mu to be lambda + diff
            mu = lambda + diff;

            %Set our convergence criteria, which we will use as a measure
            %to know when to stop applying inverse iteration below
            if i == 1
                convergenceCriteria = abs((mu - lambda)/(mu - eValues(i+1)));
            elseif i == n
                convergenceCriteria = abs((mu - lambda)/(mu - eValues(i-1)));
            else
                convergenceCriteria = max( abs((mu - lambda)/(mu - eValues(i-1))) , abs((mu - lambda)/(mu - eValues(i+1))) );
            end

            %Shift the matrix using mu
            for j = 1:n
                Ti(j,j) = Ti(j,j) - mu;
            end

            %Initialise v to be a random unit vector i.e. v0 in the inverse
            %iteration algorithm
            v = rand(n,1);
            v = v ./ norm(v);
            
            %Initialise lambdaCvgt and count, where lambdaCvgt is assessed
            %against convergenceCriteria after each iteration ant count is
            %simply the number of iterations that have been performed
            lambdaCvgt = 1000;
            count = 0;

            %Set c to be the subdiagonal of Ti, d to be the main diagonal
            %of Ti and e to be the superdiagonal of Ti
            c = diag(Ti,-1);
            d = diag(Ti);
            e = diag(Ti,1);

            while (lambdaCvgt > convergenceCriteria)
                
                %Use the function tridiagSolver to solve Ti w = v
                w = tridiagSolver(c,d,e,v);
                
                %Set vk to be the unit vector obtained from w
                v = w ./ norm(w);

                %Compute lambdaCvgt = abs(v'Tv - lambda)
                lambdaCvgt = abs(v'*T*v - lambda);
                
                %Update the iteration count
                count = count + 1;
            end
            
            %Store the eigenvector as a column of eVectors
            eVectors(:,i) = v;
        end
    end
end