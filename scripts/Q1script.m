%Form L and P for when N = 1024
[L,P] = schrodingerMatrix(1024);

%Form T and Q using L
[T,Q] = tridiagMatrix(L);

%Find the eigenvalues of T, setting tol = 1e-10
[eValues,orderedEValues] = eigenvalueFinder(T,1e-10);

%Plot the ordered eigenvalues against their index
figure()
plot(orderedEValues)
title('Graph Of E_{n} Against n')
xlabel('Eigenvalue index, n')
ylabel('Eigenvalue, E_{n}')
xlim([0 1024])

%Find the eigenvectors of T
eVectors = eigenvectorFinder(T,orderedEValues,1e-8);

%Using P and Q, find the eigenvectors of A from the eigenvectors of T
Psi = P * Q * eVectors;

%Plot the wave functions for n=1, 2, 3, 4 and 10
dx = 2*pi/1024; 
x = (0:dx:2*pi-dx)';
figure()
plot(x,Psi(:,1).*Psi(:,1))
hold on
plot(x,Psi(:,2).*Psi(:,2))
plot(x,Psi(:,3).*Psi(:,3))
plot(x,Psi(:,4).*Psi(:,4))
plot(x,Psi(:,10).*Psi(:,10))
title('Graph Of |\Psi_{n}|^{2} Against x For n=1, 2, 3, 4 and 10')
xlabel('Position, x')
ylabel('|\Psi_{n}|^{2}')
legend('n=1','n=2','n=3','n=4','n=10')
xlim([0 2*pi])
hold off

%Time how the time taken to form T using Givens rotations changes with N
NValuesQ12 = (1000:100:2500)';
timesQ12 = zeros(length(NValuesQ12),1);

for i = 1:length(NValuesQ12)
    N = NValuesQ12(i);
    L = schrodingerMatrix(N);
    
    tic;
    T = tridiagMatrix(L);
    timesQ12(i) = toc;
end

coeffQ12 = polyfit(log(NValuesQ12),log(timesQ12),1)

figure()
loglog(NValuesQ12,timesQ12)
hold on
loglog(NValuesQ12,exp(coeffQ12(2)).*NValuesQ12.^coeffQ12(1))
title('log-log Graph Of Time Against N When Forming T From L')
xlabel('Number of equispaced points, N')
ylabel('Time taken to form T from L, t')
legend('Recorded Data','t = e^c N^k')
hold off

%Time how the time taken to find the QR factorisation of T using Givens
%rotations changes with N
NValuesQ13 = (1000:100:2500)';
timesQ13 = zeros(length(NValuesQ13),1);

for i = 1:length(NValuesQ13)
    N = NValuesQ13(i);
    L = schrodingerMatrix(N);
    T = tridiagMatrix(L);
    tic;
    [Q,R] = tridiagQR(T);
    timesQ13(i) = toc;
end

coeffQ13 = polyfit(log(NValuesQ13),log(timesQ13),1)

figure()
loglog(NValuesQ13,timesQ13)
hold on
loglog(NValuesQ13,exp(coeffQ13(2)).*NValuesQ13.^coeffQ13(1))
title('log-log Graph Of Time Against N When Computing T = QR')
xlabel('Number of equispaced points, N')
ylabel('Time taken to find QR factorisation, t')
legend('Recorded Data','t = e^c N^k')
hold off

%Time how the time taken to solve Tx = b for tridiagonal T changes with N
NValuesQ14 = (1000:100:2500)';
timesQ14 = zeros(length(NValuesQ14),1);

for i = 1:length(NValuesQ14)
    N = NValuesQ14(i);
    L = schrodingerMatrix(N);
    T = tridiagMatrix(L);    
    
    c = diag(T,-1);
    d = diag(T);
    e = diag(T,1);
    
    for j = 1:10000
        v = rand(N,1);
        v = v ./ norm(v);
        tic;
        v = tridiagSolver(c,d,e,v);
        timesQ14(i) = timesQ14(i) + toc;
    end
end

timesQ14 = timesQ14./10000;

coeffQ14 = polyfit(log(NValuesQ14),log(timesQ14),1)

figure()
loglog(NValuesQ14,timesQ14)
hold on
loglog(NValuesQ14,exp(coeffQ14(2)).*NValuesQ14.^coeffQ14(1))
title('Graph Of Average Time To Solve Tx = b Against N Using tridiagSolver')
xlabel('Number of equispaced points, N')
ylabel('Time taken to solve Tx = b, t')
legend('Recorded Data','t = e^c N^k')
hold off