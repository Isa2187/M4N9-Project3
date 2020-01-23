# M4N9-Project3
Third project of Computational Linear Algebra (M4N9) module taken in 4th year. (Grade = 78.3%)

All code was done in MATLAB, and further details of the task can be found in the folder Project_Files.

This project was composed of two separate components, the first being solving the time-independent Schrodinger equation as an eigenvector problem. Using Givens rotations, a pentadiagonal matrix was transformed into a triangular matrix, before implementing the QR algorithm on this to find the eigenvalues. Having obtained these, inverse interation was performed to obtain the eigenvectors.

Next, I considerred a portfolio optimisattion problem in financial mathematics. Subject to certain constraints, the aim was to minimise a linear function, which was synonymous with solving a linear system. This was first done using LDL' factorisation, before implementing the conjugate gradient method, and comparing the performance of each.
