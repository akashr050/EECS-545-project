A = zeros(9);
A(1, 2) = 1;
A(1, 4) = 1;
A(2, 3) = 1;
A(2, 5) = 1;
A(3, 6) = 1;
A(4, 5) = 1;
A(4, 7) = 1;
A(5, 6) = 1;
A(5, 8) = 1;
A(6, 9) = 1;
A(7, 8) = 1;
A(8, 9) = 1;
A = A+A';

JT_UGMS(H, X, rho_best, th_best);