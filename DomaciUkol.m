global dataModel;
dataModel = '3D';
A = zeros(2,2,3)
A(1,1,1) = 1
A(1,1,2) = 1
A(1,2,2) = -1
A(2,1,2) = -1
A(2,2,2) = 1
A(2,2,3) = 1

b = zeros(2,3)
b(1,1) = 10
b(2,1) = 5

p = [infsup(900,1100),infsup(900,1100),infsup(900,1100)]

x1 = ilspencbauerskeel(A, b, p)
x2 = ilspenchbr(A, b, p)
x3 = ilspencresidual(A, b, p, 'RUMP')