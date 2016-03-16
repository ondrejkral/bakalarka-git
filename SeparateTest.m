function SeparateTest
global dataModel;
rng(23101993, 'twister');

dataModel = '3D';
[A, b, p] = toeplitzsystem(10,0.01,10);
x1 = ilspencbauerskeel(A,b,p)
rng(23101993, 'twister');

dataModel = 'cell';
[A, b, p] = toeplitzsystem(10,0.01,10);
x2 = ilspencbauerskeel(A,b,p)

in(x1,x2)
in(x2,x1)
end