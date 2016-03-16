function Test
global dataModel;
dataModel = '3D';
a0 = zeros(5); a2 = zeros(5); a3 = zeros(5); a4 = zeros(5); a5 = zeros(5);
a6 = zeros(5); a7 = zeros(5); a8 = zeros(5); a9 = zeros(5); a1 = zeros(5);

a1(1,1) = 1; a2(2,2) = 1; a3(3,3) = 1; a4(4,4) = 1; a5(5,5) = 1;
a6(1,1) = 1; a6(1,2) = -1; a6(2,2) = 1;  a6(2,1) = -1;
a7(2,2) = 1; a7(3,3) = 1; a7(3,2) = -1; a7(2,3) = -1;
a8(3,3) = 1; a8(4,3) = -1; a8(3,4) = -1; a8(4,4) = 1;
a9(4,4) = 1; a9(4,5) = -1; a9(5,4) = -1; a9(5,5) = 1;

b = zeros(5,10); b(1,1) = 10; b(3,1) = 10;  

p(1) = intval('1');
p(2:10) = infsup(0.99, 1.01);

% x = [infsup(7.0148,7.1671); infsup(4.1173,4.2463); infsup(5.3933,5.5158); infsup(2.1377,2.2260); infsup(1.0601,1.1217);];
% x = [infsup(6.,7.1667); infsup(4.1180,4.2456); infsup(5.3938,5.5158); infsup(2.1377,2.2260); infsup(1.0601,1.1217);];

array(:,:,1) = a0; array(:,:,2) = a1; array(:,:,3) = a2; array(:,:,4) = a3; array(:,:,5) = a4;
array(:,:,6) = a5; array(:,:,7) = a6; array(:,:,8) = a7; array(:,:,9) = a8; array(:,:,10) = a9;

format short infsup;
display('Skalna')
ilspencskalna2(array,b,p)
ilspenchbrref(array,b,p,ilspenchbr(array,b,p))
ilspencmono(array,b,p)
end