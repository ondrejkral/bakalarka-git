function TestCell
rng(23101993, 'twister');
global dataModel;
format long infsup;
r = randi([0 10000],1,1000); 
for i = 1:1000
       
    rng(r(i), 'twister');
    dataModel = '3D';
    [A1, b1, p1] = toeplitzsystem(100,0.01,100);
    x1 = ilspencresidual(A1,b1,p1, 'RUMP');

    rng(r(i), 'twister');
    dataModel = 'cell';
    [A2, b2, p2] = toeplitzsystem(100,0.01,100);
    x2 = ilspencresidual(A2,b2,p2, 'RUMP');

    if ((~isnan(x1)) & (x1 ~= x2) )
        disp('FOUND ONE!')
        %fileID = fopen('C:\Users\Ondra\Documents\MATLAB\testcell2.txt','a');
        %fprintf(fileID,'%d\n',A1); fprintf(fileID,'%d\n',b1); fprintf(fileID,'%.8f\n',p1);
        %fprintf(fileID,'%d\n',A2); fprintf(fileID,'%d\n',b2); fprintf(fileID,'%.8f\n',p2);
        %fclose(fileID);
    end
end