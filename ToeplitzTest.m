function ToeplitzTest
clc;
format long infsup;
rng(23101993, 'twister');
global dataModel;
dataModel = 'cell';
% parameters
iterations = 50;
%coefmagnitude = 0.1;
%radius = 0.001;
%matrixdim = 50; % N x N

fileID = fopen('C:\Users\Ondra\Documents\MATLAB\test20b.txt','a');
fprintf(fileID,dataModel);
fclose(fileID);

matrixdim = 1;
for l = 1:2
     matrixdim = matrixdim*10;
     radius = 0.01;
    for k = 1:3
        coefmagnitude = 1;
        for c = 1:3
            fileID = fopen('C:\Users\Ondra\Documents\MATLAB\test20b.txt','a');
            %disp('Magnitude of coeficients');
            %fprintf(fileID,'Magnitude of coefs: ');
            %fprintf(fileID,'%.4f\n',coefmagnitude);
            %coefmagnitude

            %disp('Radius: ');
            %fprintf(fileID,'Radius: ');
            %fprintf(fileID,'%d\n',radius);
            %radius

            %%disp('Toeplitz matrix dimension (square): ');
            %fprintf(fileID,'MatrixDim: ');
            %fprintf(fileID,'%d\n',matrixdim);
            %matrixdim
            %disp(' ')

            t = zeros(1,6); r = zeros(1,6);
            skips = zeros(1,6);
            solution = zeros(matrixdim,6);
            for i = 1:iterations
                disp(i);
                warning off;
                %disp('Toeplitz system creation');
                [A, b, p] = toeplitzsystem(coefmagnitude,radius,matrixdim);
                warning on;
                if (~verifyinvert(A, b, p))
                    skips(1) = skips(1) + 1;
                    skips(2) = skips(2) + 1;
                    skips(3) = skips(3) + 1;
                    skips(4) = skips(4) + 1;
                else
                    %disp('Bauer-Skeel enclousure');
                    tic;
                    x1 = ilspencbauerskeel(A, b, p);
                    t(1) = t(1) + toc;
                    r(1) = r(1) + avgradius(x1);
                    solution(:,1) = solution(:,1) + mid(x1);

                    %disp('Bauer-Skeel refinement');
                    tic;
                    r1 = ilspencbsref(A, b, p, x1);
                    t(2) = t(2) + toc;
                    r(2) = r(2) + avgradius(r1);
                    solution(:,2) = solution(:,2) + mid(r1);
                    refcheck(x1,r1, 'bauer-skeel');

                    %disp('Hans-Bliek-Rohn enclosure');
                    tic;
                    x2 = ilspenchbr(A, b, p);
                    t(3) = t(3) + toc;
                    r(3) = r(3) + avgradius(x2);
                    solution(:,3) = solution(:,3) + mid(x2);

                    %disp('Hans-Bliek-Rohn refinement');
                    tic;
                    r2 = ilspenchbrref(A, b, p, x2);
                    t(4) = t(4) + toc;
                    r(4) = r(4) + avgradius(r2);
                    solution(:,4) = solution(:,4) + mid(r2);
                    refcheck(x2,r2, 'hans-bliek-rohn');
                end
                
                if (~verifyinvert2(A,b,p))
                    skips(5) = skips(5) + 1;
                else
                    %disp('Residual method with RUMP');
                    tic
                    x3 = ilspencresidual(A, b, p, 'RUMP');
                    t(5) = t(5) + toc;
                    r(5) = r(5) + avgradius(x3);
                    solution(:,5) = solution(:,5) + mid(x3);
                end
                
                if (~verifyinvert3(A,b,p))
                    skips(6) = skips(6) + 1;
                else
                    %disp('Residual method with SKALNA');
                    tic;
                    x4 = ilspencresidual(A, b, p, 'SKALNA');
                    t(6) = t(6) + toc;
                    r(6) = r(6) + avgradius(x4);
                    solution(:,6) = solution(:,6) + mid(x4);
                end
            end

            %disp('Time:');
            %fprintf(fileID,'Time: \n');
            resulttime = zeros(1,6);
            for z = 1:6
                resulttime(z) = t(z)/(iterations - skips(z));
            end
            %fprintf(fileID,'%.8f\n',resulttime);
            %t/iterations

            %disp('Radius:')
            %fprintf(fileID,'Radius:\n ');
            resultradius = zeros(1,6);
            for z = 1:6
                resultradius(z) =r(z)/(iterations - skips(z));
            end
            %fprintf(fileID,'%.8f\n', resultradius);
            %r/iterations
            
            %disp('AproxSol:');
            %fprintf(fileID,'AproxSol:\n ');
            aproxsol = zeros(matrixdim,6);
            for z = 1:6
                aproxsol(:,z) = solution(:,z)/(iterations - skips(z));
            end
            %fprintf(fileID,'%.8f\n',aproxsol);

            %disp('Skips (not convenient M):')
            %fprintf(fileID,'Skips: ');
            %fprintf(fileID,'%d\n', skips);
            %fprintf(fileID,'--------------------\n');
            %skips
            
            for a = 1:6
                fprintf(fileID,'%.8f %.8f %d %d %d %d\n',resultradius(a),resulttime(a),matrixdim, coefmagnitude,a,skips(a));
            end
            fclose(fileID);
            coefmagnitude = coefmagnitude*10;
        end
        radius = radius*10;
    end
end
end

