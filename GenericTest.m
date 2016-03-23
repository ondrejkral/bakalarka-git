function GenericTest(model, matrixType)
clc;
format long infsup;
rng(23101993, 'twister');
global dataModel;
dataModel = model;
iterations = 10;
radiusSample = [0.05, 0.1, 0.5, 1];
coefMagMultSample = [5, 10, 15 20, 25];
matrixDimSample = [5, 10, 15, 20, 25 50, 100];

fileName = strcat('test-',datestr(datetime(),'dd-mm-yy-HH-MM'),'.txt');
fileLoc = 'C:\Users\ondre\Documents\MATLAB\';
fileAddr = strcat(fileLoc,fileName);

fileID = fopen(fileAddr,'a');

fprintf(fileID,model);
fclose(fileID);

for l = 1:length(matrixDimSample)
     matrixdim = matrixDimSample(l);
    for k = 1:length(coefMagMultSample)
        coefmagnitude = coefMagMultSample(k);
        for c = 1:length(radiusSample)
            radius = radiusSample(c);
            fileID = fopen(fileAddr,'a');
            
            t = zeros(1,6); r = zeros(1,6);
            skips = zeros(1,6);
            solution = zeros(matrixdim,6);
            
            for i = 1:iterations
                disp(i);
                warning off;
                switch(matrixType)
                    case 'Sym'
                        [A, b, p] = SymetricSystem(coefmagnitude, matrixdim, radius);
                    case 'Toeplitz'                       
                        %disp('Toeplitz system creation');
                        [A, b, p] = toeplitzsystem(coefmagnitude,radius,matrixdim);
                    otherwise
                        % invalid matrixType
                end
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
                    %refcheck(x1,r1, 'bauer-skeel');

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
                    %refcheck(x2,r2, 'hans-bliek-rohn');
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

            resulttime = zeros(1,6);
            for z = 1:6
                resulttime(z) = t(z)/(iterations - skips(z));
            end
            
            resultradius = zeros(1,6);
            
            for z = 1:6
                resultradius(z) =r(z)/(iterations - skips(z));
            end
    
            aproxsol = zeros(matrixdim,6);
            
            for z = 1:6
                aproxsol(:,z) = solution(:,z)/(iterations - skips(z));
            end
            
            for a = 1:6
                fprintf(fileID,'%.8f %.8f %d %d %d %d\n',resultradius(a),resulttime(a),matrixdim, coefmagnitude,a,skips(a));
            end
            fclose(fileID);
        end
    end
end
end

