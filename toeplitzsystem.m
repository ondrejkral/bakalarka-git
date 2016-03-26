function [ A, b, p ] = toeplitzsystem( coefmagnitude, radius, matrixmagnitude)
% Number of parameters needed.
parametersCount = 2*matrixmagnitude - 1;

% Random center of parameters. First parameter is diagonal one, so we 
% made it bigger to avoid irregularity.
random = -10 + 20*rand(parametersCount,1);
random(1) = random(1) + 10*matrixmagnitude;
% parameters vector
p = midrad(random,radius);

% random vector
% randomBCenter = -coefmagnitude + (2*coefmagnitude)*rand(matrixmagnitude,1);
randomBCenter = -10 + 20*rand(1,dimension);

global dataModel;

switch(dataModel)
    case '3D'
        
        % matrix A representation
        A = zeros(matrixmagnitude,matrixmagnitude,2*matrixmagnitude - 1);
        
        c = zeros(matrixmagnitude,1); c(1) = 1;
        r = zeros(matrixmagnitude,1);

        for i = (matrixmagnitude+1):(2*matrixmagnitude - 1)
            c = circshift(c,1);
            A(:,:,i) = toeplitz(c,r);
        end

        c = zeros(matrixmagnitude,1); c(1) = 1;
        r = zeros(matrixmagnitude,1); r(1) = 1;
        for i = 1:matrixmagnitude
            A(:,:,i) = toeplitz(c,r);
            c(1) = 0;
            r = circshift(r,1);
        end

        % vector b representation
        b = zeros(matrixmagnitude,2*matrixmagnitude-1);
        b(:,1) = randomBCenter;
               
    case 'cell'
        
        % matrix A representation
        A{1} = [ matrixmagnitude; matrixmagnitude; 0;];
        b{1} = [matrixmagnitude; 0;];
        
        for i = 1:matrixmagnitude
            A{i+1} = zeros(3,matrixmagnitude - i + 1);
            for j = 1:(matrixmagnitude - i + 1)
                A{i+1}(:,j)=[j; (j+i-1); 1;];
            end
        end
        for i = 2:matrixmagnitude
            A{matrixmagnitude + i} = zeros(3,matrixmagnitude - i + 1);
            for j = 1:(matrixmagnitude - i + 1)
                A{matrixmagnitude + i}(:,j)=[(j+i-1); j; 1;];
            end
        end
        
        % vector b representation
        randomvector = randomBCenter;
        for j = 1:matrixmagnitude 
            b{2}(:,j) =[ j; randomvector(j);] ;
        end
        for i = 3:(2*matrixmagnitude+1)
            b{i} = [1; 0;];
        end
            
end                
end