function V = ilspencgetak( A, k )
%ILSPENCGETAK Returning matrix of coeficients for k-th parameter.
%   Additional layer for testing repesentation of parameters' linear 
%   dependencies. Because of MATLAB's lazy copying this shouldn't be
%   significant performance overhead.

global dataModel;
switch(dataModel)
    case '3D'
        V = A(:,:,k);
    case 'cell'
        % obtaining matrix dimension and other meta-data
        m = A{1}(1);
        n = A{1}(2);
        
        par = A{1}(3);
        % par == 0: normal
        % par == 1: intval values, mainly because monotonicity check func
        
        switch(par)
            case 0
                V = zeros(m,n);
                Ak = A{k+1};
                for i = 1:length(Ak(1,:))
                    column = Ak(:,i); 
                    V(column(1),column(2)) = column(3);
                end
            case 1
                V = intval(zeros(m,n));
                Ak = A{k+1};
                for i = 1:length(Ak(1,:))
                    column = Ak(:,i); 
                    V(column(1).inf,column(2).inf) = column(3);
                end
            otherwise
                % rise error
        end
end
end

