function v = ilspencgetbk( b, k )
%ILSPENCGETBK Returning matrix of coeficients for k-th parameter.
%   Additional layer for testing repesentation of parameters' linear 
%   dependencies. Because of MATLAB's lazy copying this shouldn't be
%   significant performance overhead.

global dataModel;
switch(dataModel)
    case '3D'
        v = intval(b(:,k));
    case 'cell'
        % obtaining b-vector dimension and other meta-data
        m = b{1}(1);
        par = b{1}(2);
        % par == 0: normal
        % par == 1: intval values, mainly because monotonicity check func
        switch(par)
            case 0
            v = zeros(m,1);

            bk = b{k+1};
            for i = 1:length(bk(1,:))
                column = bk(:,i); 
                v(column(1)) = column(2);
            end
            
            case 1
            v = intval(zeros(m,1));

            bk = b{k+1};
            for i = 1:length(bk(1,:))
                column = bk(:,i); 
                v(column(1).inf) = column(2);
            end 
            
            otherwise
                % rise error
        end
        
        v = intval(v);
        
end

end

