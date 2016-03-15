function v = ilspenc( A, b, p, option )
%ILSPENC Parametric interval systems solver.

% ILSPENCGETAK, ILSPENCGETBK, ILSPENCMATRIXDIM, ILSPENCVECTORDIM
%   - reimplement these in case of testing new representation 
% of parameters' linear dependencies.

% Testing of conditions like spectral radius etc.
% Determining which algorithms or combination of them to use.

if (nargin < 3 || nargin > 4)
    % ERROR?
    v = NaN; 
    return;
elseif (nargin == 3)
    % no option method
else   
    switch option
        case 'TIGHTEST'
            % using all methods an intersecting results
            x1 = ilspencbauerskeel(A, b, p); x1 = ilspencbsref(A, b, p, x1);
            x2 = ilspenchbr(A, b, p); x2 = ilspenchbrref(A, b, p, x2);
            x3 = ilspencresidual(A, b, p, 'RUMP');
            
            v1 = intersect(x1,x2);
            v = intersect(v1,x3);
        otherwise
            % invalid option argument
    end
end
end

