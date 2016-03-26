function v = ilspencresidual( A, b, p, option )
%ILSPENCRESIDUAL Solving parametric equations in residual form.

% Inicialization of general variables.
v = NaN;
% x-asterisk from Theorem 4.
x = inv(ilspencmatrixcenter(A,p))*ilspencbcenter(b ,p);

[Ares, bres] = ilspencresidualform(A, b, p);

if (nargin ~= 4)
    %ERROR
    return;
else
    switch option
        case 'RUMP'
         	y = ilspencrump(Ares, bres, 15);
        case 'SKALNA'
            y = ilspencskalna(Ares, bres);
        otherwise
            % invalid option
            return;
    end
end

v = x + y;
    
end

