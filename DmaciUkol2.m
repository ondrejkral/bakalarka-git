delta = 1;
X(1) = infsup(2-delta,2+delta);
X(2) = infsup(2-delta,2+delta);

%% natural extension
y1 = (X(1)^3 - X(1)*X(2))/(X(1)^2 + X(2)^2)

%% mean value form

a = mid(X);
Xg = gradientinit(X);
DFX = (Xg(1)^3 - Xg(1)*Xg(2))/(Xg(1)^2 + Xg(2)^2);
fa = (a(1)^3 - a(1)*a(2))/(a(1)^2 + a(2)^2);
y2 = fa + DFX.dx*(X - a)'

%% slope form

a = mid(X);
Xs = slopeinit(a,X);
fa = (a(1)^3 - a(1)*a(2))/(a(1)^2 + a(2)^2);
S = (Xs(1)^3 - Xs(1)*Xs(2))/(Xs(1)^2 + Xs(2)^2);
y3 = fa + S.s*(X - a)'