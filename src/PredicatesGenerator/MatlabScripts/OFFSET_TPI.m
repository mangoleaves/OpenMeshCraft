syms xa ya za
syms xb yb zb
syms xc yc zc
syms xo yo zo
syms xp yp zp
syms xq yq zq
syms xr yr zr
syms xs ys zs
syms xt yt zt

assume([xa ya za xb yb zb xc yc zc], 'real');
assume([xo yo zo xp yp zp xq yq zq], 'real');
assume([xr yr zr xs ys zs xt yt zt], 'real');

a = [xa ya za];
b = [xb yb zb];
c = [xc yc zc];
o = [xo yo zo];
p = [xp yp zp];
q = [xq yq zq];
r = [xr yr zr];
s = [xs ys zs];
t = [xt yt zt];

% nor_abc = cross(b - a, c - a);
nor_opq = cross(p - o, q - o);
nor_rst = cross(s - r, t - r);

% We want to find the intersection point with the below formulation

% syms u v
% I = a + u * (b - a) + v * (c - a);

% If we substitute the intersection points into the plane equation,
% we will obtain the following system of linear equations.
% L (u v)^T = m

m0 = dot(o - a, nor_opq);
m1 = dot(r - a, nor_rst);
% m = [m0; m1];

L00 = dot(nor_opq, b - a);
L01 = dot(nor_opq, c - a);
L10 = dot(nor_rst, b - a);
L11 = dot(nor_rst, c - a);
% L = [L00 L01; L10 L11];

% Now solve L (u v)^T = m by Cramer's Law

% L_sub0 = [m0 L01; m1 L11];
% L_sub1 = [L00 m0; L10 m1];

% calculate determinant by hand
det_L = (L00 - L11) * (L01 - L10);
det_L_sub0 = (m0 - L11) * (L01 - m1);
det_L_sub1 = (L00 - m1) * (m0 - L10);

% the solution
% u_sol = det_L_sub0 / det_L;
% v_sol = det_L_sub1 / det_L;

% formulate as base + offset
base = o;
lambda_u = det_L_sub0 * (b - a);
lambda_v = det_L_sub1 * (c - a);
lambda = lambda_u + lambda_v
d = det_L