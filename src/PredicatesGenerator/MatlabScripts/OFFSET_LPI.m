syms xa ya za
syms xb yb zb
syms xo yo zo
syms xp yp zp
syms xq yq zq

assume([xa ya za xb yb zb xo yo zo xp yp zp xq yq zq], 'real');

% line ab
a = [xa ya za];
b = [xb yb zb];
% plane opq
o = [xo yo zo];
p = [xp yp zp];
q = [xq yq zq];

% intersection "I = a + t (b - a)"
% we have "(I- o) dot ((p - o) cross (q - o)) = 0"
ba = b - a;
po = p - o;
qo = q - o;
ao = a - o;
% normal of plane opq
normal_opq = cross(po, qo)

d = dot(ba, normal_opq)

n = - dot(ao, normal_opq)

base = a;

lambda = n * ba

% now "I = base + lambda / d"