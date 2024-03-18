syms xa ya za
syms xb yb zb
syms xp yp zp
syms xq yq zq

assume([xa ya za xb yb zb xp yp zp xq yq zq], 'real');

% line ab
a = [xa ya za];
b = [xb yb zb];
% line pq
p = [xp yp zp];
q = [xq yq zq];

% let "I = a + t (b - a)"
% we have "(I - p) cross (q - p) = 0"
% transform to "((a - p) + t (b - a)) cross (q - p) = 0"

% calculate "t = n / d"
% n = cross (a - p, q - p)
n = cross(a-p, q-p)
% d = - cross (b - a, q - p)
d = - cross(b-a, q-p)

lambda = n * (b - a)

% base is a
base = a;

% now "I = base + lambda / d"
I = base + lambda / d