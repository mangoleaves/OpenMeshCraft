syms xa ya za
syms xb yb zb
syms xp yp zp
syms xq yq zq

n = (xa-xp)*(yq-yp) - (ya-yp)*(xq-xp);
d = (yb-ya)*(xq-xp) - (xb-xa)*(yq-yp);

lambda_x = xa*d + n*(xb - xa);
lambda_y = ya*d + n*(yb - ya);
lambda_z = za*d + n*(zb - za);

simplify(lambda_x)
simplify(lambda_y)
simplify(lambda_z)
