function [tmp,res] = split(poly,var)

order0 = subs(poly,var,0);

tmp = poly - order0;

order1 = subs(expand(simplify(tmp/var)),var,0);

tmp = tmp - order1*var;

res = order0 + order1*var;
tmp = simplify(expand(tmp/var^2));


