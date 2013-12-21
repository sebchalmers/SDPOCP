function [tmp,res] = split(poly,var)

poly = expand(poly);
order0 = subs(poly,var,sym(0));
tmp = poly - order0;

% keyboard
order1 = subs(simplify(expand(tmp/var)),var,sym(0));

tmp = tmp - order1*var;

res = order0 + order1*var;
tmp = simplify(expand(tmp/var^2));


