function [tmp,res] = split0(poly,var)

order0 = subs(poly,var,sym(0));

tmp = poly - order0;

res = order0;
tmp = simplify(expand(tmp/var));


