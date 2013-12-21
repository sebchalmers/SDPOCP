function [new,k] = replace(poly,var,newvar)

new = 0;

tmp = poly;
for k = 0:100
    
    [tmp,res] = split(expand(tmp),var);
    new = new + res*newvar^k;
    
    if tmp == 0
        break
    end
end

if (k == 100) && (tmp ~= 0)
    keyboard
    error('Polynomial of degree > 100')
end