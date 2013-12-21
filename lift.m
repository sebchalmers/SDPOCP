function [poly,newVar] = lift(poly,var)

newVar = sym(1);

for j = 1:length(var)
    eval(['var1 = sym(''', char(var(j)),'_0'');'])
    
    newVar(end+1) = var1;

    poly = subs(expand(poly),var(j),var1);
end

tmpVar = newVar(2:end);

for j = 1:length(var)
    
    var1 = tmpVar(j);
    
    % Replace var1^2k
    for k = 1:100
        
        eval(['var2 = sym(''', char(var(j)),'_', num2str(2*k),''');'])

        sym(var2,'real');

        newVar(end+1) = var2;
        
        [poly,maxDeg] = replace(expand(poly),var1,var2);

        var1 = var2;
        
        if maxDeg <= 1
            break
        end
    end
end


% Replace mixed terms
status = 1;
for k = 1:100
    [poly,newVar,status] = replace_mixed(poly,newVar);
    if status ==0
        break
    end
end
    
% check = 1-isempty(symvar(jacobian(jacobian(poly,newVar),newVar)))
    