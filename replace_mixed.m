function [poly,newVar,status] = replace_mixed(poly,newVar)

status = 0;
for k = 2:length(newVar)
    [tmp,res] = split0(poly,newVar(k));
    tmpVar = symvar(tmp);
    for j = 1:length(tmpVar)
        [tmp1,res1] = split0(tmp,tmpVar(j));
        tmpVar1 = symvar(tmp1);
        if ~isempty(tmpVar1)
            eval(['new_variable = sym(''', char(newVar(k)),'_',char(tmpVar(j)),''');'])
            eval(['poly = feval(symengine,''subsex'',expand(poly),''',char([char(newVar(k)),'*',char(tmpVar(j)),'=',char(new_variable)]),''');'])
            newVar(end+1) = new_variable;
            status = 1;
            return
        end
    end
end