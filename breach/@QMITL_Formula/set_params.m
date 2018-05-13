function phi = set_params(phi, param_name, param_value)
%SET_PARAMS does ...
% 
% Synopsis: phi = set_params(phi, param_name, param_value)
% 
% Inputs:
%  - phi         : QMITL formula 
%  - param_name  : name of the parameter to change (string). Can be
%                  obtained from get_param (in particular for time
%                  interval bounds) 
%  - param_value : new value of the parameter
% 
% Ouput:
%  - phi: updated QMITL_Formula
%

st = struct;
[phi, success] = find_and_set(phi, param_name, param_value);

if(~success)
    fprintf( 'warning: parameter not found, not set ')
end

end

function [phi, success] = find_and_set(phi, param_name, param_value)

success = 0;

switch(phi.type)
    
    case 'predicate'
        if isfield(phi.params, param_name)
            success=1;
            phi.params = setfield(phi.params,param_name,param_value);
        end
        
        phi.evalfn = @(mode,traj,t,params) feval('generic_predicate',mode,traj,t,params);
        
    case 'not'
        [phi.phi, success] = find_and_set(phi.phi, param_name, param_value);
        
    case 'or'
        [phi.phi1, success] = find_and_set(phi.phi1, param_name, param_value);
        if(~success)
            [phi.phi2, success] = find_and_set(phi.phi2, param_name, param_value);
        end
        
    case 'and'
        [phi.phi1, success] = find_and_set(phi.phi1, param_name, param_value);
        if(~success)
            [phi.phi2, success] = find_and_set(phi.phi2, param_name, param_value);
        end
        
    case '=>'
        [phi.phi1, success] = find_and_set(phi.phi1, param_name, param_value);
        if(~success)
            [phi.phi2, success] = find_and_set(phi.phi2, param_name, param_value);
        end
        
    case 'always'
        [phi.phi, success] = find_and_set(phi.phi, param_name, param_value);
        
    case 'eventually'
        [phi.phi, success] = find_and_set(phi.phi, param_name, param_value);
        
    case 'until'
        [phi.phi1, success] = find_and_set(phi.phi1, param_name, param_value);
        if(~success)
            [phi.phi2, success] = find_and_set(phi.phi2, param_name, param_value);
        end
end

end

function st1 = ConcatStruct(st1,st2)

listfield2 = fieldnames(st2);
fieldvalues2 = struct2cell(st2);

for i = 1:numel(listfield2)
    st1 = setfield(st1, listfield2{i}, fieldvalues2{i});
end

end
