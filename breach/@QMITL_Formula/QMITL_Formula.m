function [phi, phistruct] = QMITL_Formula(varargin)
%QMITL_FORMULA QMITL formula constructor, should be called without output.
%
%     QMITL_Formula(id) creates a new formula named id
%
%     | QMITL_Formula(id,'phi_expression') where phi obeys the following grammar
%
%      phi_expression := predicate_expression | unary_op phi_expression | '('phi_expression')' binary_op '('phi_expression')'
%
%      predicate_expression :=   matlab_expression = matlab_expression '|' threshold = num_value, max_true_value = num_value
%                              | matlab_expression comp matlab_expression
%
%      comp := > | >= | < | <=
%
%      unary_op := not | unary_op_temp | unary_op_temp_[ti,tf]
%
%      unary_op_temp :=  ev | alw | eventually | always
%
%      binary_op  := or | and | until | until_[matlab_expession, matlab_expression]
%
%
%     | QMITL_Formula(id,'unary_op',phi) where 'unary_op' is one of 'not',
%                'eventually' or 'always' (abreviations 'ev' or 'alw')
%
%     | QMITL_Formula(id,'unary_op',phi, interval) where 'unary_op' is
%                either 'eventually' or 'always'(abreviations 'ev' or 'alw')
%
%     | QMITL_Formula(id,'binary_op',phi1, phi2) where 'binary_op' is
%                either 'or','and' or 'until'
%
%     | QMITL_Formula(id,'until',phi1,interval, phi2)
%


% test if formula already exists

if(nargin==2)
    try % check copy operation
        st = varargin{2};
        phi = evalin('base', st);
        if isa(phi,'QMITL_Formula')
            phi.id = varargin{1};
            phistruct = struct(phi);
            assignin('base', phi.id, phi);
            return;
        end
    end
end

% OK, new formula or erroneous args, let's proceed

if(numel(varargin)==0)
    varargin{1} = 'phi';
end

phi.id = varargin{1};
phi.st = '';
phi.evalfn = [];
phi.interval = [0 inf];
phi.phi = [];
phi.phi1 = [];
phi.phi2 = [];
phi.phin = [];
phi.type = '';
phi.params = struct;
phi.params_interval = struct;
varargin = varargin(2:end);

switch numel(varargin)
    case 0
        phi.st = 'true';
        phi.evalfn = @true_formula;
        phi.type = 'predicate';
        
    otherwise
        if(numel(varargin)==1 && ischar(varargin{1}))
            varargin{1} = regexprep(varargin{1},'eventually', 'ev');
            varargin{1} = regexprep(varargin{1},'always','alw');
            varargin{1} = regexprep(varargin{1},'<=', '<');
            varargin{1} = regexprep(varargin{1},'>=', '>');
        end
        
        phi = QMITL_Parse(phi,varargin{:});
end

phi=class(phi, 'QMITL_Formula');
phistruct = struct(phi);
assignin('caller', phi.id, phi);
assignin('base', phi.id, phi);

end


function phi = QMITL_Parse(phi,varargin)
%QMITL_PARSE fills the field type, phi, phi1 and phi2 
%
% Synopsis: phi = QMITL_Parse(phi, phi_str)
%
%      OR : phi = QMITL_Parse(phi, unary_op, phi0)
%      OR : phi = QMITL_Parse(phi, unary_op2, interv, phi2)
%      OR : phi = QMITL_Parse(phi, binary_op, phi1, phi2)
%      OR : phi = QMITL_Parse(phi, 'until', phi1, interv, phi2)
%      OR : phi = QMITL_Parse(phi, 'andn', [phi1, phi2, ..., phin])
%
% Inputs:
%  - phi       : is the QMITL_Formula to create
%  - phi_str   : a string describing the formula. This string must follow
%                the grammar described in the QMITL_Formula documentation
%  - phi0      : a QMITL Formula
%  - phi1      : a QMITL Formula
%  - phi2      : a QMITL Formula
%  - phin      : a QMITL Formula
%  - unary_op  : is either 'not', 'ev', 'alw', 'eventually' or 'always'
%  - unary_op2 : is either 'ev', 'alw', 'eventually' or 'always'
%  - binary_op : is either 'or', 'and' or '=>'
%  - interv    : is an interval
%
% Output:
%  - phi : a QMITL Formula structure
%

switch(numel(varargin))
    
    case 1 % Here, the formula is defined by a string. We parse this string
        if ~ischar(varargin{1})
            error('QMITL_Formula:QMITL_Parse','Invalid formula');
        end
        st = varargin{1}; %NM: cannot be replaced by remove_parenthesis(varargin{1}); because, "((a) and (b))" fails (it becomes "a) and (b")
        
        
        % deals with true and false
        
        st = regexprep(st, 'true', 'inf>0');
        st = regexprep(st, 'false', 'inf<0');
        
        % test not
        
        tokens = regexp(st, '^\s*\<not\>(.+)','tokens');
        if ~isempty(tokens)
            st1 = remove_parenthesis(tokens{1}{1});
            phi1 = QMITL_Formula([phi.id '1__'],st1);
            phi = QMITL_Parse(phi, 'not', phi1);
            return
        end
        
        % test eventually
        
        tokens = regexp(st, '^\s*\<ev\>(.+)','tokens');
        if ~isempty(tokens)
            st1 = remove_parenthesis(tokens{1}{1});
            phi1 = QMITL_Formula([phi.id '1__'],st1);
            phi = QMITL_Parse(phi, 'eventually', phi1);
            return
        end
        
        % test eventually_[ti,tf]
        tokens = regexp(st, '^\s*\<ev_\[(.+?)\](.+)','tokens');
        if ~isempty(tokens)
            st1 = remove_parenthesis(tokens{1}{2});
            I = ['[' tokens{1}{1} ']'];
            phi1 = QMITL_Formula([phi.id '1__'],st1);
            phi = QMITL_Parse(phi,'eventually',I,phi1);
            return
        end
        
        % test evp
        tokens = regexp(st, '^\s*\<evp\>(.+)','tokens');
        if ~isempty(tokens)
            st1 = remove_parenthesis(tokens{1}{1});
            phi1 = QMITL_Formula([phi.id '1__'],st1);
            phi = QMITL_Parse(phi, 'evp', phi1);
            return
        end
        
        % test evp_[ti,tf]
        tokens = regexp(st, '^\s*\<evp_\[(.+?)\](.+)','tokens');
        if ~isempty(tokens)
            st1 = remove_parenthesis(tokens{1}{2});
            I = ['[' tokens{1}{1} ']'];
            phi1 = QMITL_Formula([phi.id '1__'],st1);
            phi = QMITL_Parse(phi,'evp',I,phi1);
            return
        end
        
        
        % test always
        tokens = regexp(st, '^\s*\<alw\>(.+)','tokens');
        if ~isempty(tokens)
            st1 = remove_parenthesis(tokens{1}{1});
            phi1 = QMITL_Formula([phi.id '1__'],st1);
            phi = QMITL_Parse(phi,'always', phi1);
            return
        end
        
        % test always_[ti,tf]
        tokens = regexp(st, '^\s*\<alw_\[(.+?)\](.+)','tokens');
        if ~isempty(tokens)
            st1 = remove_parenthesis(tokens{1}{2});
            I = ['[' tokens{1}{1} ']'];
            phi1 = QMITL_Formula([phi.id '1__'],st1);
            phi = QMITL_Parse(phi,'always',I,phi1);
            return
        end
        
        % test and
        [success, st1, st2] = parenthesisly_balanced_split(st, 'and');
        if success
            phi1 = QMITL_Formula([phi.id '1__'],st1);
            phi2 = QMITL_Formula([phi.id '2__'],st2);
            phi = QMITL_Parse(phi,'and', phi1, phi2);
            return
        end
        
        % test or
        [success, st1, st2] = parenthesisly_balanced_split(st, 'or');
        if success
            phi1 = QMITL_Formula([phi.id '1__'],st1);
            phi2 = QMITL_Formula([phi.id '2__'],st2);
            phi = QMITL_Parse(phi,'or', phi1, phi2);
            return
        end
        
        % test until
        [success, st1, st2] = parenthesisly_balanced_split(st, 'until');
        interval = '[0 inf]';
        if success
            phi1 = QMITL_Formula([phi.id '1__'],st1);
            phi2 = QMITL_Formula([phi.id '2__'],st2);
            phi = QMITL_Parse(phi,'until', phi1, interval, phi2);
            return
        end
        
        
        % test implies
        [success, st1, st2] = parenthesisly_balanced_split(st, '=>');
        if success
            phi1 = QMITL_Formula([phi.id '1__'],st1);
            phi2 = QMITL_Formula([phi.id '2__'],st2);
            phi = QMITL_Parse(phi,'=>', phi1, phi2);
            return
        end
        
        
        % test until_[ti,tf]
        [success, st1, st2, interval] = parenthesisly_balanced_split_interval(st, 'until');
        if success
            phi1 = QMITL_Formula([phi.id '1__'],st1);
            phi2 = QMITL_Formula([phi.id '2__'],st2);
            phi = QMITL_Parse(phi,'until', phi1, interval, phi2);
            return
        end
        
        % test expr op expr | params
        
        % parse additional params
        
        tokens = regexp(st, '(.+)\s*\|\s*(.+)','tokens');
        
        if ~isempty(tokens)
            st = tokens{1}{1};
            param_st = tokens{1}{2};
            param_tokens = regexp(param_st,'\s*,\s*','split');
            for i=1:numel(param_tokens)
                tk2 = regexp(param_tokens{i},'\s*(.+?)\s*=(.+)','tokens');
                phi.params.(tk2{1}{1}) = eval(tk2{1}{2});
            end
        end
        
        % parse operator
        
        tokens = regexp(st, '(.+)\s*<\s*(.+)','tokens');
        if ~isempty(tokens)
            phi.type='predicate';
            phi.st = st;
            
            try
                val_ = eval(tokens{1}{2});  % NM: it's not wrong, but why doing that?
                if (val_ == 0)
                    phi.params.fn =  [ '-(' tokens{1}{1} ')' ];
                else
                    phi.params.fn = [ '-(' tokens{1}{1} '- (' tokens{1}{2} '))' ];
                end
            catch
                phi.params.fn = [ '-(' tokens{1}{1} '- (' tokens{1}{2} '))' ];
            end
            
            phi.evalfn = @(mode,traj,t,params) feval('generic_predicate',mode,traj,t,params);
            return
        end
        
        tokens = regexp(st, '(.+)\s*>\s*(.+)','tokens');
        if ~isempty(tokens)
            phi.type = 'predicate';
            phi.st = st;
            
            try
                val_ = eval(tokens{1}{2});
                if (val_ == 0)
                    phi.params.fn = tokens{1}{1};
                else
                    phi.params.fn = [ '(' tokens{1}{1} ')-(' tokens{1}{2} ')' ];
                end
            catch
                phi.params.fn = [ '(' tokens{1}{1} ')-(' tokens{1}{2} ')' ];
            end
            
            phi.evalfn = @(mode,traj,t,params) feval('generic_predicate',mode,traj,t,params);
            return
        end
        
        tokens = regexp(st, '(.+)\s*=\s*(.+)', 'tokens');
        if ~isempty(tokens)
            phi.type = 'predicate';
            phi.st = st;
            if ~isfield(phi.params,'threshold')
                phi.params.threshold = 1e-14;
            end
            if ~isfield(phi.params,'max_true_value')
                phi.params.max_true_value = 1;
            end
            if ~isfield(phi.params,'alpha')
                phi.params.alpha = 1;
            end
            phi.params.fn = [ 'fun__zero(abs(' tokens{1}{1} '-(' tokens{1}{2} ')),threshold,max_true_value,alpha)'];
            phi.evalfn = @(mode,traj,t,params) feval('generic_predicate',mode,traj,t,params);
            
            return
        end
        
        % Last possibility, the formula already exists
        
        try
            id = phi.id;
            phi = struct(evalin('base', st));
            phi.id = id;
        catch
            error('QMITL_Formula:QMITL_Parse',['Unknown predicate or malformed formula: ' st]);
        end
        
    case 2
        switch(varargin{1})
            case 'not'
                phi.type = 'not';
                phi.phi = varargin{2};
                
            case 'ev'
                phi.type = 'eventually';
                phi.phi = varargin{2};                
                phi.interval = '[0 inf]';
                
            case 'evp'    % past eventually (under construction)
                phi.type = 'evp';
                phi.phi = varargin{2};
                phi.interval = '[0 inf]';
                
            case 'alw'
                phi.type = 'always' ;
                phi.phi = varargin{2};
                phi.interval = '[0 inf]';
                
            case 'eventually'
                phi.type = 'eventually' ;
                phi.phi = varargin{2};
                phi.interval = '[0 inf]';
                
            case 'always'
                phi.type = 'always' ;
                phi.phi = varargin{2};
                phi.interval = '[0 inf]';
                
            case 'andn'
                phi.type = 'andn';
                phi.phin = varargin{2}; % array of QMITL_Formula
                
            otherwise
                phi.st = varargin{1};
                phi.evalfn = @(mode,traj,t) feval(varargin{1},mode,traj,t,varargin{2});
                phi.type = 'predicate';
        end
        
    case 3
        switch(varargin{1})
            case 'or'
                phi.type = 'or' ;
                phi.phi1 = varargin{2};
                phi.phi2 = varargin{3};
                
            case 'and'
                phi.type = 'and';
                phi.phi1 = varargin{2};
                phi.phi2 = varargin{3};
                
            case '=>'
                phi.type = '=>';
                phi.phi1 = varargin{2};
                phi.phi2 = varargin{3};
                
            case 'eventually'
                phi.type = 'eventually' ;
                phi.interval = varargin{2};
                phi.phi = varargin{3};
         
            case 'evp'
                phi.type = 'evp' ;
                phi.interval = varargin{2};
                phi.phi = varargin{3};
                
            case 'always'
                phi.type = 'always' ;
                phi.interval = varargin{2};
                phi.phi = varargin{3};
                
            case 'ev'
                phi.type = 'eventually' ;
                phi.interval = varargin{2};
                phi.phi = varargin{3};
                
            case 'alw'
                phi.type = 'always' ;
                phi.interval = varargin{2};
                phi.phi = varargin{3};
        end
        
    case 4
        switch(varargin{1})
            case 'until'
                phi.type = 'until' ;
                phi.interval =  varargin{3};
                phi.phi1 = varargin{2};
                phi.phi2 = varargin{4};
        end
    otherwise
        error('QMITL_Formula:QMITL_Parse','Too many arguments.')
end

end

function [success, st1, st2] = parenthesisly_balanced_split(st, op)

% split st into st1 op st2 where st1 and st2 are parenthesisly balanced

success = 0;
st1 = '';
st2 = '';

op = ['\<' op '\>'];
[start_idx, end_idx] = regexp(st,op);

for i = 1:numel(start_idx)
    st1 = st(1:start_idx(i)-1);
    lps = regexp(st1,'(');
    num_left_par = numel(lps);
    rps = regexp(st1,')');
    num_right_par = numel(rps);
    if(num_right_par==num_left_par)
        st1 = remove_parenthesis(st1);
        st2 = remove_parenthesis(st(end_idx(i)+1:end));
        success = 1;
        return
    end
end

end

function [success, st1, st2, interval] = parenthesisly_balanced_split_interval(st, op)

% split st into st1 op st2 where st1 and st2 are parenthesisly balanced

success = 0;
st1 = '';
st2 = '';
interval = [];

op = ['\<' op '_\[(.+?)\]\>'];
[start_idx, end_idx, ~, ~, tokens] = regexp(st,op);

for i = 1:numel(start_idx)
    st1 = st(1:start_idx(i)-1);
    lps = regexp(st1,'(');
    num_left_par = numel(lps);
    rps = regexp(st1,')');
    num_right_par = numel(rps);
    if(num_right_par==num_left_par)
        st1 = remove_parenthesis(st1);
        st2 = remove_parenthesis(st(end_idx(i)+1:end));
        interval = ['[' tokens{i}{1} ']'];
        success = 1;
        return
    end
end

end

function st = remove_parenthesis(st)

tokens = regexp(st, '^\s*\(\s*(.+)\s*\)','tokens');
if ~isempty(tokens)
    st = tokens{1}{1};
end

end

