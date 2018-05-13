function [props_names, props] = QMITL_ReadFile(fname)
%QMITL_READFILE defines formulas from a text file.
% 
% Synopsis: [props_names, prop] = QMITL_ReadFile(fname)
% 
% Input:
%  - fname the text file containing formulas. This text file fname should
%          consist of a sequence of definitions of the form:
%
%          phi1 := some formula
%          phi2 := some formula (might depend on phi1)
%          etc
%
%          blanks and comments beginning with a # are allowed eg :
%
%          phis.stl:
%
%          mu := x0[t] > 3  # what a great predicate !
%          phi1 := mu until mu # I could have found something else...
%          phi2 := alw_[0,2.3] phi1    # what about that !
%
%          -- end of phis.stl
% 
% Outputs:
%  - props_names : is a cell array describing the names of the created 
%                  formulas.
%  - props       : is a cell array containing the defined STL formulas in
%                  the same order than props_name.
% 
% Example (Lorentz84):
%  [props_names,props] = QMITL_ReadFile('oscil_prop.stl');
%  props_names
% 
%See also QMITL_Formula RecoverFormula
%

fid = fopen(fname,'r');

if(fid==-1)
    error('QMITL_ReadFile:OpeningError',['Couldn''t open file ' fname]);
end

tline = fgetl(fid);

current_id = '';
current_formula = '';
num_line = 0;
props_names = {};
props = {};

while ischar(tline)
    num_line = num_line+1;
    
    % first, dismiss comments  (anything starting with a #) and starting spaces
    tline = regexprep(tline, '^\s*','');
    tline = regexprep(tline, '\s*$','');
    
    if regexp(tline, '^\#')
        tline = '';
    end
    
    tokens = regexp(tline, '\s*(.*?)\#(.+)|\s*(.*)','tokens');
    if ~isempty(tokens)
        tline = tokens{1}{1};
    else
        tline = '';
    end
    
    if ~isempty(tline)
        % second, checks if we're starting the def. of a new formula
        tokens = regexp(tline, '(\w+)\s*:=(.*)','tokens');
        if ~isempty(tokens)
            % ok try wrapping up what we have so far before starting a new formula
            if ~isempty(current_id)
                try
                    %fprintf(['assign ' current_id ':=' current_formula '\n']);
                    props = [props, {QMITL_Formula(current_id, current_formula)}]; %#ok<*AGROW>
                    %assignin('caller', current_id, current_formula);
                    props_names = [props_names, {current_id}]; %#ok<AGROW>
                catch err
                    fprintf(['ERROR: Problem with formula ' current_id ' before line ' ...
                        int2str(num_line-1) '\n']);
                    rethrow(err);
                end
            end % we're done : if current_id is empty this is our first formula
            
            % test the new id
            current_id = tokens{1}{1};
            try
                assignin('base', current_id, 0);
            catch %#ok<CTCH>
                error('QMITL_ReadFile:IdError',[current_id ' on line ' int2str(num_line) ' is not a valid id.']);
            end
            
            % start definition of formula
            current_formula = tokens{1}{2};
            
        else % we're continuing the definition of a formula
            
            if isempty(current_id)
                error('QMITL_ReadFile:MissingId',['On line ' int2str(num_line) ', no id given yet.']);
            end
            current_formula = [current_formula ' ' tline]; %#ok<AGROW>
        end
    end
    
    tline = fgetl(fid);
end

% last formula
try
    %fprintf(['assign ' current_id ':=' current_formula '\n']);
    props = [props, {QMITL_Formula(current_id, current_formula)}];
    %assignin('caller', current_id, phi);
    props_names = [props_names, {current_id}];
catch err
    error('QMITL_ReadFile:Error',['Problem with formula ' current_id ' after line '  ...
        int2str(num_line-1) ': ' err.message ]);
end

end
