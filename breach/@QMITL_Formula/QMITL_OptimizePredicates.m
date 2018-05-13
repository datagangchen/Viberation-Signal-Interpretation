function phi = QMITL_OptimizePredicates(S, phi)
%QMITL_OPTIMIZEPREDICATE pre-parses predicates in formula to make its
% evaluation faster. The function extracts all predicates present in phi
% and optimizes them.
% 
% Synopsis: phi = QMITL_OptimizePredicates(S, phi)
% 
% Inputs:
%   - S   : can be a system (Sys) or a parameter set (usually P)
%   - phi : the formula to optimize
% 
% Output:
%  - phi : the optimized version of the formula
% 
% Example (Lorentz84):
%  CreateSystem;
%  P = CreateParamSet(Sys);
%  P = ComputeTraj(Sys,P,0:0.1:10);
%  phi = QMITL_Formula('phi','x0[t]>2');
%  tic ; for ii=1:2000 ; QMITL_Eval(Sys,phi,P,P.traj,0:0.1:10) ; end ; toc  % takes around 5.5s on Core 2 Duo E8400
%  phi = QMITL_OptimizePredicates(Sys,phi);
%  tic ; for ii=1:2000 ; QMITL_Eval(Sys,phi,P,P.traj,0:0.1:10) ; end ; toc  % takes around 4.4s on Core 2 Duo E8400
% 
%See also QMITL_Formula SEvalProp QMITL_OptimizeFormula
%

if strcmp(phi.type, 'predicate')
    phi = parse_predicate(S, phi);
    return;
end
if ~isempty(phi.phi)
    phi.phi = QMITL_OptimizePredicates(S,phi.phi);
end
if ~isempty(phi.phi1)
    phi.phi1 = QMITL_OptimizePredicates(S,phi.phi1);
end
if ~isempty(phi.phi2)
    phi.phi2 = QMITL_OptimizePredicates(S,phi.phi2);
end
if ~isempty(phi.phin)
    for ii=1:numel(phi.phin)
        phi.phin(ii) = QMITL_OptimizePredicates(S,phi.phin(ii));
    end
end

end
