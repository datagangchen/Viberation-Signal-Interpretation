function [val, tau] = QMITL_Eval(Sys, phi, P, trajs, varargin)
%QMITL_EVAL computes the satisfaction function of a property for one or
% many traces.
% 
% Synopsis: [val, tau] = QMITL_Eval(Sys, phi, P, trajs[, taus][, method])
% 
% Inputs:
%  - Sys    : the system
%  - phi    : a QMITL Formula
%  - P      : is a parameter set which contains one parameter vector only
%             used for properties parameters.
%  - trajs  : is a structure with fields X and time. It may contains many
%             trajectories. In this case, all will be checked wrt the
%             property parameter described in P.
%  - taus   : (Optional, default=traj.time for each traj in trajs) is the
%             time, possibly an array, when to eval the satisfaction of the
%             property. All time points not belonging to traj.time will be
%             linearly interpolated.
%  - method : (Optional, default='thom') is a string indicating the method
%             to use to evaluate the formula. It is either 'classic' or
%             'thom'. The main difference between the two methods is that
%             the classic method uses a fixed time step while 'thom' method
%             uses a variable time step.
% 
% Outputs:
%  - val : is a cell array of dimension 1 x numel(trajs). Each cell
%          contains an line array describing the evaluation of phi at each
%          time of tau. If numel(trajs) is 1, val is the content of its
%          only cell to avoid a useless cell array (so, it is a line
%          array).
%  - tau : is a cell array of dimension 1 x numel(trajs). It contains time
%          points at which the formula is evaluated or interpolated. If the
%          parameter taus is provided, tau is equal to taus. If
%          numel(trajs) is 1, tau contains the content of the only cell
%          (this avoid a useless cell array), and thus becomes a line
%          array.
% 
% Example (Lorentz84):
%  CreateSystem;
%  P = CreateParamSet(Sys, 'F', [10, 20]);
%  P = SetParam(P, 'theta', 2);
%  P = ComputeTraj(Sys, P, 0:.01:10);
%  phi = QMITL_Formula('phi','ev_[0,1] (x0[t]>theta)');
%  [val,tau] = QMITL_Eval(Sys, phi, P, P.traj(P.traj_ref(1)), 0)
%  [val,tau] = QMITL_Eval(Sys, phi, P, P.traj(P.traj_ref(1)), [3,7])
% 
%See also SEvalProp QMITL_Formula
%

switch nargin
    case 4
        [val, tau] = QMITL_EvalThom(Sys, phi, P, trajs);
    case 5
        if ~ischar(varargin{1})
            [val, tau] = QMITL_EvalThom(Sys, phi, P, trajs, varargin{1});
        elseif strcmpi(varargin{1},'classic')
            [val, tau] = QMITL_EvalClassic(Sys, phi, P, trajs);
        elseif strcmpi(varargin{1},'online')
            [val, tau] = QMITL_EvalClassicOnline(Sys, phi, P, trajs);
        else
            [val, tau] = QMITL_EvalThom(Sys, phi, P, trajs);
        end
    case 6
        if strcmpi(varargin{2},'classic')
            [val, tau] = QMITL_EvalClassic(Sys, phi, P, trajs, varargin{1});
        elseif strcmpi(varargin{1},'online')
            [val, tau] = QMITL_EvalClassicOnline(Sys, phi, P, trajs);
        else
            [val, tau] = QMITL_EvalThom(Sys, phi, P, trajs, varargin{1});
        end
end

end
