function varargout = sqrtm(varargin)
%SQRTM (overloaded)

% Author Johan Löfberg
% $Id: sqrtm_internal.m,v 1.1 2006/11/15 10:03:57 joloef Exp $

switch class(varargin{1})

    case 'double' % What is the numerical value of this argument (needed for displays etc)
        % SHOULD NEVER HAPPEN, THIS SHOULD BE CAUGHT BY BUILT-IN
        varargout{1} = sqrt(varargin{1});
%        error('Overloaded SDPVAR/SQRTM CALLED WITH DOUBLE. Report error')

    case 'sdpvar' % Overloaded operator for SDPVAR objects.
        if length(varargin{1}) == 1
            varargout{1} = yalmip('addEvalVariable',mfilename,varargin{1});
        else
            error('SQRTM is only implemented for scalar argument')
%             y = [];
%             for i = 1:length(varargin{1})
%                 y = [y;yalmip('addEvalVariable',mfilename,extsubsref(varargin{1},i))];
%             end
%             varargout{1} = y;
        end
        
    case 'char' % YALMIP sends 'model' when it wants the epigraph or hypograph
        switch varargin{1}
            case {'graph','milp'}
                t = varargin{2};
                X = varargin{3};
                
                % This is different from so called extended operators
                % Just do it!
                F = SetupEvaluationVariable(varargin{:});
                
                % Now add your own code, such as domain constraints
                F = F + set(X > eps);
                
                % Let YALMIP know about convexity etc                
                varargout{1} = F;
                varargout{2} = struct('convexity','none','monotonicity','none','definiteness','none');
                varargout{3} = X;                
                
            otherwise
                error('SDPVAR/SQRTM called with CHAR argument?');
        end
    otherwise
        error('SDPVAR/SQRTM called with CHAR argument?');
end
