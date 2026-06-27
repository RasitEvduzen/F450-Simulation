function out = quadControl(name)
%QUADCONTROL  Controller factory. Selects a controller by name and returns it.
%   C = quadControl()            -> default 'cascade' controller
%   C = quadControl('cascade')   -> multi-rate cascaded P/PID controller
%   C = quadControl('indi')      -> incremental nonlinear dynamic inversion
%
%   Every controller lives in its own file (quadControlCascade.m,
%   quadControlINDI.m, ...) and exposes the same interface:
%     C.name, C.mixer(P), C.initState(P), C.step(x, ref, st, P, Binv, run)
%   so the top-level scripts can switch controller by changing only this name.

    if nargin < 1 || isempty(name), name = 'cascade'; end
    switch lower(name)
        case 'cascade', out = quadControlCascade();
        case 'indi',    out = quadControlINDI();
        otherwise, error('quadControl: unknown controller "%s"', name);
    end
end