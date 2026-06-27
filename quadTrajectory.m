function out = quadTrajectory()
%QUADTRAJECTORY  Reference trajectory generators.
%   out.build(WP, segT, yawT, dt)        waypoint mission: fly straight, then
%                                         turn in place (LSPB + cosine yaw)
%   out.lemniscate(opt, dt)               Gerono lemniscate (figure-8) with the
%                                         nose kept tangent to the path
%
%   build returns [t, REF, REFd]; lemniscate returns [t, REF, REFd, YAWRATE]
%   where REF = Nx4 [x y z yaw], REFd = Nx4 velocity, YAWRATE = Nx1 yaw-rate
%   feedforward (rad/s). For lemniscate the heading tracks the velocity vector
%   (yaw = atan2(vy, vx)), like a coordinated turn.

    out.build      = @buildTraj;
    out.lemniscate = @buildLemniscate;
end

function [T, REF, REFd, YAWRATE] = buildLemniscate(opt, dt)
%BUILDLEMNISCATE  Lemniscate of Gerono in the x-y plane with a z ramp.
%   x = A sin(w t),  y = (B/2) sin(2 w t),  z = z0 + (z1-z0) t/Ttot
%   Heading is tangent to the path: yaw = atan2(ydot, xdot).
    if nargin < 1 || isempty(opt), opt = struct; end
    A    = getf(opt,'A',4.0);        % x amplitude [m]
    B    = getf(opt,'B',4.0);        % y span (lobe height) [m]
    z0   = getf(opt,'z0',-2.0);      % start altitude [m]
    z1   = getf(opt,'z1',-4.0);      % end altitude [m]
    Tloop= getf(opt,'Tloop',30.0);   % time for one full lemniscate [s]
    laps = getf(opt,'laps',1);       % number of loops

    Ttot = Tloop*laps;  w = 2*pi/Tloop;
    T = (0:dt:Ttot)';
    x  = A*sin(w*T);          xd = A*w*cos(w*T);
    y  = (B/2)*sin(2*w*T);    yd = (B/2)*2*w*cos(2*w*T);
    z  = z0 + (z1-z0)*(T/Ttot);  zd = (z1-z0)/Ttot*ones(numel(T),1);
    yaw = atan2(yd, xd);

    REF  = [x, y, z, yaw];
    REFd = [xd, yd, zd, zeros(numel(T),1)];
    YAWRATE = gradient(unwrap(yaw), dt);   % yaw-rate feedforward
end

function [T,P,Pd] = buildTraj(WP, segT, yawT, TS)
    T=[]; P=[]; Pd=[]; toff=0;
    for i = 1:size(WP,1)-1
        qi = WP(i,:); qf = WP(i+1,:);

        % --- 1) translation segment (heading held at qi yaw) ---
        t = (0:TS:segT(i))'; n = numel(t);
        p = zeros(n,4); pd = zeros(n,4);
        for k = 1:3
            if abs(qf(k)-qi(k)) < 1e-9, p(:,k)=qi(k); continue; end
            V = 1.5*abs(qf(k)-qi(k))/segT(i);
            [pp,pv] = SCurveTrajectory(qi(k), qf(k), segT(i), V, TS);
            p(:,k)=pp; pd(:,k)=pv;
        end
        p(:,4) = qi(4);                      % yaw held while translating
        if i > 1, t(1)=[]; p(1,:)=[]; pd(1,:)=[]; end
        T=[T; t+toff]; P=[P; p]; Pd=[Pd; pd]; toff = toff + segT(i);

        % --- 2) in-place yaw turn on arrival (position held at qf) ---
        dpsi = qf(4)-qi(4); dpsi = mod(dpsi+pi,2*pi)-pi;
        if abs(dpsi) > 1e-6
            tw = (TS:TS:yawT)'; nw = numel(tw);
            pw = zeros(nw,4); vw = zeros(nw,4);
            pw(:,1)=qf(1); pw(:,2)=qf(2); pw(:,3)=qf(3);   % hold position
            for j = 1:nw
                f = tw(j)/yawT;
                pw(j,4) = qi(4) + dpsi*(0.5 - 0.5*cos(pi*f));
            end
            T=[T; tw+toff]; P=[P; pw]; Pd=[Pd; vw]; toff = toff + yawT;
        end
    end
end

function [p,pd] = SCurveTrajectory(qi, qf, tf, V, TS)
    t  = (0:TS:tf)';
    V  = abs(V)*sign(qf-qi);
    tb = (qi - qf + V*tf)/V;
    a  = V/tb;
    p  = zeros(length(t),1); pd = p;
    for i = 1:length(t)
        if t(i) <= tb
            p(i)=qi + a/2*t(i)^2;                          pd(i)=a*t(i);
        elseif t(i) <= (tf-tb)
            p(i)=(qf+qi-V*tf)/2 + V*t(i);                  pd(i)=V;
        else
            p(i)=qf - a/2*tf^2 + a*tf*t(i) - a/2*t(i)^2;   pd(i)=a*tf - a*t(i);
        end
    end
end

function v = getf(s,f,d), if isfield(s,f), v=s.(f); else, v=d; end, end