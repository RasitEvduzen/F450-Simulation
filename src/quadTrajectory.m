function out = quadTrajectory()
%QUADTRAJECTORY  Reference trajectory generators.
%   out.build(WP, segT, yawT, dt, useDouble)  waypoint mission: fly straight,
%                                         then turn in place. useDouble (opt,
%                                         default false) picks the per-axis
%                                         translation profile: false = 3-segment
%                                         trapezoidal, true = 7-segment S-curve
%                                         (finite jerk, smoother accel).
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

function [T,P,Pd] = buildTraj(WP, segT, yawT, TS, useDouble)
    if nargin < 5 || isempty(useDouble), useDouble = false; end
    T=[]; P=[]; Pd=[]; toff=0;
    for i = 1:size(WP,1)-1
        qi = WP(i,:); qf = WP(i+1,:);

        % --- 1) translation segment (heading held at qi yaw) ---
        t = (0:TS:segT(i))'; n = numel(t);
        p = zeros(n,4); pd = zeros(n,4);
        for k = 1:3
            if abs(qf(k)-qi(k)) < 1e-9, p(:,k)=qi(k); continue; end
            V = 1.5*abs(qf(k)-qi(k))/segT(i);     % mid-range speed (Ta = tf/3)
            if useDouble
                Tj = segT(i)/8;                   % safe: Tj < Ta/2 = tf/6
                [pp,pv] = DoubleSCurve(qi(k), qf(k), segT(i), V, Tj, TS);
            else
                [pp,pv] = SCurve(qi(k), qf(k), segT(i), V, TS);
            end
            p(1:n,k)=pp(1:n); pd(1:n,k)=pv(1:n);
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

function [p,pd] = SCurve(qi, qf, tf, V, TS)
% 3-segment trapezoidal-velocity profile (position p [m], velocity pd [m/s]).
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

function [p,pd] = DoubleSCurve(qi, qf, tf, V, Tj, TS)
% 7-segment constant-jerk S-curve (position p [m], velocity pd [m/s]).
% Position cubic, velocity quadratic, acceleration piecewise-linear (C2),
% jerk piecewise-constant. Feasibility: V in [h/tf,2h/tf], Tj in (0,Ta/2].
    h = abs(qf-qi); sgn = sign(qf-qi);
    t = (0:TS:tf)'; N = numel(t);
    if h == 0, p = repmat(qi,N,1); pd = zeros(N,1); return; end

    Ta   = tf - h/V;             % accel-phase duration [s]
    amax = V/(Ta - Tj);          % peak acceleration     [m/s^2]
    J    = amax/Tj;              % jerk magnitude         [m/s^3]
    Tv   = max(tf - 2*Ta, 0);    % cruise duration        [s]

    t1 = Tj; t2 = Ta-Tj; t3 = Ta; t4 = Ta+Tv; t5 = Ta+Tv+Tj; t6 = tf-Tj;
    v1 = J*Tj^2/2;          p1 = J*Tj^3/6;
    d2 = t2-t1;             v2 = v1+amax*d2;  p2 = p1+v1*d2+amax*d2^2/2;
                                              p3 = p2+v2*Tj+amax*Tj^2/2-J*Tj^3/6;
    p4 = p3+V*Tv;
    v5 = V-J*Tj^2/2;        p5 = p4+V*Tj-J*Tj^3/6;
    d6 = Ta-2*Tj;           v6 = v5-amax*d6;  p6 = p5+v5*d6-amax*d6^2/2;

    p = zeros(N,1); pd = zeros(N,1);
    m1 = t<=t1; m2 = t>t1 & t<=t2; m3 = t>t2 & t<=t3; m4 = t>t3 & t<=t4;
    m5 = t>t4 & t<=t5; m6 = t>t5 & t<=t6; m7 = t>t6;
    a=t(m1);     pd(m1)=J.*a.^2/2;            p(m1)=J.*a.^3/6;
    a=t(m2)-t1;  pd(m2)=v1+amax*a;            p(m2)=p1+v1*a+amax*a.^2/2;
    a=t(m3)-t2;  pd(m3)=v2+amax*a-J*a.^2/2;   p(m3)=p2+v2*a+amax*a.^2/2-J*a.^3/6;
    a=t(m4)-t3;  pd(m4)=V;                    p(m4)=p3+V*a;
    a=t(m5)-t4;  pd(m5)=V-J*a.^2/2;           p(m5)=p4+V*a-J*a.^3/6;
    a=t(m6)-t5;  pd(m6)=v5-amax*a;            p(m6)=p5+v5*a-amax*a.^2/2;
    a=t(m7)-t6;  pd(m7)=v6-amax*a+J*a.^2/2;   p(m7)=p6+v6*a-amax*a.^2/2+J*a.^3/6;

    p  = qi + sgn*p;
    pd = sgn*pd;
end

function v = getf(s,f,d), if isfield(s,f), v=s.(f); else, v=d; end, end