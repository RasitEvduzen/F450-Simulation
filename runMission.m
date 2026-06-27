%% F450 quadcopter simulator
%  Modules:
%     quadParams      - configuration (model, gains, visuals)
%     quadDynamics    - nonlinear model (Eq. 3.101) + RK4
%     quadControl     - multi-rate cascaded controller + mixer
%     quadTrajectory  - waypoint + lemniscate trajectory generators
%     quadVisualize   - STL animation + analysis figures
%     quatUtils       - quaternion helpers
%
%  MODE:
%     'mission'    - fly a square loop (take off, stop-turn at waypoints, land)
%     'lemniscate' - fly a Gerono lemniscate (figure-8) with the nose tangent
%
%  (system identification lives in its own script: runSysId.m)
% Written By: Rasit Evduzen
% Date: 26.06.2026
%%
clc; clear; close all;

mode = 'lemniscate';          % 'mission' (square) or 'lemniscate'

% --- load modules ---
P    = quadParams();
dyn  = quadDynamics();
ctrl = quadControl('cascade');
traj = quadTrajectory();
vis  = quadVisualize();

switch mode
    case 'mission',    runMissionMode(P, dyn, ctrl, traj, vis);
    case 'lemniscate', runLemniscateMode(P, dyn, ctrl, traj, vis);
    otherwise, error('unknown mode "%s"', mode);
end


%% ------------------------------------------------------------------------
function runMissionMode(P, dyn, ctrl, traj, vis)
    % mission waypoints [x y z yaw], NED (z<0 = up); yaw = heading after arrival
    WP = [ 0 0  0     0;        % on the ground
           0 0 -2     0;        % take off to 2 m
           4 0 -2     pi/2;     % go forward (+x), face +y
           4 4 -2     pi;       % go right (+y), face -x
           0 4 -2    -pi/2;     % go back (-x), face -y
           0 0 -2     0;        % close loop over start, turn to 0 deg
           0 0  0     0 ];      % land straight down
    segT = [5 10 10 10 10 5];
    yawT = 2.0;

    [t, REF, REFd] = traj.build(WP, segT, yawT, P.dt);
    N = numel(t);

    Binv = ctrl.mixer(P);
    st   = ctrl.initState(P);
    x = zeros(17,1); x(1:4)=P.Om_hover; x(5)=1;
    X = zeros(17,N); propAng = zeros(4,N); pa = zeros(4,1);

    LOG.vel_sp  = zeros(3,N);
    LOG.rate_sp = zeros(3,N);
    LOG.qd      = zeros(4,N);
    LOG.T_mot   = zeros(4,N);
    LOG.torque  = zeros(3,N);
    LOG.Fz      = zeros(1,N);

    n_rate = max(1, round(P.f_phys/P.f_rate));
    n_att  = max(1, round(P.f_phys/P.f_att));
    n_pos  = max(1, round(P.f_phys/P.f_pos));

    for k = 1:N
        ref.pos_sp = REF(k,1:3)';
        ref.vel_ff = REFd(k,1:3)';
        ref.yaw_sp = REF(k,4);

        run.pos  = mod(k-1, n_pos)  == 0;
        run.att  = mod(k-1, n_att)  == 0;
        run.rate = mod(k-1, n_rate) == 0;

        [delta, st] = ctrl.step(x, ref, st, P, Binv, run);
        x = dyn.rk4(x, delta, P, P.dt);

        X(:,k) = x;
        pa = pa + P.dir.*x(1:4)*P.dt;
        propAng(:,k) = pa;

        LOG.vel_sp(:,k)  = st.vel_sp;
        LOG.rate_sp(:,k) = st.rate_sp;
        LOG.qd(:,k)      = st.qd;
        LOG.T_mot(:,k)   = st.T_mot;
        LOG.torque(:,k)  = st.torque;
        LOG.Fz(k)        = st.Fz_cmd;
    end

    vis.figures(t, X, REF, P);
    vis.states(t, X, REF, LOG, P);
    vis.motors(t, LOG, P);
    vis.animate(t, X, REF, propAng, P);
end


%% ------------------------------------------------------------------------
function runLemniscateMode(P, dyn, ctrl, traj, vis)
    % Gerono lemniscate (figure-8) in x-y with a z ramp; the nose stays tangent
    % to the path (yaw = atan2(vy,vx)). The trajectory yaw-rate is fed forward
    % so the controller does not lag behind the continuously turning heading.
    opt = struct('A',4.0,'B',4.0,'z0',-2.0,'z1',-4.0,'Tloop',30.0,'laps',2);
    [t, REF, REFd, YAWRATE] = traj.lemniscate(opt, P.dt);
    N = numel(t);

    Binv = ctrl.mixer(P);
    st   = ctrl.initState(P);
    % start on the path, already pointing along it
    x = zeros(17,1); x(1:4)=P.Om_hover;
    x(5:8) = [cos(REF(1,4)/2); 0; 0; sin(REF(1,4)/2)];
    x(15:17) = REF(1,1:3)';
    X = zeros(17,N); propAng = zeros(4,N); pa = zeros(4,1);

    LOG.vel_sp  = zeros(3,N);  LOG.rate_sp = zeros(3,N);
    LOG.qd      = zeros(4,N);  LOG.T_mot   = zeros(4,N);
    LOG.torque  = zeros(3,N);  LOG.Fz      = zeros(1,N);

    n_rate = max(1, round(P.f_phys/P.f_rate));
    n_att  = max(1, round(P.f_phys/P.f_att));
    n_pos  = max(1, round(P.f_phys/P.f_pos));

    for k = 1:N
        ref.pos_sp = REF(k,1:3)';
        ref.vel_ff = REFd(k,1:3)';
        ref.yaw_sp = REF(k,4);
        ref.yaw_rate_ff = YAWRATE(k);

        run.pos  = mod(k-1, n_pos)  == 0;
        run.att  = mod(k-1, n_att)  == 0;
        run.rate = mod(k-1, n_rate) == 0;

        [delta, st] = ctrl.step(x, ref, st, P, Binv, run);
        x = dyn.rk4(x, delta, P, P.dt);

        X(:,k) = x;
        pa = pa + P.dir.*x(1:4)*P.dt;  propAng(:,k) = pa;
        LOG.vel_sp(:,k)  = st.vel_sp;   LOG.rate_sp(:,k) = st.rate_sp;
        LOG.qd(:,k)      = st.qd;       LOG.T_mot(:,k)   = st.T_mot;
        LOG.torque(:,k)  = st.torque;   LOG.Fz(k)        = st.Fz_cmd;
    end

    vis.figures(t, X, REF, P);
    vis.states(t, X, REF, LOG, P);
    vis.motors(t, LOG, P);
    vis.animate(t, X, REF, propAng, P);
end