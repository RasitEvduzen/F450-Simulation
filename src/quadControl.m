function out = quadControl(name)
%QUADCONTROL  Pluggable, multi-rate cascaded controller for the quadcopter.
%   C = quadControl()            -> default 'cascade' controller
%   C = quadControl('cascade')   -> cascaded P/PID controller (this base model)
%
%   Common interface (controller-agnostic):
%     C.name                            controller identifier
%     C.mixer(P)                        -> Binv (control allocation, shared)
%     C.initState(P)                    -> controller state (integrators + held outputs)
%     C.step(x, ref, st, P, Binv, run)  -> [delta, st]
%
%   Multi-rate: the cascade runs three loops (position/velocity, attitude,
%   rate), each slower than the physics step. `run` flags {pos, att, rate}
%   say which loops fire on a given physics tick; loops that don't fire reuse
%   their held output (zero-order hold). Each loop uses its own period
%   (st.dt_pos / dt_att / dt_rate) for integration and differentiation.
%
%   To add a controller (dynamic inversion, LQR, MPC): implement a step with
%   the same signature and register it in the switch below. The mixer and
%   rotor-command mapping are shared.

    if nargin < 1 || isempty(name), name = 'cascade'; end
    out.name      = name;
    out.mixer     = @buildMixer;
    out.rotorCmd  = @rotorCommand;
    switch lower(name)
        case 'cascade'
            out.initState = @cascadeInit;
            out.step      = @cascadeStep;
        % case 'dyninv'   % future: dynamic inversion
        % case 'lqr'      % future: LQR
        % case 'mpc'      % future: MPC
        otherwise
            error('quadControl: unknown controller "%s"', name);
    end
end

%% control allocation (mixer) ----
function Binv = buildMixer(P)
    kq = P.cQ/P.cT; B = zeros(4,4);
    for i = 1:4
        tau = cross(P.pos(i,:)', P.n);
        B(1,i) = P.n(3);
        B(2,i) = tau(1);
        B(3,i) = tau(2);
        B(4,i) = -P.dir(i)*P.n(3)*kq;
    end
    Binv = pinv(B);
end

%% [collective; torque] -> per-rotor PWM ----
function [delta, T_mot] = rotorCommand(coll, torque, P, Binv)
    Fz_cmd = coll/P.hover_thr*P.m*P.g;
    T_mot  = max(Binv*[Fz_cmd; torque], 0.01);
    Om_cmd = sqrt(T_mot/P.cT);
    delta  = max(min((Om_cmd + P.cQ*Om_cmd.^2/P.kT)/P.Ehat, 1), 0);
end

%% Cascade controller (multi-rate): position -> velocity -> attitude -> rate
function st = cascadeInit(P)
    st.rate_int = zeros(3,1);
    st.vel_int  = zeros(3,1);
    st.w_prev   = zeros(3,1);
    st.v_prev   = zeros(3,1);
    % held outputs (zero-order hold between loop updates)
    st.coll     = -0.5;                 % held collective (hover)
    st.qd       = [1;0;0;0];            % held desired attitude
    st.rate_sp  = zeros(3,1);           % held rate setpoint
    st.delta    = P.delta_hover*ones(4,1);  % held motor command (hover)
    st.vel_sp   = zeros(3,1);           % logged velocity setpoint
    st.T_mot    = P.m*P.g/4*ones(4,1);  % logged per-motor thrust (hover)
    st.torque   = zeros(3,1);           % logged body torque command
    st.Fz_cmd   = -P.m*P.g;             % logged collective thrust (hover, N)
    % per-loop periods
    st.dt_pos = 1/P.f_pos;
    st.dt_att = 1/P.f_att;
    st.dt_rate= 1/P.f_rate;
end

function [delta, st] = cascadeStep(x, ref, st, P, Binv, run)
    Q = quatUtils();
    q = x(5:8); w = x(9:11); pos = x(15:17); vel = x(12:14);

    % --- POSITION (P) -> velocity sp ; VELOCITY (PID) -> accel sp ; bridge ---
    if run.pos
        dt = st.dt_pos;
        vel_sp = (ref.pos_sp - pos).*P.Kp_pos + ref.vel_ff;
        vel_sp(1:2) = max(min(vel_sp(1:2), P.vmax_xy), -P.vmax_xy);
        vel_sp(3)   = max(min(vel_sp(3),   P.vmax_dn), -P.vmax_up);
        st.vel_sp = vel_sp;            % logged for analysis

        vel_err = vel_sp - vel;
        vel_dot = (vel - st.v_prev)/dt;  st.v_prev = vel;
        acc_sp  = vel_err.*P.Kp_vel + st.vel_int - vel_dot.*P.Kd_vel;
        st.vel_int = st.vel_int + vel_err.*P.Ki_vel*dt;
        st.vel_int(3) = max(min(st.vel_int(3), P.g), -P.g);

        z_sf   = -P.g + acc_sp(3);
        body_z = [-acc_sp(1); -acc_sp(2); -z_sf]; body_z = body_z/norm(body_z);
        thr_z  = acc_sp(3)*(P.hover_thr/P.g) - P.hover_thr;
        st.coll = min(thr_z/body_z(3), -0.12);
        qz = Q.between([0;0;-1], -body_z/norm(body_z));
        st.qd = Q.mul(qz, Q.yaw(ref.yaw_sp));  st.qd = st.qd/norm(st.qd);
    end

    % --- ATTITUDE (quaternion P) -> rate sp ; + yaw-rate feedforward ---
    if run.att
        st.rate_sp = attitudeControl(q, st.qd, P, Q);
        if isfield(ref,'yaw_rate_ff')
            st.rate_sp(3) = st.rate_sp(3) + ref.yaw_rate_ff;  % nose-tangent FF
        end
    end

    % --- RATE (PID + FF) -> torque ; allocation ---
    % runs at its own rate; between updates the previous delta is reused
    if run.rate
        dt = st.dt_rate;
        rate_err  = st.rate_sp - w;
        ang_accel = (w - st.w_prev)/dt;  st.w_prev = w;
        torque = P.Kp_rate.*rate_err + st.rate_int - P.Kd_rate.*ang_accel + P.Kff_rate.*st.rate_sp;
        i_factor = max(0, 1 - (rate_err/deg2rad(400)).^2);
        st.rate_int = st.rate_int + i_factor.*P.Ki_rate.*rate_err*dt;
        st.rate_int = max(min(st.rate_int, P.int_lim), -P.int_lim);
        st.torque = torque;                          % logged: [tau_x; tau_y; tau_z]
        st.Fz_cmd = st.coll/P.hover_thr*P.m*P.g;      % logged: collective thrust [N]
        [st.delta, st.T_mot] = rotorCommand(st.coll, torque, P, Binv);
    end
    delta = st.delta;
end

function rate_sp = attitudeControl(q, qd, P, Q)
    q = q/norm(q); qd = qd/norm(qd);
    e_z = Q.dcm_z(q); e_z_d = Q.dcm_z(qd);
    qd_red = Q.between(e_z, e_z_d);
    if abs(qd_red(2)) > 1-1e-5 || abs(qd_red(3)) > 1-1e-5
        qd_full = qd;
    else
        qd_full = Q.mul(qd_red, q);
    end
    q_mix = Q.canon(Q.mul(Q.inv(qd_full), qd));
    q_mix(1) = max(min(q_mix(1),1),-1);
    q_mix(4) = max(min(q_mix(4),1),-1);
    qd_m = Q.mul(qd_full, [cos(P.yaw_w*acos(q_mix(1))); 0; 0; sin(P.yaw_w*asin(q_mix(4)))]);
    qe = Q.canon(Q.mul(Q.inv(q), qd_m));
    rate_sp = max(min(2*qe(2:4).*P.Kp_att, P.rate_max), -P.rate_max);
end