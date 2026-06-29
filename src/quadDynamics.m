function out = quadDynamics()
%QUADDYNAMICS  Nonlinear flight-mechanical model + RK4 integrator.
%   Returns a struct of handles:
%     out.deriv(x,delta,P)   -> state derivative (Stephan Eq. 3.101 + battery)
%     out.rk4(x,delta,P,dt)  -> one RK4 step (with quaternion renormalise)
%     out.vbat(x,P)          -> battery terminal voltage E = Ehat(Q) - Rbat*Isig
%
%   STATE (18): x = [Omega(4); q(4); omega_body(3); v_g(3); r_g(3); Qbat(1)]
%     Omega      - rotor speeds [rad/s]                 (propulsion state)
%     q          - attitude quaternion (body->geodetic), scalar-first
%     omega_body - body angular rate (p,q,r) [rad/s]
%     v_g, r_g   - velocity / position in geodetic NED  [m/s, m]
%     Qbat       - battery state of charge [Ah]          (energy state, NEW)
%   Input (4): delta in [0,1] - per-rotor PWM command (duty cycle).
%
%   ----------------------------------------------------------------------
%   EQUATIONS OF MOTION (Stephan "Multicopter Flight Control")
%   ----------------------------------------------------------------------
%   The plant couples three dynamic subsystems, each with its own time scale:
%   the rigid body (slow, ~s), the rotors (fast, ~ms), and the battery
%   (very slow, ~min). They are integrated together as one 18-state ODE.
%
%   (A) MOTOR / ESC  (Eqs. 3.31-3.36).  Per rotor i:
%         Emot_i = delta_i * E           ESC modulates battery voltage (3.32a)
%         Imot_i = (Emot_i - ke*Omega_i)/Rmot   armature current, quasi-static
%                                               (dI/dt=0, L neglected, Eq. 3.33)
%         MR_i   = kT * Imot_i           motor torque (3.34)
%         JR*Omdot_i = MR_i - MF - Q_i   rotor spin-up (3.31/3.36)
%       Combined (substituting Imot, with ke = kT):
%         JR*Omdot_i = (kT/Rmot)*(E*delta_i - kT*Omega_i) - MF - cQ*Omega_i^2
%       E here is the live battery terminal voltage (see C), so a sagging /
%       discharging pack directly weakens every rotor.
%
%   (B) BATTERY  (LiPo discharge, Eqs. 3.52-3.54).
%         Imot_i  = (E*delta_i - kT*Omega_i)/Rmot       per-motor current
%         Isig    = sum_i delta_i*Imot_i                pack current drawn (3.42b)
%         Ehat(Q) = E0 - Epol*(Q0/Q) + Eexp*exp(-(Q0-Q)/Qexp)   open-circuit (3.53)
%         E       = Ehat(Q) - Rbat*Isig                 terminal voltage w/ sag
%         Qdot    = -Isig/3600                          Coulomb counting (3.54)
%                                                       (/3600: A -> Ah per second)
%       So Q falls as charge is drawn; Ehat(Q) drops along the LiPo curve;
%       the sag term Rbat*Isig adds an instantaneous load-dependent dip.
%
%   (C) RIGID BODY  (Eq. 3.101).
%         per-rotor thrust/drag:  T_i = cT*Omega_i^2,  Q_i = cQ*Omega_i^2
%         F_Sigma   = sum_i n_i*T_i
%         tau_Sigma = sum_i p_i x (n_i*T_i) - th_i*n_i*Q_i
%         Kinematics:   qdot   = 1/2 * Theta(q) * omega
%         Rotational:   Jb*omdot = tau_Sigma
%                                  - omega x ( Jb*omega + sum_i th_i*JR*n_i*Omega_i )
%         Translation:  m*vdot  = R(q)*F_Sigma + m*g
%         Position:     rdot    = v
%       n_i is the rotor thrust axis, p_i its position, th_i = +/-1 its spin
%       direction, g = [0;0;+g] (NED, gravity points +z / down).
%
%   ----------------------------------------------------------------------
%   MODELLING ASSUMPTIONS (what we keep, what we simplify, and why)
%   ----------------------------------------------------------------------
%     1. Rigid body, fixed mass and inertia Jb (no payload change).
%     2. Near-hover / low-speed aerodynamics: undisturbed inflow V_inf = 0,
%        so in-plane force H_i = 0 and T_i, Q_i are quadratic polars in
%        Omega_i (Eq. 3.104). High-speed airframe aero is out of scope here.
%     3. Out of ground effect (cGE = 1).
%     4. No airframe aerodynamic drag (C = 0), no wind, no disturbances.
%     5. Rotors horizontal, thrust axis vertical in body: n_i = [0;0;-1];
%        rotor spin axis +z in the gyroscopic term (consistent with n_i).
%     6. Gyroscopic momentum coupling (omega x sum th_i*JR*n_i*Omega_i) is KEPT;
%        the spin-up reaction torque (sum th_i*JR*n_i*Omdot_i) is OMITTED
%        (JR small, numerically stiff at 1 kHz; see the note at the wdot line).
%     7. MOTOR (NEW vs. baseline): the lumped tuned form has been replaced by
%        the physical Eq. 3.36 with real kT, ke, Rmot. Armature inductance L
%        is neglected (electrical time constant L/Rmot ~0.3 ms is ~1000x faster
%        than the rotor, so current is quasi-static, Eq. 3.33). ESC is ideal
%        (no extra lag); the rotor's own first-order lag JR*Rmot/kT^2 already
%        captures the dominant actuator delay.
%     8. BATTERY (NEW vs. baseline): the constant-voltage assumption is dropped.
%        The pack now discharges (Q is a state); terminal voltage follows the
%        LiPo curve Ehat(Q) with an internal-resistance sag Rbat*Isig. This
%        makes thrust fade as the battery drains (endurance is now emergent).

    out.deriv = @quadDeriv;
    out.rk4   = @rk4;
    out.vbat  = @batteryVoltage;
    out.current = @packCurrent;
end

function xdot = quadDeriv(x, delta, P)
    Q  = quatUtils();
    Om = x(1:4); q = x(5:8); w = x(9:11); v = x(12:14); Qbat = x(18);
    q  = q/norm(q);
    delta = min(max(delta,0),1);              % PWM physically in [0,1]

    % --- Battery: open-circuit voltage on the LiPo curve, then the load --
    % current it must supply, then the terminal voltage including I*R sag.
    Ehat = ocvFromCharge(Qbat, P);            % open-circuit voltage Ehat(Q)
    % per-motor armature current (quasi-static, Eq. 3.33); clamp >=0 since the
    % ESC cannot sink current (no active braking)
    Imot = max((Ehat*delta - P.ke*Om)/P.Rmot, 0);
    Isig = sum(delta.*Imot);                  % pack current drawn (Eq. 3.42b)
    Isig = min(Isig, P.bat.Imax);             % ESC/pack current limit (C-rating)
    Ebat = Ehat - P.bat.Rbat*Isig;            % terminal voltage with sag (3.53)

    % --- Motor dynamics (Eq. 3.36), driven by the LIVE terminal voltage ----
    % JR*Omdot = (kT/Rmot)*(Ebat*delta - kT*Omega) - MF - cQ*Omega^2
    Qd    = P.cQ*Om.^2;                        % rotor drag torque Q_i
    MR    = (P.kT/P.Rmot)*(Ebat*delta - P.ke*Om);   % electromagnetic torque
    Omdot = (MR - P.MF - Qd)/P.JR;

    % --- per-rotor thrust, body force/torque, rotor angular momentum -------
    T = P.cT*Om.^2;
    Fsum = [0;0;0]; tau = [0;0;0]; rotorAng = [0;0;0];
    for i = 1:4
        Fi   = P.n*T(i);
        Fsum = Fsum + Fi;
        tau  = tau + cross(P.pos(i,:)', Fi) - P.dir(i)*P.n*Qd(i);
        rotorAng = rotorAng + P.dir(i)*P.JR*Om(i)*[0;0;1];
    end

    % --- rotational dynamics with gyroscopic coupling --------------------------
    % wdot = Jb \ ( tau - omega x (Jb*omega + h_R) - tau_spin )
    %
    % SPIN-UP REACTION TORQUE (tau_spin) IS OMITTED ON PURPOSE.
    % Full term (would be subtracted inside the bracket below):
    %     tau_spin = sum_i th_i * JR * n_i * Omdot_i
    %     Omdot_i  = ( (kT/Rmot)*(Ebat*delta_i - ke*Omega_i) - MF - cQ*Omega_i^2 ) / JR
    % With n_i = [0;0;-1] it acts only on yaw: ~0 for a collective change
    % (CCW/CW cancel) and nonzero for a yaw manoeuvre.
    %
    % Why it is left out:
    %  - Numerically stiff. JR is tiny (6e-5), so when delta changes Omdot spikes
    %    to ~1e4 rad/s^2 for a few ms. That makes tau_spin ~0.5 Nm, ~10x the real
    %    aerodynamic yaw torque (~0.05 Nm), and at the 1 kHz RK4 step it drives the
    %    closed loop unstable (the vehicle tumbles).
    %  - Two physical fixes were tried and neither held: (a) a first-order ESC lag
    %    on delta (smooths the Omdot spike but, at tau_esc ~ the step size, did not
    %    cure the instability); (b) operator splitting (freeze tau_spin across the
    %    RK4 sub-steps) - it removed the sub-step inconsistency but the loop still
    %    diverged in full closed-loop missions.
    %  - Effect is small anyway: JR is small, the term is zero in steady flight and
    %    only transient during manoeuvres. Stephan likewise drops it in the control
    %    design model, so omitting it is the standard, defensible choice here.
    %  - To re-enable: compute tau_spin (per-step, frozen) and subtract it below.
    %    A stiff/implicit integrator or a slower-acceleration motor model would be
    %    needed first.
    wdot = P.Jb \ (tau - cross(w, P.Jb*w + rotorAng));

    % --- kinematics, translation (NED, gravity +z), position ---------------
    qdot = 0.5*Q.Omega(q)*w;
    vdot = (Q.toR(q)*Fsum)/P.m + [0;0;P.g];
    rdot = v;

    % --- battery Coulomb counting: Qdot = -Isig, in Ah (Isig is in A) ------
    Qdot = -Isig/3600;

    xdot = [Omdot; qdot; wdot; vdot; rdot; Qdot];
end

% Open-circuit voltage of the LiPo pack as a function of state of charge Q [Ah]
%   Ehat(Q) = E0 - Epol*(Q0/Q) + Eexp*exp(-(Q0-Q)/Qexp)   (Stephan Eq. 3.53)
% Q is floored to a small positive value to keep the 1/Q term finite near empty.
function Ehat = ocvFromCharge(Qbat, P)
    b    = P.bat;
    Qsafe = max(Qbat, 1e-3);
    Ehat  = b.E0 - b.Epol*(b.Q0/Qsafe) + b.Eexp*exp(-(b.Q0-Qsafe)/b.Qexp);
end

% Battery open-circuit voltage Ehat(Q), exposed for the controller's proactive
% voltage compensation (it scales PWM by the live pack voltage) and for logging.
% The load-dependent sag Rbat*Isig is applied inside quadDeriv, not here, since
% it depends on the instantaneous current the controller is about to command.
function Ebat = batteryVoltage(x, P)
    Ebat = ocvFromCharge(x(18), P);
end

% Pack current I_Sigma drawn from the battery at the current operating point,
% for logging. Same quasi-static current as quadDeriv (Eq. 3.33 / 3.42b):
%   Imot_i = max((Ehat*delta_i - ke*Omega_i)/Rmot, 0),  Isig = sum_i delta_i*Imot_i
function Isig = packCurrent(x, delta, P)
    Om   = x(1:4);
    delta = min(max(delta(:),0),1);
    Ehat = ocvFromCharge(x(18), P);
    Imot = max((Ehat*delta - P.ke*Om)/P.Rmot, 0);
    Isig = min(sum(delta.*Imot), P.bat.Imax);
end

function x = rk4(x, delta, P, dt)
    k1 = quadDeriv(x,           delta, P);
    k2 = quadDeriv(x + dt/2*k1, delta, P);
    k3 = quadDeriv(x + dt/2*k2, delta, P);
    k4 = quadDeriv(x + dt*k3,   delta, P);
    x  = x + dt/6*(k1 + 2*k2 + 2*k3 + k4);
    x(5:8) = x(5:8)/norm(x(5:8));        % renormalise quaternion
    x(18)  = max(x(18), 0);             % charge cannot go negative (empty pack)
end