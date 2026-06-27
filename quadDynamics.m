function out = quadDynamics()
%QUADDYNAMICS  Nonlinear flight-mechanical model + RK4 integrator.
%   Returns a struct of handles:
%     out.deriv(x,delta,P)   -> state derivative (Stephan Eq. 3.101)
%     out.rk4(x,delta,P,dt)  -> one RK4 step (with quaternion renormalise)
%
%   State (17): x = [Omega(4); q(4); omega_body(3); v_g(3); r_g(3)]
%     Omega      - rotor speeds [rad/s]
%     q          - attitude quaternion (body->geodetic), scalar-first
%     omega_body - body angular rate (p,q,r) [rad/s]
%     v_g, r_g   - velocity / position in geodetic NED [m/s, m]
%   Input (4): delta in [0,1] - per-rotor PWM command.
%
%   EQUATIONS OF MOTION (Stephan "Multicopter Flight Control", Eq. 3.101):
%
%     Motor:        JR*Omdot_i = kT*(Ehat*delta_i - Omega_i) - Q_i - MF
%     Kinematics:   qdot       = 1/2 * Theta(q) * omega
%     Rotational:   Jb*omdot   = tau_Sigma
%                                - omega x ( Jb*omega + sum_i th_i*JR*n_i*Omega_i )
%     Translation:  m*vdot     = R(q)*F_Sigma + m*g
%     Position:     rdot       = v
%
%   with per-rotor thrust/drag and the resulting body force/torque
%     T_i = cT*Omega_i^2,   Q_i = cQ*Omega_i^2,
%     tau_Sigma = sum_i  p_i x (n_i*T_i)  -  th_i*n_i*Q_i,
%     F_Sigma   = sum_i  n_i*T_i,
%   where n_i is the rotor thrust axis, p_i its position, th_i = +/-1 its spin
%   direction, and g = [0;0;+g] (NED, gravity points +z / down).
%
%   MODELLING ASSUMPTIONS:
%     1. Rigid body, fixed mass and inertia Jb (no fuel/payload change).
%     2. Near-hover / low-speed flight: undisturbed rotor inflow V_inf = 0,
%        so the in-plane force H_i = 0 and T_i, Q_i reduce to quadratic
%        polars in Omega_i (Eq. 3.104).
%     3. Out of ground effect (cGE = 1).
%     4. No airframe/rotor aerodynamic drag (C = C_(.) = 0) and no wind
%        (v_w = omega_w = 0); no external disturbance torques/forces.
%     5. Rotors horizontal, thrust axis vertical in body frame: n_i = [0;0;-1].
%     6. Rotor spin axis taken as +z in the gyroscopic term (consistent with 5).
%     7. Gyroscopic momentum coupling (omega x sum th_i*JR*n_i*Omega_i) is KEPT;
%        the spin-up reaction torque (sum th_i*JR*n_i*Omdot_i) is OMITTED,
%        consistent with the control-design model (JR small).
%     8. Battery voltage Ehat held constant (no discharge dynamics, Eq. 3.103).

    out.deriv = @quadDeriv;
    out.rk4   = @rk4;
end

function xdot = quadDeriv(x, delta, P)
    Q  = quatUtils();
    Om = x(1:4); q = x(5:8); w = x(9:11); v = x(12:14);
    q  = q/norm(q);

    % motor dynamics: JR*Omdot = kT*(Ehat*delta - Omega) - Q - MF
    Qd    = P.cQ*Om.^2;                       % rotor drag torque Q_i
    Omdot = (P.kT*(P.Ehat*delta - Om) - Qd - P.MF)/P.JR;

    % per-rotor thrust T_i, summed body force F_Sigma and torque tau_Sigma,
    % and rotor angular momentum sum_i th_i*JR*n_i*Omega_i (gyroscopic term)
    T = P.cT*Om.^2;
    Fsum = [0;0;0]; tau = [0;0;0]; rotorAng = [0;0;0];
    for i = 1:4
        Fi   = P.n*T(i);
        Fsum = Fsum + Fi;
        tau  = tau + cross(P.pos(i,:)', Fi) - P.dir(i)*P.n*Qd(i);
        rotorAng = rotorAng + P.dir(i)*P.JR*Om(i)*[0;0;1];
    end

    % rotational dynamics with gyroscopic coupling (reaction torque omitted)
    wdot = P.Jb \ (tau - cross(w, P.Jb*w + rotorAng));

    % kinematics, translation (NED, gravity +z), position
    qdot = 0.5*Q.Omega(q)*w;
    vdot = (Q.toR(q)*Fsum)/P.m + [0;0;P.g];
    rdot = v;

    xdot = [Omdot; qdot; wdot; vdot; rdot];
end

function x = rk4(x, delta, P, dt)
    k1 = quadDeriv(x,           delta, P);
    k2 = quadDeriv(x + dt/2*k1, delta, P);
    k3 = quadDeriv(x + dt/2*k2, delta, P);
    k4 = quadDeriv(x + dt*k3,   delta, P);
    x  = x + dt/6*(k1 + 2*k2 + 2*k3 + k4);
    x(5:8) = x(5:8)/norm(x(5:8));   % renormalise quaternion
end