function P = quadParams()
%QUADPARAMS  All configuration for the F450 quadcopter simulator: physical
%   model, baseline controller gains, simulation rates and visual settings.
%   Everything tunable lives here.

    %% Physical model (Stephan "Multicopter Flight Control", Eq. 3.101)
    P.m   = 1.5;                      % mass [kg]
    P.g   = 9.81;                     % gravity [m/s^2]
    P.Jb  = diag([0.02 0.02 0.04]);   % body inertia [kg m^2]

    P.JR   = 6e-5;                    % rotor + motor inertia [kg m^2]

    %% Motor / ESC (BLDC, Stephan Eqs. 3.31-3.36)
    %  Motor dynamics:  JR*Omdot = (kT/Rmot)*(E*delta - ke*Omega) - MF - Q
    %  ESC voltage:     E_mot = delta*E,  E = battery terminal voltage Ehat(Q)
    %  Armature current is taken quasi-static (Eq. 3.33): the full electrical
    %  dynamics is  L*dI/dt = E*delta - ke*Omega - Rmot*I, but the electrical
    %  time constant L/Rmot (~0.3 ms) is ~1000x faster than the rotor (~600 ms),
    %  so dI/dt = 0 is assumed and current is algebraic. L is kept as a parameter
    %  for documentation only (=0); making it >0 turns I into a stiff state that
    %  needs an implicit solver, not the 1 kHz RK4 used here.
    P.kT   = 0.02894;                 % torque constant [Nm/A] (kV ~ 330 RPM/V)
    P.ke   = P.kT;                    % back-EMF constant [V/(rad/s)]; ke = kT (BLDC), edit if real data differs
    P.Rmot = 0.08;                    % motor resistance [Ohm]
    P.Lmot = 0.0;                     % armature inductance [H]; modeled as 0 (quasi-static current, see above)
    P.MF   = 0.0;                     % motor friction torque [Nm]
    P.cT   = 1.5e-5;                  % thrust coeff:  T = cT*Omega^2
    P.cQ   = 2.5e-7;                  % drag   coeff:  Q = cQ*Omega^2

    %% Battery (LiPo discharge model, Stephan Eqs. 3.52-3.54)
    %  Test pack: 8S / 5Ah / 100C. State of Charge Q [Ah] is an extra plant
    %  state; terminal voltage follows the open-circuit curve Ehat(Q) minus the
    %  internal-resistance sag Rbat*I_Sigma. Coulomb counting: Qdot = -I_Sigma.
    %    Ehat(Q) = E0 - Epol*(Q0/Q) + Eexp*exp(-(Q0-Q)/Qexp)
    %    E = Ehat(Q) - Rbat*I_Sigma
    P.bat.NS      = 8;                % series cells
    P.bat.Q0      = 5.0;             % capacity / full charge [Ah]
    P.bat.Crating = 100;            % C-rating
    P.bat.Imax    = P.bat.Q0*P.bat.Crating;   % max discharge current [A] (=500)
    P.bat.E0      = P.bat.NS*3.8;    % nominal voltage [V] (plateau, ~3.8 V/cell)
    P.bat.Epol    = P.bat.NS*0.045;  % polarisation voltage [V] (end-of-charge sag)
    P.bat.Eexp    = P.bat.NS*0.30;   % exponential-zone overshoot [V] (full-charge peak)
    P.bat.Qexp    = 0.35;           % exponential-zone inverse decay [Ah]
    P.bat.Rbat    = 0.010;          % internal resistance [Ohm] (voltage sag)
    P.bat.Q_init  = P.bat.Q0;        % initial charge [Ah] (start full)
    
    %% Geometry: F450 X-config
    % M1 front-right (CCW), M2 rear-left (CCW), M3 front-left (CW), M4 rear-right (CW)
    P.arm = 0.225;                    % rotor radius from CoG [m] (450 mm span)
    az    = deg2rad([45 225 315 135]);              % azimuths for M1..M4
    P.pos = [P.arm*cos(az'), P.arm*sin(az'), zeros(4,1)];
    P.dir = [+1; +1; -1; -1];         % spin direction (+1 CCW, -1 CW)
    P.n   = [0; 0; -1];               % thrust axis (body -z = up)

    %% Hover trim (derived)
    P.Om_hover    = sqrt(P.m*P.g/(4*P.cT));
    % full-charge open-circuit voltage, then invert the motor model (Eq. 3.36,
    % stationary) for the PWM that holds Om_hover: used only as the held initial
    % command; the controller recomputes delta each step at the live voltage.
    Ehat_full     = P.bat.E0 - P.bat.Epol + P.bat.Eexp;   % Ehat at Q = Q0
    P.delta_hover = (P.kT*P.Om_hover + P.Rmot*(P.cQ*P.Om_hover^2 + P.MF)/P.kT)/Ehat_full;

    %% Controller gains (baseline cascade)
    P.Kp_rate = [0.15; 0.15; 0.20];   % rate loop
    P.Ki_rate = [0.20; 0.20; 0.10];
    P.Kd_rate = [0.003; 0.003; 0.0];
    P.Kff_rate= [0; 0; 0];
    P.int_lim = 0.30;

    P.Kp_att  = [6.5; 6.5; 2.8];      % attitude loop
    P.rate_max= deg2rad([220; 220; 200]);
    P.yaw_w   = 0.4;

    P.Kp_pos  = [1.5; 1.5; 2.0];      % position loop
    P.Kp_vel  = [2.5; 2.5; 4.0];      % velocity loop
    P.Ki_vel  = [0.4; 0.4; 2.0];
    P.Kd_vel  = [0.2; 0.2; 0.0];
    P.hover_thr = 0.5;
    P.vmax_xy = 12;                   % velocity limits [m/s]
    P.vmax_up = 3;
    P.vmax_dn = 1.5;

    %% Simulation (multi-rate)
    % Physics integrates at f_phys; each loop runs at its own rate and holds
    % its output zero-order between updates, like a real flight stack.
    P.f_phys = 1000;                  % physics / integration [Hz]
    P.f_rate = 1000;                  % rate (inner) loop
    P.f_att  = 250;                   % attitude loop
    P.f_pos  = 50;                    % position/velocity (outer) loop
    P.dt     = 1/P.f_phys;            % integration step [s]

    %% Trajectory
    % Per-axis translation profile for the waypoint mission (out.build):
    %   false = 3-segment trapezoidal velocity (acceleration steps)
    %   true  = 7-segment S-curve (finite jerk, smoother acceleration)
    P.useDoubleScurve = true;

    %% Visualisation
    P.stl.base    = 'stl/Base.STL';
    P.stl.propCW  = 'stl/MotorPropCW.STL';
    P.stl.propCCW = 'stl/MotorPropCCW.STL';
    P.stl.scale   = 1e-3;             % STL mm -> model m
    P.stl.flip    = [1 0 0; 0 -1 0; 0 0 -1];  % mount prop the right way up
    P.stl.propZ   = -0.0185;          % display-only prop height above body [m]
    P.vis.viewR   = 0.5;              % follow-cam half-window [m]
    P.vis.triLen  = 0.25;             % body triad arm length [m]
    P.vis.step    = 200;              % draw one frame every N physics steps
    P.vis.view    = [135 22];         % [azimuth elevation]
end