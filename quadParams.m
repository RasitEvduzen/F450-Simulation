function P = quadParams()
%QUADPARAMS  All configuration for the F450 quadcopter simulator: physical
%   model, baseline controller gains, simulation rates and visual settings.
%   Everything tunable lives here.

    %% Physical model (Stephan "Multicopter Flight Control", Eq. 3.101)
    P.m   = 1.5;                      % mass [kg]
    P.g   = 9.81;                     % gravity [m/s^2]
    P.Jb  = diag([0.02 0.02 0.04]);   % body inertia [kg m^2]

    P.JR   = 6e-5;                    % rotor + motor inertia [kg m^2]
    P.kT   = 0.003;                   % motor torque constant (tau_mot = JR/kT)
    P.Ehat = 1031.3;                  % voltage term (tuned so hover delta ~0.5)
    P.MF   = 0.0;                     % motor friction
    P.cT   = 1.5e-5;                  % thrust coeff:  T = cT*Omega^2
    P.cQ   = 2.5e-7;                  % drag   coeff:  Q = cQ*Omega^2

    %% Geometry: F450 X-config
    % M1 front-right (CCW), M2 rear-left (CCW), M3 front-left (CW), M4 rear-right (CW)
    P.arm = 0.225;                    % rotor radius from CoG [m] (450 mm span)
    az    = deg2rad([45 225 315 135]);              % azimuths for M1..M4
    P.pos = [P.arm*cos(az'), P.arm*sin(az'), zeros(4,1)];
    P.dir = [+1; +1; -1; -1];         % spin direction (+1 CCW, -1 CW)
    P.n   = [0; 0; -1];               % thrust axis (body -z = up)

    %% Hover trim (derived)
    P.Om_hover    = sqrt(P.m*P.g/(4*P.cT));
    P.delta_hover = (P.Om_hover + P.cQ*P.Om_hover^2/P.kT)/P.Ehat;

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

    %% Visualisation
    P.stl.base    = 'Base.STL';
    P.stl.propCW  = 'MotorPropCW.STL';
    P.stl.propCCW = 'MotorPropCCW.STL';
    P.stl.scale   = 1e-3;             % STL mm -> model m
    P.stl.flip    = [1 0 0; 0 -1 0; 0 0 -1];  % mount prop the right way up
    P.stl.propZ   = -0.0185;          % display-only prop height above body [m]
    P.vis.viewR   = 0.5;              % follow-cam half-window [m]
    P.vis.triLen  = 0.25;             % body triad arm length [m]
    P.vis.step    = 200;              % draw one frame every N physics steps
    P.vis.view    = [135 22];         % [azimuth elevation]
end
