%% F450 quadcopter system identification
%  Modules:
%     quadParams      - configuration (model, gains, visuals)
%     quadDynamics    - nonlinear model (Eq. 3.101) + RK4
%     quadControl     - multi-rate cascaded controller + mixer
%     quadVisualize   - STL animation + analysis figures
%     quadSysId       - closed-loop chirp excitation + ARX identification
%     quatUtils       - quaternion helpers
%
%  Closed-loop sys-id: take off, hold hover, and excite roll/pitch/yaw/altitude
%  in turn with a chirp added on top of the hold command. Each axis is fit with
%  a batch-LS ARX model and reported (I/O, fit, poles, Bode, coherence, TF).
%  STL files (Base / MotorPropCW / MotorPropCCW) must be on the path.
% Written By: Rasit Evduzen
% Date: 26.06.2026
%%
clc; clear; close all;
addpath('src');

% --- load modules ---
P    = quadParams();
dyn  = quadDynamics();
ctrl = quadControl('cascade');
vis  = quadVisualize();
sid  = quadSysId();

% --- run the full sys-id flight (excite, estimate, report, animate) ---
sid.run(P, dyn, ctrl, vis);