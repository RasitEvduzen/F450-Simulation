function out = quadSysId()
%QUADSYSID  Closed-loop system identification with external chirp excitation.
%   out.run(P, dyn, ctrl, vis)    -> full sys-id flight: take off, hover, and
%                                    excite roll/pitch/yaw/altitude in turn;
%                                    estimate + report each, then animate.
%   out.chirp(P, opt)             -> [t, exc] log-swept excitation
%   out.estimate(u, y, na, nb, nk)-> ARX model (batch least squares)
%   out.coherence(u, y, dt, nfft) -> [f, Cxy] magnitude-squared coherence
%   out.report(axisName, model, u, y, P) -> per-axis analysis figure
%
%   Scheme: the vehicle flies closed-loop (full position hold at hover) while
%   a chirp is added on top of one setpoint. Identifying in closed loop
%   without external excitation fails because the input is correlated with
%   the controller's response to noise; the injected chirp breaks this
%   correlation. Each input->output path is fit with a batch-LS ARX model:
%       A(q) y[k] = B(q) u[k].
%
%   Axes (inner -> outer):
%     roll  : torque tau_x -> roll rate p
%     pitch : torque tau_y -> pitch rate q
%     yaw   : torque tau_z -> yaw rate r
%     alt   : collective thrust Fz -> vertical velocity w
%   (altitude uses thrust->w, a single integrator + motor lag, rather than
%    thrust->z which is a double integrator and ill-conditioned to identify.)

    out.run       = @runSysId;
    out.chirp     = @chirpExcitation;
    out.estimate  = @arxEstimate;
    out.coherence = @coherence;
    out.report    = @report;
end

%% ------------------------------------------------------------------------
function runSysId(P, dyn, ctrl, vis)
    axes = { 'roll',     struct('Tend',15,'f0',0.5,'f1',40,'amp',deg2rad(40));
             'pitch',    struct('Tend',15,'f0',0.5,'f1',40,'amp',deg2rad(40));
             'yaw',      struct('Tend',15,'f0',0.5,'f1',30,'amp',deg2rad(45));
             'altitude', struct('Tend',15,'f0',1.0,'f1',20,'amp',0.05) };
    z_hover  = -2.0;  T_take = 4.0;  T_settle = 1.5;

    Binv = ctrl.mixer(P);
    x = zeros(17,1); x(1:4)=P.Om_hover; x(5)=1;
    s.rate_int = zeros(3,1); s.w_prev = zeros(3,1);
    s.vel_int  = zeros(3,1); s.v_prev = zeros(3,1);
    pos_sp = [0; 0; z_hover];
    Xall = [];  phases = {};      % per-step state log and animation labels

    % take off + settle
    for k = 1:round(T_take/P.dt)
        [delta, s] = hoverStep(x, pos_sp, 0, '', 'settle', P, Binv, s);
        x = dyn.rk4(x, delta, P, P.dt);  Xall = [Xall, x];
        phases{end+1} = 'Take-off / hover';
    end

    % excite each axis
    for a = 1:size(axes,1)
        name = axes{a,1};
        [~, exc] = chirpExcitation(P, axes{a,2});  Nc = numel(exc);
        for k = 1:round(T_settle/P.dt)      % settle
            [delta, s] = hoverStep(x, pos_sp, 0, '', 'settle', P, Binv, s);
            x = dyn.rk4(x, delta, P, P.dt);  Xall = [Xall, x];
            phases{end+1} = sprintf('Settle (before %s)', name);
        end
        U = zeros(Nc,1);  Y = zeros(Nc,1);  % chirp segment
        for k = 1:Nc
            [delta, s, uk, yk] = hoverStep(x, pos_sp, exc(k), name, 'chirp', P, Binv, s);
            x = dyn.rk4(x, delta, P, P.dt);  Xall = [Xall, x];
            U(k) = uk;  Y(k) = yk;
            phases{end+1} = sprintf('Sys-ID: %s axis', name);
        end
        model = arxEstimate(U, Y, 2, 2, 1);
        report(name, model, U, Y, P);
    end

    % animate the whole sys-id flight (single wide view, active phase shown).
    % the chirp wobble is small and fast, so use a high frame rate; the drone
    % holds position, so a single follow view is clearer than a dual view.
    Nall = size(Xall,2);  t = (0:Nall-1)'*P.dt;
    REF  = [zeros(Nall,2), z_hover*ones(Nall,1), zeros(Nall,1)];
    propAng = zeros(4,Nall); pa = zeros(4,1);
    for k = 1:Nall, pa = pa + P.dir.*Xall(1:4,k)*P.dt;  propAng(:,k) = pa; end
    Panim = P;
    Panim.vis.step   = 50;          %  step
    Panim.vis.single = true;        % single wide view (no whole-trajectory)
    Panim.vis.viewR  = 0.7;
    vis.animate(t, Xall, REF, propAng, Panim, phases);
end

%% ---- closed-loop hover step with chirp injection -----------------------
% Full position hold at pos_sp; the chirp is ADDED on top of the hold command
% so the drone stays in place while being excited. On a rate axis it adds to
% the rate setpoint (input=torque, output=body rate); on 'altitude' it adds to
% the collective (input=Fz, output=vertical velocity w). State s carries the
% integrators and previous values between calls.
function [delta, s, uk, yk] = hoverStep(x, pos_sp, exc, ax, phase, P, Binv, s)
    Q = quatUtils();
    q = x(5:8); w = x(9:11); vel = x(12:14); pos = x(15:17);
    uk = 0; yk = 0;  idxmap = struct('roll',1,'pitch',2,'yaw',3);

    % position -> velocity -> acceleration -> attitude (level + altitude hold)
    vel_sp = (pos_sp - pos).*P.Kp_pos;
    vel_sp(1:2) = max(min(vel_sp(1:2), P.vmax_xy), -P.vmax_xy);
    vel_sp(3)   = max(min(vel_sp(3),   P.vmax_dn), -P.vmax_up);
    vel_err = vel_sp - vel;
    vel_dot = (vel - s.v_prev)/P.dt;  s.v_prev = vel;
    acc_sp  = vel_err.*P.Kp_vel + s.vel_int - vel_dot.*P.Kd_vel;
    s.vel_int = s.vel_int + vel_err.*P.Ki_vel*P.dt;
    s.vel_int(3) = max(min(s.vel_int(3), P.g), -P.g);

    z_sf   = -P.g + acc_sp(3);
    body_z = [-acc_sp(1); -acc_sp(2); -z_sf];  body_z = body_z/norm(body_z);
    thr_z  = acc_sp(3)*(P.hover_thr/P.g) - P.hover_thr;
    coll   = min(thr_z/body_z(3), -0.12);
    qz = Q.between([0;0;-1], -body_z/norm(body_z));
    qd = Q.mul(qz, [1;0;0;0]);  qd = qd/norm(qd);
    rate_sp = attitudeHold(q, qd, P, Q);

    if strcmp(phase,'chirp')                 % inject chirp on top of hold
        if strcmp(ax,'altitude'), coll = coll + exc;
        else, idx = idxmap.(ax);  rate_sp(idx) = rate_sp(idx) + exc; end
    end

    rate_err  = rate_sp - w;
    ang_accel = (w - s.w_prev)/P.dt;  s.w_prev = w;
    torque = P.Kp_rate.*rate_err + s.rate_int - P.Kd_rate.*ang_accel;
    i_factor = max(0, 1 - (rate_err/deg2rad(400)).^2);
    s.rate_int = s.rate_int + i_factor.*P.Ki_rate.*rate_err*P.dt;
    s.rate_int = max(min(s.rate_int, P.int_lim), -P.int_lim);

    Fz = coll/P.hover_thr*P.m*P.g;
    Tm = max(Binv*[Fz; torque], 0.01);
    Om = sqrt(Tm/P.cT);
    delta = max(min((Om + P.cQ*Om.^2/P.kT)/P.Ehat, 1), 0);

    if strcmp(phase,'chirp')
        if strcmp(ax,'altitude'), uk = Fz; yk = x(14);
        else, idx = idxmap.(ax);  uk = torque(idx);  yk = x(8+idx); end
    end
end

function rate_sp = attitudeHold(q, qd, P, Q)
    q = q/norm(q); qd = qd/norm(qd);
    e_z = Q.dcm_z(q); e_z_d = Q.dcm_z(qd);
    qd_red = Q.between(e_z, e_z_d);
    if abs(qd_red(2)) > 1-1e-5 || abs(qd_red(3)) > 1-1e-5
        qd_full = qd; else, qd_full = Q.mul(qd_red, q); end
    q_mix = Q.canon(Q.mul(Q.inv(qd_full), qd));
    q_mix(1) = max(min(q_mix(1),1),-1); q_mix(4) = max(min(q_mix(4),1),-1);
    qd_m = Q.mul(qd_full, [cos(P.yaw_w*acos(q_mix(1))); 0; 0; sin(P.yaw_w*asin(q_mix(4)))]);
    qe = Q.canon(Q.mul(Q.inv(q), qd_m));
    rate_sp = max(min(2*qe(2:4).*P.Kp_att, P.rate_max), -P.rate_max);
end

%% ---- chirp / estimation / coherence ------------------------------------
function [t, exc] = chirpExcitation(P, opt)
    if nargin < 2, opt = struct; end
    dt   = P.dt;
    Tend = getf(opt,'Tend',15);  f0 = getf(opt,'f0',0.5);
    f1   = getf(opt,'f1',40);    amp = getf(opt,'amp',deg2rad(60));
    t   = (0:dt:Tend)';
    phi = 2*pi*(f0*t + (f1-f0)/(2*Tend)*t.^2);
    exc = amp*sin(phi);
end

function model = arxEstimate(u, y, na, nb, nk)
    if nargin < 5, nk = 1; end
    u = u(:); y = y(:); N = numel(y);
    p0 = max(na, nb+nk-1);  rows = N - p0;
    Phi = zeros(rows, na+nb);  Y = zeros(rows,1);
    for i = 1:rows
        k = i + p0;
        Phi(i,1:na)     = -y(k-1:-1:k-na)';
        Phi(i,na+1:end) =  u(k-nk:-1:k-nk-nb+1)';
        Y(i)            =  y(k);
    end
    theta = Phi \ Y;
    model.A     = [1, theta(1:na)'];
    model.B     = [zeros(1,nk), theta(na+1:end)'];
    model.na    = na;  model.nb = nb;  model.nk = nk;
    model.poles = roots(model.A);
    model.fit   = 100*(1 - norm(Y-Phi*theta)/norm(Y-mean(Y)));
end

function [f, Cxy] = coherence(u, y, dt, nfft)
    if nargin < 4, nfft = 2048; end
    u = u(:) - mean(u);  y = y(:) - mean(y);
    win = hann(nfft);  step = nfft - floor(nfft/2);
    Puu = 0; Pyy = 0; Puy = 0;
    for st = 1:step:(numel(u)-nfft)
        U = fft(u(st:st+nfft-1).*win);  Y = fft(y(st:st+nfft-1).*win);
        Puu = Puu + abs(U).^2;  Pyy = Pyy + abs(Y).^2;  Puy = Puy + conj(U).*Y;
    end
    Cxy = abs(Puy).^2 ./ (Puu.*Pyy + 1e-12);
    nb2 = floor(nfft/2)+1;  Cxy = Cxy(1:nb2);
    f = (0:nb2-1)'/(nfft*dt);
end

%% ---- per-axis analysis figure ------------------------------------------
function report(axisName, model, u, y, P)
    u = u(:); y = y(:); N = numel(y);  t = (0:N-1)'*P.dt;
    isRate = ~strcmpi(axisName,'altitude');
    sc = 1; yunit = 'm/s'; uunit = 'N';
    if isRate, sc = 180/pi; yunit = 'deg/s'; uunit = 'N\cdotm'; end

    % one-step-ahead prediction
    p0 = max(model.na, model.nb+model.nk-1);
    yhat = zeros(N,1);  yhat(1:p0) = y(1:p0);
    a = model.A(2:end);  b = model.B(model.nk+1:end);
    for k = p0+1:N
        yhat(k) = -a*y(k-1:-1:k-model.na) + b*u(k-model.nk:-1:k-model.nk-model.nb+1);
    end

    % frequency response
    f = logspace(-1, log10(0.5/P.dt), 500);
    z = exp(1j*2*pi*f*P.dt);
    H = polyval(model.B,1./z) ./ polyval(model.A,1./z);
    [fc, Cxy] = coherence(u, y, P.dt);

    cRef=[0 0 0]; cAct=[0.85 0.15 0.15]; cIn=[0.20 0.40 0.70];
    figure('Name',sprintf('Sys-ID: %s',axisName),'Color','w', ...
           'units','normalized','outerposition',[0 0 1 1]);
    tl = tiledlayout(3,2,'Padding','compact','TileSpacing','compact');

    % (1) input
    nexttile; hold on; grid on; box on;
    plot(t,u,'-','Color',cIn,'LineWidth',1.0);
    ylabel(sprintf('input [%s]',uunit)); xlabel('time [s]');
    title('Excitation input (chirp)'); legend({'u'},'Location','best'); xlim([0 t(end)]);

    % (2) output + one-step prediction
    nexttile; hold on; grid on; box on;
    plot(t,sc*y,'-','Color',cRef,'LineWidth',1.0);
    plot(t,sc*yhat,'-','Color',cAct,'LineWidth',1.2);
    ylabel(sprintf('output [%s]',yunit)); xlabel('time [s]');
    title(sprintf('Output & ARX prediction  (fit %.1f%%)',model.fit));
    legend({'measured','ARX model'},'Location','best'); xlim([0 t(end)]);

    % (3) Bode magnitude
    nexttile; hold on; grid on; box on;
    semilogx(f,20*log10(abs(H)),'-','Color',cAct,'LineWidth',1.4);
    set(gca,'XScale','log'); ylabel('|H| [dB]'); xlabel('frequency [Hz]');
    title('Bode - magnitude');

    % (4) Bode phase
    nexttile; hold on; grid on; box on;
    semilogx(f,rad2deg(unwrap(angle(H))),'-','Color',cIn,'LineWidth',1.4);
    set(gca,'XScale','log'); ylabel('phase [deg]'); xlabel('frequency [Hz]');
    title('Bode - phase');

    % (5) coherence
    nexttile; hold on; grid on; box on;
    semilogx(fc, Cxy, '-','Color',cAct,'LineWidth',1.4);
    set(gca,'XScale','log'); ylabel('\gamma^2'); xlabel('frequency [Hz]');
    title('Input-output coherence'); ylim([0 1.05]); xlim([fc(2) fc(end)]);

    % (6) identified transfer function (text)
    nexttile; axis off;
    txt = tfText(axisName, model, P, uunit, yunit);
    text(0.02, 0.98, txt, 'Units','normalized', 'VerticalAlignment','top', ...
         'FontName','FixedWidth', 'FontSize',11, 'Interpreter','none');

    title(tl, sprintf('System identification: %s axis  (ARX %dp%dz, delay %d)', ...
          axisName, model.na, model.nb, model.nk), 'FontWeight','bold');
end

% build a readable transfer-function description block
function txt = tfText(axisName, model, P, uunit, yunit)
    poly2s = @(c) strjoin(arrayfun(@(i) sprintf('%+.4g q^-%d',c(i),i-1), ...
                 1:numel(c), 'UniformOutput',false), ' ');
    pls = sprintf('%.4f  ', model.poles);
    Ts  = P.dt;
    txt = {
        sprintf('Identified discrete TF   H(q) = B(q)/A(q),   Ts = %.4f s', Ts)
        sprintf('input u [%s]  ->  output y [%s]', uunit, yunit)
        ''
        ['A(q) = ' poly2s(model.A)]
        ['B(q) = ' poly2s(model.B)]
        ''
        sprintf('poles : %s', pls)
        sprintf('fit   : %.2f %%  (one-step prediction)', model.fit)
        ''
        'Notes:'
        '  pole near z=1   -> integrator (rate/velocity)'
        sprintf('  pole near %.3f -> motor lag (z = exp(-Ts/tau_mot))', exp(-Ts/(P.JR/P.kT)))
    };
    txt = strjoin(txt, newline);
end

function v = getf(s,f,d), if isfield(s,f), v=s.(f); else, v=d; end, end