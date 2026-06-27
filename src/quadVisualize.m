function out = quadVisualize()
%QUADVISUALIZE  Flight visualisation: STL animation + analysis figures.
%   out.animate(t, X, REF, propAng, P)  live STL flight, dual view
%                                        (left follow-cam, right whole trajectory)
%   out.figures(t, X, REF, P)           position / altitude / error / heading
%   out.states (t, X, REF, LOG, P)      velocity / attitude / rate vs command
%   out.motors (t, LOG, P)              motor thrusts + mixer commands
%
%   STL technique: stlread -> transform point cloud (R*scale*pts + pos)
%   -> trisurf -> camlight.

    out.animate = @animate;
    out.figures = @figures;
    out.states  = @states;
    out.motors  = @motors;
end

function animate(t, X, REF, propAng, P, phases)
%ANIMATE  STL flight animation. Optional `phases` is a 1xN cellstr labelling
%   each step (e.g. the active sys-id axis); when given, the camera title
%   shows the current phase so you can tell which axis is being excited.
%   If P.vis.single is true, only the follow-cam is drawn (single wide view);
%   otherwise a dual view (follow-cam + whole trajectory) is shown.
    if nargin < 6, phases = {}; end
    single = isfield(P.vis,'single') && P.vis.single;
    Q = quatUtils();
    baseGeo = stlread(P.stl.base);
    propCW  = stlread(P.stl.propCW);
    propCCW = stlread(P.stl.propCCW);
    propSel = {propCCW, propCCW, propCW, propCW};   % M1,M2 CCW ; M3,M4 CW
    % motor positions come straight from P.pos (same order/numbering as the
    % model & mixer); only the display z-offset is added.
    motorPos = P.pos + [0 0 P.stl.propZ];

    N = numel(t);
    fig = figure('Name','F450 Flight','Color','w','units','normalized','outerposition',[0 0 1 1]);
    % wide-view fixed limits (whole mission, equal scale)
    allP = [REF(:,1:3); X(15:17,:)'];
    mn=min(allP,[],1); mx=max(allP,[],1); ctr=(mn+mx)/2; sp=max(mx-mn)+1.0;
    wideLim=[ctr(1)-sp/2 ctr(1)+sp/2 ctr(2)-sp/2 ctr(2)+sp/2 ctr(3)-sp/2 ctr(3)+sp/2];

    nv = 1 + ~single;        % 1 view (single) or 2 (dual)

    % --- follow cam (left, or the whole figure when single) ---
    if single, axF = axes; else, axF = subplot(1,2,1); end
    hold(axF,'on'); grid(axF,'on'); axis(axF,'equal');
    view(axF,P.vis.view(1),P.vis.view(2)); set(axF,'ZDir','reverse');
    xlabel(axF,'N [m]'); ylabel(axF,'E [m]'); zlabel(axF,'Up [m]');
    title(axF,'Follow camera'); camlight(axF,'headlight'); lighting(axF,'gouraud');
    hTitleF = get(axF,'Title');
    plot3(axF,REF(:,1),REF(:,2),REF(:,3),'--','Color',[0 0 0],'LineWidth',1);
    hTrailF=plot3(axF,nan,nan,nan,'-','Color',[0.85 0.15 0.15],'LineWidth',1.4);
    plot3(axF,REF(1,1),REF(1,2),REF(1,3),'r.','MarkerFaceColor','r','MarkerSize',20,'HandleVisibility','off');

    % --- wide view (right) -- dual mode only ---
    if ~single
    axW = subplot(1,2,2); hold(axW,'on'); grid(axW,'on'); axis(axW,'equal');
    view(axW,P.vis.view(1),P.vis.view(2)); set(axW,'ZDir','reverse'); axis(axW,wideLim);
    xlabel(axW,'N [m]'); ylabel(axW,'E [m]'); zlabel(axW,'Up [m]');
    title(axW,'Whole trajectory'); camlight(axW,'headlight'); lighting(axW,'gouraud');
    plot3(axW,REF(:,1),REF(:,2),REF(:,3),'--','Color',[0 0 0],'LineWidth',1.2);
    hTrailW=plot3(axW,nan,nan,nan,'-','Color',[0.85 0.15 0.15],'LineWidth',1.4);
    % ground marks for waypoints
    plot3(axW,REF(1,1),REF(1,2),REF(1,3),'r.','MarkerFaceColor','r','MarkerSize',20,'HandleVisibility','off');
    end

    hB=gobjects(1,2); hP=gobjects(4,2); hT=gobjects(3,2);
    trail=nan(3,N); triCols=[1 0 0;0 1 0;0 0 1];

    for k=1:P.vis.step:N
        q=X(5:8,k); pos=X(15:17,k); R=Q.toR(q); trail(:,k)=pos;
        Pb=(R*(baseGeo.Points'*P.stl.scale))'+pos';
        for col=1:nv
            if col==1, ax=axF; else, ax=axW; end
            if ishandle(hB(col)), delete(hB(col)); end
            hB(col)=trisurf(baseGeo.ConnectivityList,Pb(:,1),Pb(:,2),Pb(:,3), ...
                    'FaceColor',[0.22 0.22 0.26],'EdgeColor','none','Parent',ax,'HandleVisibility','off');
            for i=1:4
                ca=cos(propAng(i,k)); sa=sin(propAng(i,k));
                Rspin=[ca -sa 0; sa ca 0; 0 0 1]; gp=propSel{i};
                Plocal=(Rspin*P.stl.flip*(gp.Points'*P.stl.scale))';
                Pw=(R*(Plocal+motorPos(i,:))')'+pos';
                if ishandle(hP(i,col)), delete(hP(i,col)); end
                c=[0.85 0.25 0.25]; if i==1||i==2, c=[0.2 0.45 0.9]; end  % CCW (M1,M2) blue, CW (M3,M4) red
                hP(i,col)=trisurf(gp.ConnectivityList,Pw(:,1),Pw(:,2),Pw(:,3), ...
                        'FaceColor',c,'EdgeColor','none','Parent',ax,'HandleVisibility','off');
            end
            for a=1:3
                v=R(:,a)*P.vis.triLen;
                if ishandle(hT(a,col)), delete(hT(a,col)); end
                hT(a,col)=plot3(ax,[pos(1) pos(1)+v(1)],[pos(2) pos(2)+v(2)], ...
                        [pos(3) pos(3)+v(3)],'Color',triCols(a,:),'LineWidth',2.5,'HandleVisibility','off');
            end
        end
        set(hTrailF,'XData',trail(1,1:k),'YData',trail(2,1:k),'ZData',trail(3,1:k));
        if ~single
            set(hTrailW,'XData',trail(1,1:k),'YData',trail(2,1:k),'ZData',trail(3,1:k));
        end
        if ~isempty(phases)
            set(hTitleF,'String',phases{k});
        end
        axis(axF,[pos(1)-P.vis.viewR pos(1)+P.vis.viewR pos(2)-P.vis.viewR ...
                  pos(2)+P.vis.viewR pos(3)-P.vis.viewR pos(3)+P.vis.viewR]);
        drawnow limitrate;
    end
end

function figures(t, X, REF, ~)
    Q = quatUtils();
    N = numel(t);
    err = X(15:17,:)' - REF(:,1:3);
    eul = zeros(3,N); for k=1:N, eul(:,k)=Q.toEulZYX(X(5:8,k)); end

    % colour convention: reference = black (dashed), actual = red (solid)
    cRef=[0 0 0]; cAct=[0.85 0.15 0.15];

    fig=figure('Name','Tracking analysis','Color','w','units','normalized','outerposition',[0 0 1 1]);
    tl=tiledlayout(fig,2,2,'Padding','compact','TileSpacing','compact');

    nexttile; hold on; grid on; box on;
    plot(t,REF(:,1),'--','Color',cRef,'LineWidth',1.2);
    plot(t,X(15,:),'-','Color',cAct,'LineWidth',1.6);
    plot(t,REF(:,2),'--','Color',cRef,'LineWidth',1.2,'HandleVisibility','off');
    plot(t,X(16,:),'-','Color',cAct,'LineWidth',1.6,'HandleVisibility','off');
    ylabel('position [m]'); title('Horizontal position (x, y)');
    legend({'command','actual'},'Location','best'); xlim([0 t(end)]);

    nexttile; hold on; grid on; box on;
    plot(t,-REF(:,3),'--','Color',cRef,'LineWidth',1.2);
    plot(t,-X(17,:),'-','Color',cAct,'LineWidth',1.6);
    ylabel('altitude [m]'); title('Altitude');
    legend({'command','actual'},'Location','best'); xlim([0 t(end)]);

    nexttile; hold on; grid on; box on;
    en=vecnorm(err,2,2);
    plot(t,en,'-','Color',cAct,'LineWidth',1.2);
    ylabel('error [m]'); xlabel('time [s]');
    title(sprintf('Position error  (RMS %.3f m,  max %.3f m)',sqrt(mean(en.^2)),max(en)));
    legend({'error norm'},'Location','best');
    xlim([0 t(end)]); ylim([0 max(0.05,max(en)*1.15)]);

    nexttile; hold on; grid on; box on;
    plot(t,rad2deg(unwrap(REF(:,4))),'--','Color',cRef,'LineWidth',1.2);
    plot(t,rad2deg(unwrap(eul(3,:))),'-','Color',cAct,'LineWidth',1.6);
    ylabel('heading [deg]'); xlabel('time [s]'); title('Heading (yaw)');
    legend({'command','actual'},'Location','best'); xlim([0 t(end)]);

    title(tl,'F450 Trajectory Tracking','FontWeight','bold');
end

function states(t, X, REF, LOG, ~)
%STATES  3x3 figure: velocity / attitude / body-rate references vs actual.
    Q = quatUtils();
    N = numel(t);

    % actual signals
    vel_act = X(12:14,:);                    % vx,vy,vz [m/s]
    eul_act = zeros(3,N); rate_act = X(9:11,:);
    for k=1:N, eul_act(:,k)=Q.toEulZYX(X(5:8,k)); end
    % commanded signals
    vel_cmd = LOG.vel_sp;                     % vx,vy,vz sp
    eul_cmd = zeros(3,N);
    for k=1:N, eul_cmd(:,k)=Q.toEulZYX(LOG.qd(:,k)); end
    rate_cmd= LOG.rate_sp;

    cRef=[0 0 0]; cAct=[0.85 0.15 0.15];
    figure('Name','States: velocity / attitude / rate','Color','w', ...
           'units','normalized','outerposition',[0 0 1 1]);
    tl=tiledlayout(3,3,'Padding','compact','TileSpacing','compact');

    vlab={'v_x [m/s]','v_y [m/s]','v_z [m/s]'};
    for i=1:3
        nexttile; hold on; grid on; box on;
        plot(t,vel_cmd(i,:),'--','Color',cRef,'LineWidth',1.2);
        plot(t,vel_act(i,:),'-','Color',cAct,'LineWidth',1.4);
        ylabel(vlab{i}); xlim([0 t(end)]);
        title(['Velocity  ' vlab{i}]);
        legend({'command','actual'},'Location','best');
    end
    alab={'roll \phi [deg]','pitch \theta [deg]','yaw \psi [deg]'};
    for i=1:3
        nexttile; hold on; grid on; box on;
        % yaw is a continuously turning angle -> unwrap; roll/pitch stay raw
        if i==3, cmd = unwrap(eul_cmd(i,:)); act = unwrap(eul_act(i,:));
        else,    cmd = eul_cmd(i,:);          act = eul_act(i,:);  end
        plot(t,rad2deg(cmd),'--','Color',cRef,'LineWidth',1.2);
        plot(t,rad2deg(act),'-','Color',cAct,'LineWidth',1.4);
        ylabel(alab{i}); xlim([0 t(end)]);
        title(['Attitude  ' alab{i}]);
        legend({'command','actual'},'Location','best');
    end
    rlab={'p [deg/s]','q [deg/s]','r [deg/s]'};
    for i=1:3
        nexttile; hold on; grid on; box on;
        plot(t,rad2deg(rate_cmd(i,:)),'--','Color',cRef,'LineWidth',1.2);
        plot(t,rad2deg(rate_act(i,:)),'-','Color',cAct,'LineWidth',1.4);
        ylabel(rlab{i}); xlabel('time [s]'); xlim([0 t(end)]);
        title(['Body rate  ' rlab{i}]);
        legend({'command','actual'},'Location','best');
    end
    title(tl,'Velocity / Attitude / Body-rate tracking','FontWeight','bold');
end

function motors(t, LOG, ~)
%MOTORS  Two-panel figure.
%   LEFT  - per-motor thrusts laid out like the drone seen from above
%           (forward = up, right = right), total thrust on top.
%   RIGHT - mixer command: collective thrust F_z and body torques
%           tau_x (roll), tau_y (pitch), tau_z (yaw).
    T    = LOG.T_mot;  Ttot = sum(T,1);
    tau  = LOG.torque; Fz   = LOG.Fz;

    % single consistent colour for every panel (motors, thrust, torques)
    cLine=[0.20 0.40 0.70];
    cMot=cLine; cTot=cLine;
    cFz =cLine; cTx =cLine; cTy =cLine; cTz =cLine;
    te=t(end);

    figure('Name','Motor thrusts & mixer commands','Color','w', ...
           'units','normalized','outerposition',[0 0 1 1]);

    % ---------- LEFT: motor thrust map ----------
    % total thrust (top, spans left half)
    axes('Position',[0.06 0.74 0.40 0.18]); hold on; grid on; box on;
    plot(t,Ttot,'-','Color',cTot,'LineWidth',1.5);
    ylabel('thrust [N]'); title('Total thrust to motors'); xlim([0 te]);
    legend({'\Sigma T'},'Location','best');

    % physical 2x2 layout (view from above):
    %   front-left (M3)  front-right (M1)
    %   rear-left  (M2)  rear-right  (M4)
    mp = { [0.06 0.41 0.18 0.22], 3, 'M3  front-left (CW)';
           [0.28 0.41 0.18 0.22], 1, 'M1  front-right (CCW)';
           [0.06 0.10 0.18 0.22], 2, 'M2  rear-left (CCW)';
           [0.28 0.10 0.18 0.22], 4, 'M4  rear-right (CW)'};
    for r=1:4
        axes('Position',mp{r,1}); hold on; grid on; box on;
        plot(t,T(mp{r,2},:),'-','Color',cMot,'LineWidth',1.2);
        title(mp{r,3}); ylabel('T [N]'); xlim([0 te]);
        legend({sprintf('M%d',mp{r,2})},'Location','best');
        if r>=3, xlabel('time [s]'); end
    end

    % ---------- RIGHT: mixer commands ----------
    cmd = { [0.56 0.74 0.40 0.18], 'Collective thrust  F_z [N]',        Fz,       cFz, 'F_z';
            [0.56 0.52 0.40 0.18], 'Roll torque  \tau_x [N\cdotm]',     tau(1,:), cTx, '\tau_x';
            [0.56 0.30 0.40 0.18], 'Pitch torque  \tau_y [N\cdotm]',    tau(2,:), cTy, '\tau_y';
            [0.56 0.08 0.40 0.18], 'Yaw torque  \tau_z [N\cdotm]',      tau(3,:), cTz, '\tau_z'};
    for r=1:4
        axes('Position',cmd{r,1}); hold on; grid on; box on;
        plot(t,cmd{r,3},'-','Color',cmd{r,4},'LineWidth',1.3);
        title(cmd{r,2}); xlim([0 te]);
        legend({cmd{r,5}},'Location','best');
        if r==4, xlabel('time [s]'); end
    end
end