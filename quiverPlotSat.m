function quiverPlotSat(ani, varargin)

clear fig2;

global t x y z xdot ydot zdot q1 q2 q3 q4 uF ...
        Thrust_dir A_transl movieName movie_rot_rate

try
    cd('C:\Users\zkahl\Dropbox (MIT)\Thesis\Optimal GEO\GPOPS, etc');
catch it
    fprintf('Couldn''t cd to Dropbox Thesis directory.\n');
end%try

% Solve unforced Clohessy Wiltshire dynamics to determine unforced
% trajectory
syms X(T) Y(T) Z(T) XDOT(T) YDOT(T) ZDOT(T)
STATE = [X; Y; Z; XDOT; YDOT; ZDOT];
IC = STATE(0) == [x(1); y(1); z(1); xdot(1); ydot(1); zdot(1)];
eqn = diff(STATE) == A_transl*STATE;
[xSol(T), ~, ySol(T), ~, zSol(T), ~] = dsolve(eqn,IC);

unforcedTraj_x = double(xSol(t));
unforcedTraj_y = double(ySol(t));
unforcedTraj_z = double(zSol(t));

% Initialize figure and storage variables
fig2 = figure('units','normalized','outerposition',[0 0 1 1]);

arrow_x = [];
arrow_y = [];
arrow_z = [];

limits = zeros(1,6);

limits(1) = 1.1*min([x; unforcedTraj_x]);
limits(2) = 1.1*max([x; unforcedTraj_x]);
limits(3) = 1.1*min([y; unforcedTraj_y]);
limits(4) = 1.1*max([y; unforcedTraj_y]);
limits(5) = 1.1*min([z; unforcedTraj_z]);
limits(6) = 1.1*max([z; unforcedTraj_z]);

arrow_scale = 0.02*max(abs(limits));

% default view
v_loc = [-37.5 30];
if ~isempty(varargin)
    if varargin{1} == 1
        fprintf('Default 3D view (az:-37.5, el:30) set.\n');
    elseif varargin{1} == 2
        v_loc = [0,90];
        fprintf('Default 2D view (top-down) set.\n');
    elseif varargin{1} == 3
        v_loc = [0,0];
        fprintf('Trailing-the-target view set.\n');
    elseif varargin{1} == 4
        v_loc = [90,0];
        fprintf('Looking-back-at-Earth view set.\n');
    elseif varargin{1} == 'a'
        v_loc = [60,90];
        view(v_loc);
        fprintf('Az: 60, El: 60 view set.\n');
    elseif varargin == 'b'
        v_loc = [150,60];
        view(v_loc);
        fprintf('Az: 150, El: 60 view set.\n');
    elseif varargin == 'c'
        v_loc = [240,60];
        view(v_loc);
        fprintf('Az: 240, El: 60 view set.\n');
    elseif varargin == 'd'
        v_loc = [330,60];
        view(v_loc);
        fprintf('Az: 330, El: 60 view set.\n');
    elseif varargin == 'e'
        v_loc = [60,30];
        view(v_loc);
        fprintf('Az: 60, El: 30 view set.\n');
    elseif varargin == 'f'
        v_loc = [150,30];
        view(v_loc);
        fprintf('Az: 150, El: 30 view set.\n');
    end%if
end%if

if ani == 1
    vidString = strcat(movieName,'.avi');
    vidObj = VideoWriter(vidString);
    open(vidObj); 
    ax = gca; vaz = v_loc(1); vel = v_loc(2);
    for j = 1:1:numel(t)
        th_arrow_rsw = alt_quatrotate([q1(j) q2(j) q3(j) q4(j)],Thrust_dir);
        th_arrow_r(j) = arrow_scale.*th_arrow_rsw(1);
        th_arrow_s(j) = arrow_scale.*th_arrow_rsw(2);
        th_arrow_w(j) = arrow_scale.*th_arrow_rsw(3);

        v_arrow_r(j) = xdot(j);
        v_arrow_s(j) = ydot(j);
        v_arrow_w(j) = zdot(j);

        ax.View = [vaz vel];
        if vel ~= 0 && vel ~=90 % Rotate the plot if 3D view
            vaz = vaz + movie_rot_rate;
        end%if
        th = quiver3(ax, x(j),y(j),z(j),uF(j)*th_arrow_r(j),uF(j)*th_arrow_s(j),uF(j)*th_arrow_w(j), ...
            'rs','AutoScale','off', ...
            'LineWidth',1.5, ...
            'MarkerFaceColor','r',...
            'ShowArrowHead','on', ...
            'MaxHeadSize',1.5);
        hold on; % [vaz,vel] = view;
        vl = quiver3(ax, x(j),y(j),z(j),v_arrow_r(j),v_arrow_s(j),v_arrow_w(j), ...
            'bs','AutoScale','off', ...
            'LineWidth',1.5, ...
            'ShowArrowHead','on', ...
            'MaxHeadSize',1.5);
        unf_traj = plot3(unforcedTraj_x(1:j),unforcedTraj_y(1:j),unforcedTraj_z(1:j),'b--');
        traj = plot3(x(1:j),y(1:j),z(1:j),'c--');
        axis(limits);
        if j==1
            start = plot3(x(j),y(j),z(j),'go','MarkerSize',10,'MarkerFaceColor','g');
            tgt = plot(0,0,'mo','MarkerSize',10,'MarkerFaceColor','m');
            title('Optimal Rendezvous Visualization');
            legend('Thrust Vector','Velocity Vector','Unforced Trajectory','Forced Trajectory','Chaser Start','Target');
            xlabel('R [m]'); ylabel('S [m]'); zlabel('W [m]');
        end%if
        currFrame = getframe(fig2); delete(th); delete(vl);%delete(traj);
        writeVideo(vidObj,currFrame);
    end%for
    close(vidObj);
    
elseif ani == 0
    
    ax = gca; vaz = v_loc(1); vel = v_loc(2);
    for j = 1:numel(t)                
        th_arrow_rsw = alt_quatrotate([q1(j) q2(j) q3(j) q4(j)],Thrust_dir);

        th_arrow_r(j) = arrow_scale.*th_arrow_rsw(1);
        th_arrow_s(j) = arrow_scale.*th_arrow_rsw(2);
        th_arrow_w(j) = arrow_scale.*th_arrow_rsw(3);

        v_arrow_r(j) = xdot(j);
        v_arrow_s(j) = ydot(j);
        v_arrow_w(j) = zdot(j);

        ax.View = [vaz vel];
        hold on; grid on;
        quiver3(x(j),y(j),z(j),uF(j)*th_arrow_r(j),uF(j)*th_arrow_s(j),uF(j)*th_arrow_w(j), ...
            'rs','AutoScale','off', ...
            'LineWidth',1.5, ...
            'MarkerFaceColor','r',...
            'ShowArrowHead','on', ...
            'MaxHeadSize',1.5);
        quiver3(x(j),y(j),z(j),v_arrow_r(j),v_arrow_s(j),v_arrow_w(j), ...
            'bs','AutoScale','off', ...
            'LineWidth',1.5, ...
            'ShowArrowHead','on', ...
            'MaxHeadSize',1.5);
        unf_traj = plot3(unforcedTraj_x(1:j),unforcedTraj_y(1:j),unforcedTraj_z(1:j),'b--');
        traj = plot3(x(1:j),y(1:j),z(1:j),'c--');
        axis(limits);
        if j==1 || j == numel(t)
            start = plot3(x(j),y(j),z(j),'go','MarkerSize',10,'MarkerFaceColor','g');
            tgt = plot(0,0,'mo','MarkerSize',10,'MarkerFaceColor','m');
            title('Optimal Rendezvous Visualization');
            legend('Thrust Vector','Velocity Vector','Unforced Trajectory','Forced Trajectory','Chaser Start','Target');
            xlabel('R [m]'); ylabel('S [m]'); zlabel('W [m]');
        end%if
    end%for
    
end%if








end%function