function [ps_0, rs_0] = loadICs(str)
global a e
if nargin == 0
    str = 'default';
end

switch str
    case 'qs'
        x_0     = 0; %7; %20;                      % Relative Position [m]
        y_0     = -10; %10;                          % "         "
        z_0     = 0; %5;                           % "         "

        xdot_0  = 0;%-2; %2;                            % Relative Velocity [m/s]
        ydot_0  = 0;%3; %1;                             % "         "
        zdot_0  = 0;%5; %-1;                            % "         "

        q_1_0   = 0;    % Quaternion rotating chaser body frame to LVLH frame
        q_2_0   = 0;    % "         "
        q_3_0   = sqrt(2)/2;    % "         "
        q_4_0   = sqrt(2)/2;    % "         "

        w_x_0   = 0;    % Body Angular Rates [rad/s]
        w_y_0   = 0;    % "         "
        w_z_0   = 0;    % "         "
    case 'q'
        x_0     = -10; %7; %20;                      % Relative Position [m]
        y_0     = 15; %10;                          % "         "
        z_0     = 4; %5;                           % "         "

        xdot_0  = 1;%-2; %2;                            % Relative Velocity [m/s]
        ydot_0  = -1;%3; %1;                             % "         "
        zdot_0  = 2;%5; %-1;                            % "         "

        q_1_0   = 0;    % Quaternion rotating chaser body frame to LVLH frame
        q_2_0   = 0;    % "         "
        q_3_0   = 0;    % "         "
        q_4_0   = 1;    % "         "

        w_x_0   = 0;    % Body Angular Rates [rad/s]
        w_y_0   = 0;    % "         "
        w_z_0   = 0;    % "         "
    case 'sp'
        x_0     = 0; %7; %20;                      % Relative Position [m]
        y_0     = -5000; %10;                          % "         "
        z_0     = 0; %5;                           % "         "

        xdot_0  = 0;%-2; %2;                            % Relative Velocity [m/s]
        ydot_0  = 0;%3; %1;                             % "         "
        zdot_0  = 0;%5; %-1;                            % "         "

        q_1_0   = 0;    % Quaternion rotating chaser body frame to LVLH frame
        q_2_0   = 0;    % "         "
        q_3_0   = 0;    % "         "
        q_4_0   = 1;    % "         "

        w_x_0   = 0;    % Body Angular Rates [rad/s]
        w_y_0   = 0;    % "         "
        w_z_0   = 0;    % "         "
    case 'sp_nl'
        x_0     = -a*(1-e)+sqrt(a^2*(1-e)^2-5000^2); %7; %20;                      % Relative Position [m]
        y_0     = -5000; %10;                          % "         "
        z_0     = 0; %5;                           % "         "

        xdot_0  = 0;%-2; %2;                            % Relative Velocity [m/s]
        ydot_0  = 0;%3; %1;                             % "         "
        zdot_0  = 0;%5; %-1;                            % "         "

        q_1_0   = 0;    % Quaternion rotating chaser body frame to LVLH frame
        q_2_0   = 0;    % "         "
        q_3_0   = 0;    % "         "
        q_4_0   = 1;    % "         "

        w_x_0   = 0;    % Body Angular Rates [rad/s]
        w_y_0   = 0;    % "         "
        w_z_0   = 0;    % "         "
    case 'big'
        x_0     = -2108206.85; %7; %20;                      % Relative Position [m]
        y_0     = -2000000; %10;                          % "         "
        z_0     = 0; %5;                           % "         "

        xdot_0  = 0;%-2; %2;                            % Relative Velocity [m/s]
        ydot_0  = 1854.3156;%3; %1;                             % "         "
        zdot_0  = 5000;%5; %-1;                            % "         "

        q_1_0   = 0;    % Quaternion rotating chaser body frame to LVLH frame
        q_2_0   = 0;    % "         "
        q_3_0   = 0;    % "         "
        q_4_0   = 1;    % "         "

        w_x_0   = 0;    % Body Angular Rates [rad/s]
        w_y_0   = 0;    % "         "
        w_z_0   = 0;    % "         "
    case 'spo'
        x_0     = 5000; %7; %20;                      % Relative Position [m]
        y_0     = 0; %10;                          % "         "
        z_0     = 0; %5;                           % "         "

        xdot_0  = 0;%-2; %2;                            % Relative Velocity [m/s]
        ydot_0  = -2*n*x0;%3; %1;                             % "         "
        zdot_0  = 0;%5; %-1;                            % "         "

        q_1_0   = 0;    % Quaternion rotating chaser body frame to LVLH frame
        q_2_0   = 0;    % "         "
        q_3_0   = 0;    % "         "
        q_4_0   = 1;    % "         "

        w_x_0   = 0;    % Body Angular Rates [rad/s]
        w_y_0   = 0;    % "         "
        w_z_0   = 0;    % "         "
    case 1
        x_0     = -100; %7; %20;                      % Relative Position [m]
        y_0     = -110; %10;                          % "         "
        z_0     = -250; %5;                           % "         "

        xdot_0  = 0;%-2; %2;                            % Relative Velocity [m/s]
        ydot_0  = 0;%3; %1;                             % "         "
        zdot_0  = 0;%5; %-1;                            % "         "

        q_1_0   = 0;    % Quaternion rotating chaser body frame to LVLH frame
        q_2_0   = 0;    % "         "
        q_3_0   = 0;    % "         "
        q_4_0   = 1;    % "         "

        w_x_0   = 0;    % Body Angular Rates [rad/s]
        w_y_0   = 0;    % "         "
        w_z_0   = 0;    % "         "
    case 2
        x_0     = 100; %7; %20;                      % Relative Position [m]
        y_0     = 110; %10;                          % "         "
        z_0     = 250; %5;                           % "         "

        xdot_0  = -2; %2;                            % Relative Velocity [m/s]
        ydot_0  = 3; %1;                             % "         "
        zdot_0  = 5; %-1;                            % "         "

        q_1_0   = 0;    % Quaternion rotating chaser body frame to LVLH frame
        q_2_0   = 0;    % "         "
        q_3_0   = 0;    % "         "
        q_4_0   = 1;    % "         "

        w_x_0   = 0;    % Body Angular Rates [rad/s]
        w_y_0   = 0;    % "         "
        w_z_0   = 0;    % "         "
    case 3
        x_0     = 100; %7; %20;                      % Relative Position [m]
        y_0     = 0; %10;                          % "         "
        z_0     = 160; %5;                           % "         "

        xdot_0  = 0;%-2; %2;                            % Relative Velocity [m/s]
        ydot_0  = 0;%3; %1;                             % "         "
        zdot_0  = 0;%5; %-1;                            % "         "

        q_1_0   = 0;    % Quaternion rotating chaser body frame to LVLH frame
        q_2_0   = 0;    % "         "
        q_3_0   = 0;    % "         "
        q_4_0   = 1;    % "         "

        w_x_0   = 0;    % Body Angular Rates [rad/s]
        w_y_0   = 0;    % "         "
        w_z_0   = 0;    % "         "
    case 'tum'
        x_0     = 100; %7; %20;                      % Relative Position [m]
        y_0     = 110; %10;                          % "         "
        z_0     = 250; %250; %5;                           % "         "

        xdot_0  = 0.1;%-2; %2;                            % Relative Velocity [m/s]
        ydot_0  = 0.1;%3; %1;                             % "         "
        zdot_0  = 0.1;%5; %-1;                            % "         "

        q_1_0   = 0;    % Quaternion rotating chaser body frame to LVLH frame
        q_2_0   = 0;    % "         "
        q_3_0   = 0;    % "         "
        q_4_0   = 1;    % "         "

        w_x_0   = 1;    % Body Angular Rates [rad/s]
        w_y_0   = 1;    % "         "
        w_z_0   = 1;    % "         "
    otherwise
        x_0     = 100; %7; %20;                      % Relative Position [m]
        y_0     = 110; %10;                          % "         "
        z_0     = 250; %250; %5;                           % "         "

        xdot_0  = 0.1;%-2; %2;                            % Relative Velocity [m/s]
        ydot_0  = 0.1;%3; %1;                             % "         "
        zdot_0  = 0.1;%5; %-1;                            % "         "

        q_1_0   = 0;    % Quaternion rotating chaser body frame to LVLH frame
        q_2_0   = 0;    % "         "
        q_3_0   = 0;    % "         "
        q_4_0   = 1;    % "         "

        w_x_0   = 0;    % Body Angular Rates [rad/s]
        w_y_0   = 0;    % "         "
        w_z_0   = 0;    % "         "
end

ps_0 = [x_0, y_0, z_0, xdot_0, ydot_0, zdot_0]';     % Pos/Vel IC vector
q_0 = [q_1_0, q_2_0, q_3_0, q_4_0]';                 % Chaser quat IC vector
w_0 = [w_x_0, w_y_0, w_z_0]';                        % Chaser body rate IC
rs_0 = [q_0; w_0];                                   % Rotation IC vector













end