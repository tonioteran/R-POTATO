% R_POTATO_Main
%% Script for GEO Maneuver Optimization -- Zack's Thesis 
%  2Lt Zachary K. Funke
%  Rev 0.9  --  9 Feb 2017
%
%  Optimal Rendezvous of an active chaser satellite with a passive target
%  In circular orbit

clc; %close all;
clearvars -except s;

try
    cd('C:\Users\Zack Funke\Dropbox (MIT)\Thesis\Optimal GEO\GPOPS, etc');
    fprintf('*\n* Path automatically set to Thesis Optimal GEO folder (Zack''s Cadet Laptop).\n*\n');
catch it % don't throw an error if this isn't Zack's Computer
    try 
        cd('C:\Users\zkahl\Dropbox (MIT)\Thesis\Optimal GEO\GPOPS, etc');
        fprintf('*\n* Path automatically set to Thesis Optimal GEO folder (Zack''s New Laptop).\n*\n');
    catch it
        fprintf('Note: path could not be automatically set, or this isn''t Zack''s Computer. \nIf you encounter an error, you may need to make sure all .m files are called from the same location.\n\n');
    end%try
end%try

global a e t n s x y z xdot ydot zdot q1 q2 q3 q4 wx wy wz u uF uFx uFy uFz uTx uTy uTz p ...
    Thrust_dir A_transl movieName movie_rot_rate


               % ----- % User-Defined Options % ----- %
               
                             % Optimization Flags and Parameters %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                  Incl_Att_Dynamics = 0;             %  <----- Set to 1 to include quaternions and body rates for chaser.
                %------------------------------------%----------------------------------------------------------------------------
                    Min_Time        = 0;             %  <----- Objective is to minimize time of maneuver.
                        Sched_Pol   = 0;             %  <----- Use 'Scheduling Policy' approach to Min Time problem (regulate, then wait).
                        Hybrid_Opt  = 0;             %  <----- Use Min Time control for limiting coordinate, and Min Energy for others. 
                %------------------------------------%--------------------------------------------------------------------------------------
                    Min_Fuel        = 0;             %  <----- Objective is to minimize |u|.
                        One_Phase   = 0;             %  <----- Attempt to find solution on one phase
                        Three_Phase = 1;             %  <----- Bang - Off - Bang
                        Bang_Width  = 10;            %  <----- Max width of min-fuel control period [s]
                        Decoup_Bang = 0;             %  <----- Decouples bang on-times for in-plane vs out-of-plane coordinates 
                %------------------------------------%-----------------------------------------------------------------------------------------
                    Quadratic_Cost  = 0;             %  <----- LQR cost function of the form 1/2*(x^T*Q*x + u^T*R*u)
                        Q_R_ratio   = 1e-10;         %  <----- This is the desired relative weight of state regulation vs. control energy expenditure
                %------------------------------------%-----------------------------------------------------------------------------------------------    
                    Min_Energy      = 1;             %  <----- Objective is to minimize 0.5*u^2.
                %------------------------------------%-------------------------------------------------------------------------------------------------  
                    Target_Point    = 1;             %  <----- Selecting this will cause solver to regulate state to the origin (1 phase, unless Min Fuel)
                    Avoid_FOV       = 0;             %  <----- Enforce no-fly zone in nadir FOV of tgt (only for Target_Point)
                          FOV       = 60;            %  <----- FOV in degrees
                          FOV_buff  = 30;            %  <----- Desired 'buffer' cone around no-fly zone
                                                     %
                    Target_Ellipse  = 0;             %  <----- Select this to add a second phase for ellipse.
                        alpha       = 100;           %  <----- Design the semi-major axis of terminal ellipse (CW only) [m]
                        zeta        = 60;            %  <----- 'inclination' of free-orbit ellipse (CW only) [deg]
                        gamma       = 90;            %  <----- angle between free-orbit line of nodes and the -X axis, measured clockwise (dir. of relative motion) (CW only) [deg]
                        Vbar_Offset = 0;             %  <----- Design the offset of the center of the ellipse with respect to the LVLH origin [m]
                    coast           = 1;             %  <----- If set, this prohibits thruster firing during the ellipse phase (using above parameters for endoint constraint.) 
                        range_path  = 1;             %  <----- If coast = 0, enforce a constant range path constraint for second phase.
                            range   = 500;           %  <----- Value of range for the constant range path constraint. [m]
                        period_path = 0;             %  <----- If coast = 0, enforce a period-matching path constraint for second phase.
                        center_path = 1;             %  <----- If coast = 0, enforce a path constraint on the center of the second phase's periodic trajectory.
                            center  = 0;             %  <----- Value of the center enforced by the centering path constraint. [m]
                                                     %
                    Period_Match    = 0;             %  <----- If this flag is set, the chaser's goal is ONLY to null secular relative dynamics. 
                        ThrustTime  = 20;            %  <----- If the objective is minimum-time with attitude dynamics, specify the max duration of the thrusting phase
                %------------------------------------%----------------------------------------------------------------------------------------------------------------------------------
                    Nonlinear       = 0;             % <----- Instead of HCW eqns, use full nonlinear Keplerian equations
                        e           = 0.7;           % <----- Target / reference orbit eccentricity.
                        num_nu_pts  = 1e5;           % <----- Number of points to be calculated for the true anomaly lookup table
                        de_ON       = 0;             % <--- Set this flag to constrain the size of the periodic trajectory. WARNING: increases computation time!
                            x_size  = 50;          % <----- Desired maximum radial offset of periodic trajectory (if target-centered). For positive de = e_ch - e_tgt, specify a positive value, & vice versa
                        dw_ON       = 1;             % <--- Set this flag to constrain the center of the periodic trajectory. WARNING: increased computation time!
                            dw      = 0;             % <----- Desired differential argument of perigee (analogous to Vbar_Offset). CAUTION: highly sensitve, use small values
                        di_ON       = 0;             % <--- Set this flag to constrain the inclination of the periodic trajectory. WARNING: increases computation time!
                            di      = 0;             % <----- Desired differential inclination (analagous to zeta). CAUTION: highly sensitive, use small values
                        dG_ON       = 0;             % <--- Set this flag to constrain the azimuth of the local ascending node of the periodic trajectory. WARNING: increases computation time!
                            dG      = 0;             % <----- Desired differential ascending node azimuth (analagous to gamma). CAUTION: highly senstive, use small values                              
                                                     % 
                    enforce_CW_ellipse = 0;          % <----- IF coast = 0, set path constraint to enforce tgt-centered CW ellipse for 1 phase
                        CW_alpha    = 10000;           % <--- Set the semimajor axis of the CW ellipse to enforce (Your IC will be overwritten to start along this ellipse)
                        CW_zeta     = 60;           % <--- Set the inclination (about Y axis) of the CW ellipse to enforce                            
                                                    %
                    
                    Max_Time        = 1.0*86164;         %  1 sidereal day (for Min_Energy only)
                    Mnvr_Epoch      = 0.*Max_Time;     % Time since perigee (Nonlinear only) 
%                     Max_Time        = 12*3600;          %  1 hour (for Min_Energy only)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            % Solver Parameters *                 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    MeshTolerance   = 1e-4;         % Solver mesh tolerance (Try 1e-6 for rendezvous or 1e-3 for attitude dynamics or ellipse phasing)
                    AutomaticBounds = 0;            % Set this to have GPOPS-II scale the problem based on the given bounds (good place to start, but sometimes finds local minima)
                    MaxIPOPTiters   = 2000;         % Set the maximum number of IPOPT iterations per mesh refinement (default is 2000, but try 1000 for harder problems)
                    MaxMeshIters    = 10;           % Set the maximum number of GPOPS-II adaptive mesh refinements (start low if you're not sure problem will work)
                    ColPointsMin    = 4;           % Set the minimum number of collocation points. Try 10 for most problems, but use 4 for Nonlinear && Target_Ellipse.
                    ColPointsMax    = 10;           % Set the maximum number of collocatio
                    MaxRange        = 20000;        % LVLH coordinate maximum for problem [m] (ensure the solution isn't artificially bounded!)
                    MaxVel          = 100;          % LVLH coordinate velocity max for problem [m/s] (100 m/s is very conservative)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            % Physical Properties *
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    HallEffectThrust  = 1;           %  <----- Max thrust is 0.5 N
                    HydrazineMonoProp = 0;           %  <----- Max thrust is 22 N
                    LEO = 0;                         %  Target semimajor axis: 6,378,537 m
                    GEO = 1;                         %  Target semimajor axis: 42,164,137 m
                    m       = 100;                   %  Chaser mass [kg]
                                                     %
                    I       = [17   0    0  ; ...    %  Chaser moment of inertia [kg-m^2]
                               0   17    0  ; ...    % (Based on 1m x 1m x 1m uniform 100kg cube)
                               0    0   17 ];        %  Note: Continuous needs to be modified if 
                                                     %  non-principal inertia tensor is used.
                                                     %
                    Tmax    = 1;                     %  Maximum torque in each axis [N-m]
                                                     %
                    Thrust_dir = [1; 0; 0];          %  In which body-axis direction does the 
                                                     %  thruster produce force?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            % Initial Conditions %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
             Load_ICs    = 0;                        % If 1, load ICs from loadICs.m and ignore below
                load_str = '';                       % This string is the argument for loadICs.m (can be empty)
 Load_Prev_Finalstate    = 0;                        % If 1, load the finalstate from the workspace of the previous run and ignore below      
                                                     %
                 x_0     = 100;               % Relative Position [m]
                 y_0     = -110;                 % "         "
                 z_0     = -250;                     % "         "
                                                     %
                 xdot_0  = 0.1;%-2; %2;          % Relative Velocity [m/s]
                 ydot_0  = 0.1;%3; %1;                 % "         "
                 zdot_0  = 0.2;%5; %-1;                % "         "
                                                     %
                 q_1_0   = 0;                        % Quaternion resolving chaser body frame vector into Hill's frame
                 q_2_0   = 0;                        % "         "
                 q_3_0   = 0;                        % "         "
                 q_4_0   = 1;                        % "         "
                                                     %
                 w_x_0   = 0;                        % Body Angular Rates [rad/s]
                 w_y_0   = 0;                        % "         "
                 w_z_0   = 0;                        % "         "
                                                     %
                                                     % Note: Position ICs
                                                     % will be overwritten if Nonlinear and enforce_CW_ellipse specified
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            % Data Visualization %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
                    Plot3D_Trajectory = 1;            %  <----- Plot 3D solution trajectory
                                                      %
                    MakeAnimation     = 0;            %  <----- If 1, this generates an .avi movie of the trajectory located in the current directory.
                        View          = 1;            %  <--- Default 3D view with rotation
%                         View         = 2;            % <-- top-down 2D view
%                         View         = 3;            % <-- trailing-target 2D view
%                         View         = 4;            % <-- look-back-at-Earth 2D view
                   movie_rot_rate  = 2;                % <-- degrees of movie rotation per frame (3D views only)
                   movieName = 'quiverPlotSat_pres2';  % 
                                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            
                            
if sum([Target_Point Target_Ellipse Period_Match (Nonlinear && enforce_CW_ellipse)]) ~= 1 
    error('Need to specify ONE type of maneuver desired. Target_Point: State regulation (rendezvous). Target_Ellipse: Transfer to specified circumnavigation ellipse. Period_Match: Null secular dynamics.')
end%if
if sum([Min_Time Min_Energy Min_Fuel Quadratic_Cost]) ~= 1
    error('Specify ONE objective for optimization. (Min_Time, Min_Energy, Min_Fuel, Quadratic_Cost)')
end%if

% --------------------- Define constants ----------------------------------
Re = 6378.137 * 1000;           % Earth's mean equatorial radius [m]
mu = 3.986004419e14;            % Earth's gravitational parameter [m^3/s^2]

% --------------------- Problem-specific parameters -----------------------
if LEO
    h   = 400 * 1000;           % LEO altitude [m]
elseif GEO
    h   = 35786 * 1000;         % GEO altitude [m]
end%if

if HallEffectThrust
    Fmax    = 0.5;              % Maximum HET thrust in each axis [N]
elseif HydrazineMonoProp
    Fmax    = 22;               % Maximum MPT thrust in each axis [N]
end%if

% --------------------- Calculated parameters -----------------------------
a  = Re + h;                    % Calculate semi-major axis
n  = sqrt(mu/(a^3));            % Calculate constant orbital rate
b  = Fmax/m;                    % This value denotes the maximum control 
                                % input acceleration  
                                
% --------------------- PACK INITIAL CONDITIONS FOR STATE -----------------

if Load_ICs
    [ps_0,rs_0] = loadICs(load_str);
elseif Load_Prev_Finalstate
    ps_0 = s(end,1:6)';
    q_0 = [q_1_0, q_2_0, q_3_0, q_4_0]';        % Chaser quat IC vector
    w_0 = [w_x_0, w_y_0, w_z_0]';               % Chaser body rate IC
    rs_0 = [q_0; w_0];
    if Incl_Att_Dynamics
        try
        rs_0 = s(end,7:13);
        catch norotationhistory
            warning('Previous state did not include attitude, using user-specified values');
        end%trycatch
    end%if
elseif Nonlinear && enforce_CW_ellipse 
    if CW_zeta == 0
        ps_0 = [0, CW_alpha, 0, 0.5*n*CW_alpha, 0, 0]';
    else
        ps_0 = [0, CW_alpha, 0, 0.5*n*CW_alpha, 0, n*CW_alpha/2*tand(CW_zeta)]';
    end%if
    q_0 = [q_1_0, q_2_0, q_3_0, q_4_0]';
    w_0 = [w_x_0, w_y_0, w_z_0]';   
    rs_0 = [q_0; w_0];  
else% Use user-specified initial conditions
    ps_0 = [x_0, y_0, z_0, xdot_0, ydot_0, zdot_0]';     % Pos/Vel IC vector
    q_0 = [q_1_0, q_2_0, q_3_0, q_4_0]';        % Chaser quat IC vector
    w_0 = [w_x_0, w_y_0, w_z_0]';                       % Chaser body rate IC
    rs_0 = [q_0; w_0];                          % Rotation IC vector
end%if 
                                                

if e > 0 
    if Nonlinear
        e_str = num2str(e); num_str = num2str(num_nu_pts);
        nu_str = strcat('nu_table_e_',e_str,'_num_',num_str,'.mat');
        fprintf('\n\n* Non-zero eccentricity specified. Attempting to load existing lookup table. Searching directory for .mat file... \n')
        try 
            load(nu_str);
            fprintf('\n* Success! True anomaly lookup table already exists for e = %s and num_nu_pts = %s. \n\n* Proceeding to optimization loop... \n', e_str, num_str);
        catch 
            fprintf('\n* Failed to find existing lookup table. Generating true anomaly look-up table with %g entries...',num_nu_pts)
            if mod(num_nu_pts,2) > 0
                num_nu_pts = num_nu_pts + 1; % Use an even number of table entries for best results
            end%if
            nu_table = generateTrueAnomTable(n,e,num_nu_pts);
            fprintf('\n* Done! Saving table for future use.\n');
            save(nu_str,'nu_table','-mat');
        end%trycatch
    else
        e = 0;
        warning('Must specify ''Nonlinear'' if you want the target to have nonzero eccentricity. Value for e was set to zero because ''Nonlinear'' was not specified.');
    end%if
end%if



% % --------------------- Pos/Vel Dynamics (Clohessy-Wiltshire) -------------
% 
% 
%     A_transl = [ 0        0      0        1      0       0; ...
%                  0        0      0        0      1       0; ...
%                  0        0      0        0      0       1; ...
%                  3*n^2    0      0        0     2*n      0; ...
%                  0        0      0      -2*n     0       0; ...
%                  0        0    -n^2       0      0       0];
%              
%     B_transl =  b * R_0 * Thrust_dir;
%    
%     
%     
% % --------------------- Rotation Dynamics (Quaternions) -------------
% 
%     A_rot = [Q_0, zeros(4,3)  ; ...
%              zeros(3,7)      ];
%     
%     B_rot = Tmax*inv(I);
%     
% % --------------------- Combined Dynamics ---------------------------
% 
%     A_0 = [A_transl, zeros(6,7); ...
%            zeros(7,6), A_rot  ];
%      
%     B_0 = [zeros(3,4)          ; ...   
%            B_transl, zeros(3,3); ...
%            zeros(4,4)          ; ...
%            zeros(3,1), B_rot  ];
%      
%      % Note: These dynamics are not LTI. These are the A and B matrices at
%      % time t=0.

%%%%%%%%%%%%%%%%%%%%% End User-Modifiable Section %%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%% Setup Calls to Optimization Loop %%%%%%%%%%%%%%%%%%%%%%%%%%%%    

    
if Min_Time && Target_Point && ~Incl_Att_Dynamics %------------------- Min Time 3-DOF (Scheduling Policy or Hybrid Optimal) ------------------------
    time_constr = 0;
    num_phases = 1;
    
    if  ~Sched_Pol && ~Hybrid_Opt
        % Attempt to solve the min-time problem without a scheduling policy
        % or hybrid optimization. Will result in singular ambiguity for
        % 3-DOF problem.
        warning('Min-Time specified for 3-DOF problem, but Sched_Pol and Hybrid_Opt are off. Result may feature scheduling ambiguity (singular class of solutions).')
        InPlane_Only = 0; OutOfPlane_Only = 0;
        run('R_POTATO_Main');
        t = output.result.interpsolution.phase.time;
        s = output.result.interpsolution.phase.state;
        u = output.result.interpsolution.phase.control;
        p = output.result.interpsolution.phase.costate;
            x = s(:,1);
            y = s(:,2);
            z = s(:,3);
            xdot = s(:,4);
            ydot = s(:,5);
            zdot = s(:,6);
        uFx = u(:,1);
        uFy = u(:,2);
        uFz = u(:,3);

        % Calculate the Hamiltonian at each time step
        H = [];
        for i = 1:numel(t)
            xd                 = xdot(i);
            yd                 = ydot(i);
            zd                 = zdot(i);
            % Clohessy-Wiltshire Equations of Motion (with forcing term)
            xddot              = 2*n.*ydot(i) + 3*n^2.*x(i) + b.*uFx(i);
            yddot              = -2*n.*xdot(i)              + b.*uFy(i); 
            zddot              = -n^2.*z(i)                 + b.*uFz(i);
            s_dot = [xd; yd; zd; xddot; yddot; zddot];
            H(i) = 1 + p(i,:)*s_dot;
        end%for

       % Plot translation state histories (x,y,z, xd,yd,zd)
        figure();
        ax1 = subplot(1,3,1);
        plot(t,s(:,1:6)); hold on;
        plot(t, zeros(size(t)), '--k');
        if num_phases > 1
            hold on
            plot([t1(end), t1(end)], [-0.5*max(max(s)),0.5*max(max(s))],'--k')
            plot([t2(end), t2(end)], [-0.5*max(max(s)),0.5*max(max(s))],'--k')
            legend('X','Y','Z', 'X dot', 'Y dot', 'Z dot', 'Phase Separation Marker');
        else
            legend('X','Y','Z', 'X dot', 'Y dot', 'Z dot');
        end%if
        grid on;
        xlim([0 t(end)]);
        title('Translation State Histories')
        xlabel('Time (s)');
        ylabel('Translation State (m, m/s)');

        % Plot thruster control history (uFx, uFy, uFz)
        ax2 = subplot(1,3,2);
        grid on; hold on;
        line1 = plot(t,  Fmax*ones(size(u(:,1))), 'k', 'LineWidth',2);
        line2 = plot(t, -Fmax*ones(size(u(:,1))), 'k', 'LineWidth',2);
        plot(t, zeros(size(t)), '--k');
        force_x = plot(t, Fmax.*uFx);
        force_y = plot(t, Fmax.*uFy);
        force_z = plot(t, Fmax.*uFz);
        if num_phases > 1
            ph_mk = plot([t1(end), t1(end)], [-0.5*Fmax,0.5*Fmax],'--k');
            ph_mk = plot([t2(end), t2(end)], [-0.5*Fmax,0.5*Fmax],'--k');
            legend([line1 force_x force_y force_z ph_mk], 'Max/Min Thrust Limits','Radial Thrust','Transverse Thrust','Cross-Track Thrust', 'Phase Separation Marker');
        else 
            legend([line1 force_x force_y force_z], 'Max/Min Thrust Limits','Radial Thrust','Transverse Thrust','Cross-Track Thrust');
        end%if
        title('Thruster Control History')
        ylim([-1.1*Fmax 1.1*Fmax])
        xlim([0 t(end)]);
        xlabel('Time (s)');
        ylabel('Control Thrust (N)');

        % Plot Hamilton w.r.t. time
        ax3 = subplot(1,3,3);
        plot(t, H); hold on;
        plot(t, zeros(size(t)), '--k');
        grid on;
        title('Hamiltonian')
        xlabel('Time');
        ylabel('Hamiltonian');
        xlim([0 t(end)]);
        % Make it easier to zoom in synchronously on plots 
        linkaxes([ax1, ax2, ax3],'x');
   
        % Call the function that draws/animates trajectory
        % quiverPlotSat(MakeAnimation);
        return; % Don't execute the rest of the cases.
    end%if ~Sched_Pol && ~Hybrid_Opt
    
    % Otherwise:
    % Find time-optimal control for XY dynamics
    InPlane_Only = 1; OutOfPlane_Only = 0;
    run('R_POTATO_Main');

    t_xy = output.result.interpsolution.phase.time;
    s_xy = output.result.interpsolution.phase.state;
        x = s_xy(:,1);
        y = s_xy(:,2);
        xdot = s_xy(:,3);
        ydot = s_xy(:,4);
    u_xy = output.result.interpsolution.phase.control;
        uFx = u_xy(:,1);
        uFy = u_xy(:,2);
    p_xy = output.result.interpsolution.phase.costate;

    % Find time-optimal control for Z dynamics
    InPlane_Only = 0; OutOfPlane_Only = 1;
    run('R_POTATO_Main');

    t_z = output.result.interpsolution.phase.time;
    s_z = output.result.interpsolution.phase.state;
        z = s_z(:,1);
        zdot = s_z(:,2);
    u_z = output.result.interpsolution.phase.control;
        uFz = u_z(:,1);
    p_z = output.result.interpsolution.phase.costate;

    % Which set of controls took the longest to regulate?
    [time_constr, I] = max([t_xy(end) t_z(end)]);
    if I == 1
        t_max = t_xy;
    elseif I == 2
        t_max = t_z;
    end%if
    
    if Hybrid_Opt && ~Sched_Pol

        % Re-run the faster dynamics with fixed final time
        if I == 1
            Min_Time = 0;
            Min_Energy = 1;
            InPlane_Only = 0; OutOfPlane_Only = 1;
            run('R_POTATO_Main');

            t_z = output.result.interpsolution.phase.time;
            s_z = output.result.interpsolution.phase.state;
                z = s_z(:,1);
                zdot = s_z(:,2);
            u_z = output.result.interpsolution.phase.control;
                uFz = u_z(:,1);
            p_z = output.result.interpsolution.phase.costate;
   
            % Calculate the Hamiltonian for Z (min energy)
            H_Z = [];
            for i = 1:numel(t_z)
                zddot              = -n^2.*z(i)  + b.*uFz(i);
                H_Z(i) = 0.5*u_z(i,:)^2 + p_z(i,:)*[zdot(i); zddot];
            end%for

            % Calculate the Hamiltonian for XY (min time)
            H_XY = [];
%             if X_Actuation;
            for i = 1:numel(t_xy)         % Case: X,Y,Z actuation
                xddot              = 2*n.*ydot(i) + 3*n^2.*x(i) + b.*uFx(i);
                yddot              = -2*n.*xdot(i)              + b.*uFy(i);
                H_XY(i) = 1 + p_xy(i,:)*[xdot(i); ydot(i); xddot; yddot];
            end%for
%             else
%                 for i = 1:numel(t_xy)
%                 s_dot = A_XY*s_xy(i,:)' + b.*B_XY(:,2)*u_xy(i,2)';      % Case: Y,Z actuation
%                 H_XY(i) = 1 + p_xy(i,:)*s_dot;
%                 end%for
%             end%if

            % Reset flags for future use
            Min_Energy = 0;
            Min_Time = 1;

        elseif I == 2
            Min_Time = 0;
            Min_Energy = 1;
            InPlane_Only = 1; OutOfPlane_Only = 0;
            run('R_POTATO_Main');

            t_xy = output.result.interpsolution.phase.time;
            s_xy = output.result.interpsolution.phase.state;
                x = s_xy(:,1);
                y = s_xy(:,2);
                xdot = s_xy(:,3);
                ydot = s_xy(:,4);
            u_xy = output.result.interpsolution.phase.control;
                uFx = u_xy(:,1);
                uFy = u_xy(:,2);
            p_xy = output.result.interpsolution.phase.costate;
            
            % Calculate the Hamiltonian for XY (min energy)
            H_XY = [];
%             if X_Actuation
            for i = 1:numel(t_xy)               % Case: X,Y,Z actuation
                xddot              = 2*n.*ydot(i) + 3*n^2.*x(i) + b.*uFx(i);
                yddot              = -2*n.*xdot(i)              + b.*uFy(i); 
                H_XY(i) = 0.5*uFx(i)^2 + 0.5*uFy(i)^2 + p_xy(i,:)*[xdot(i); ydot(i); xddot; yddot];
            end%for
%             else
%                 for i = 1:numel(t_xy)
%                     s_dot = A_XY*s_xy(i,:)' + b.*B_XY(:,2)*u_xy(i,2)';                  % Case: Y,Z actuation
%                     H_XY(i) = 0.5*u_xy(i,1)^2 + 0.5*u_xy(i,2)^2 + p_xy(i,:)*s_dot;
%                 end%for
%             end%if

            % Calculate the Hamiltonian for Z (min time)
            H_Z = [];
            for i = 1:numel(t_z)
                zddot              = -n^2.*z(i)  + b.*uFz(i);
                H_Z(i) = 1 + p_z(i,:)*[zdot(i); zddot];
            end%for

            % Reset flags for future use
            Min_Energy = 0;
            Min_Time = 1;
            
        end%if I == 1 or I == 2
        
    elseif Sched_Pol && ~Hybrid_Opt
        % Compute Min-Time Hamiltonian for all coordinates

        % Calculate the Hamiltonian for Z (min time)
        H_Z = [];
        for i = 1:numel(t_z)
            zddot              = -n^2.*z(i)  + b.*uFz(i);
            H_Z(i) = 1 + p_z(i,:)*[zdot(i); zddot];
        end%for
        
        % Calculate the Hamiltonian for XY (min time)
        H_XY = [];
%             if X_Actuation;
        for i = 1:numel(t_xy)         % Case: X,Y,Z actuation
            xddot              = 2*n.*ydot(i) + 3*n^2.*x(i) + b.*uFx(i);
            yddot              = -2*n.*xdot(i)              + b.*uFy(i);
            H_XY(i) = 1 + p_xy(i,:)*[xdot(i); ydot(i); xddot; yddot];
        end%for    
    elseif  (Sched_Pol && Hybrid_Opt)
        error('Error: Problem with initialization flags (Sched_Pol and/or Hybrid_Opt). If using the Min Time formulation, please specify whether to use a ''regulate, then wait'' scheduling policy or Min Time / Min Fuel hybrid optimization.')
    end%if for Hybrid Optimization, elseif Scheduling Policy
    
elseif Min_Time && (Target_Ellipse || Period_Match) && ~Incl_Att_Dynamics %-------- 3-DOF Min Time to Ellipse ----------
    time_constr = 0;
    num_phases = 2;
    InPlane_Only = 0; OutOfPlane_Only = 0;
    run('R_POTATO_Main');
    
elseif Min_Time && Incl_Att_Dynamics %------- Min Time with Attitude ------
    time_constr = 0;
    if Target_Point
        num_phases = 1;
    elseif Target_Ellipse %|| Period_Match 
        num_phases = 2;
    elseif Period_Match
        num_phases = 2;
    end%if
    InPlane_Only = 0; OutOfPlane_Only = 0;
    run('R_POTATO_Main');
    
elseif Min_Energy %------------------ Min Energy ---------------------
    if Target_Ellipse || Period_Match
        num_phases = 2;
    elseif Target_Point || (Nonlinear && enforce_CW_ellipse)
        num_phases = 1;
    end%if
    time_constr = Max_Time;    
    % Find energy-optimal control 
    InPlane_Only = 0; OutOfPlane_Only = 0;
    run('R_POTATO_Main');

elseif Min_Fuel %------------------ Min Fuel -------------------------

    time_constr = Max_Time;   
    % Find fuel-optimal control 
    if One_Phase
        num_phases = 1;
    else % if Three_Phase or not specified
        num_phases = 3;
    end%if
    if Decoup_Bang && Three_Phase
        
        % Run only XY dynamics:
        InPlane_Only = 1; OutOfPlane_Only = 0;
        run('R_POTATO_Main');
        t1_xy = output.result.interpsolution.phase(1).time;
        t2_xy = output.result.interpsolution.phase(2).time;
        t3_xy = output.result.interpsolution.phase(3).time;
        t_xy = [t1_xy; t2_xy(2:end,:); t3_xy(2:end,:)];
        s1_xy = output.result.interpsolution.phase(1).state;
        s2_xy = output.result.interpsolution.phase(2).state;
        s3_xy = output.result.interpsolution.phase(3).state;
        s_xy = [s1_xy; s2_xy(2:end,:); s3_xy(2:end,:)];
        x = s_xy(:,1);
        y = s_xy(:,2);
        xdot = s_xy(:,3);
        ydot = s_xy(:,4);
        u1_xy = output.result.interpsolution.phase(1).control;
        u2_xy = output.result.interpsolution.phase(2).control;
        u3_xy = output.result.interpsolution.phase(3).control;
        u_xy = [u1_xy; u2_xy(2:end,:); u3_xy(2:end,:)];
        uFx = u_xy(:,1);
        uFy = u_xy(:,2);
        p1_xy = output.result.interpsolution.phase(1).costate;
        p2_xy = output.result.interpsolution.phase(2).costate;
        p3_xy = output.result.interpsolution.phase(3).costate;
        p_xy = [p1_xy; p2_xy(2:end,:); p3_xy(2:end,:)];
        
        % Run only Z dynamics:
        InPlane_Only = 0; OutOfPlane_Only = 1;
        run('R_POTATO_Main');
        t1_z = output.result.interpsolution.phase(1).time;
        t2_z = output.result.interpsolution.phase(2).time;
        t3_z = output.result.interpsolution.phase(3).time;
        t_z = [t1_z; t2_z(2:end,:); t3_z(2:end,:)];
        s1_z = output.result.interpsolution.phase(1).state;
        s2_z = output.result.interpsolution.phase(2).state;
        s3_z = output.result.interpsolution.phase(3).state;
        s_z = [s1_z; s2_z(2:end,:); s3_z(2:end,:)];
        z = s_z(:,1);
        zdot = s_z(:,2);
        u1_z = output.result.interpsolution.phase(1).control;
        u2_z = output.result.interpsolution.phase(2).control;
        u3_z = output.result.interpsolution.phase(3).control;
        u_z = [u1_z; u2_z(2:end,:); u3_z(2:end,:)];
        uFz = u_z(:,1);
        p1_z = output.result.interpsolution.phase(1).costate;
        p2_z = output.result.interpsolution.phase(2).costate;
        p3_z = output.result.interpsolution.phase(3).costate;
        p_z = [p1_z; p2_z(2:end,:); p3_z(2:end,:)];
        
        
        % Compute Min-Fuel Hamiltonian for all coordinates

        % Calculate the Hamiltonian for Z (min fuel)
        H_Z = [];
        for i = 1:numel(t_z)
            zddot              = -n^2.*z(i)  + b.*uFz(i);
            H_Z(i) = abs(uFz(i)) + p_z(i,:)*[zdot(i); zddot];
        end%for
        
        % Calculate the Hamiltonian for XY (min time)
        H_XY = [];
        for i = 1:numel(t_xy)   
            xddot              = 2*n.*ydot(i) + 3*n^2.*x(i) + b.*uFx(i);
            yddot              = -2*n.*xdot(i)              + b.*uFy(i);
            H_XY(i) = abs(uFx(i)) + abs(uFy(i)) + p_xy(i,:)*[xdot(i); ydot(i); xddot; yddot];
        end%for    
        
    else
        InPlane_Only = 0; OutOfPlane_Only = 0;
        run('R_POTATO_Main');
    end%if
    
elseif Quadratic_Cost %------------------ Quadratic Cost -------------------------
    if Target_Ellipse || Period_Match
        num_phases = 2;
    elseif Target_Point
        num_phases = 1;
    end%if
    time_constr = 0;    % <------------- This is set to 24 hrs - CHANGE ME if there is a different tf desired!
    % Find optimal control 
    InPlane_Only = 0; OutOfPlane_Only = 0;
    run('R_POTATO_Main');

end%if

%% Solution Unpacking and Data Visualization

% Gather the solution from the GPOPS output
if num_phases > 1 && Min_Fuel
    switch Decoup_Bang
        case 0 % X,Y,Z dynamics all have same bang on-time
            t1 = output.result.interpsolution.phase(1).time;
            t2 = output.result.interpsolution.phase(2).time;
            t3 = output.result.interpsolution.phase(3).time;
            t = [t1; t2(2:end,:); t3(2:end,:)];
            s1 = output.result.interpsolution.phase(1).state;
            s2 = output.result.interpsolution.phase(2).state;
            s3 = output.result.interpsolution.phase(3).state;
            s = [s1; s2(2:end,:); s3(2:end,:)];
            u1 = output.result.interpsolution.phase(1).control;
            u2 = output.result.interpsolution.phase(2).control;
            u3 = output.result.interpsolution.phase(3).control;
            u = [u1; u2(2:end,:); u3(2:end,:)];
            p1 = output.result.interpsolution.phase(1).costate;
            p2 = output.result.interpsolution.phase(2).costate;
            p3 = output.result.interpsolution.phase(3).costate;
            p = [p1; p2(2:end,:); p3(2:end,:)];
        case 1 % Already unpacked values above.
    end%switchcase
elseif ~Incl_Att_Dynamics && (Target_Ellipse || Period_Match)
    t1 = output.result.interpsolution.phase(1).time;
    t2 = output.result.interpsolution.phase(2).time;
    t = [t1; t2];
    s1 = output.result.interpsolution.phase(1).state;
    s2 = output.result.interpsolution.phase(2).state;
    s = [s1; s2];
    u1 = output.result.interpsolution.phase(1).control;
    u2 = output.result.interpsolution.phase(2).control;
    u = [u1; u2];
    p1 = output.result.interpsolution.phase(1).costate;
    p2 = output.result.interpsolution.phase(2).costate;
    p = [p1; p2];
elseif Incl_Att_Dynamics && Period_Match
    t1 = output.result.interpsolution.phase(1).time;
    t2 = output.result.interpsolution.phase(2).time;
    t = [t1; t2];
    s1 = output.result.interpsolution.phase(1).state;
    s2 = output.result.interpsolution.phase(2).state;
    s = [s1; s2];
    u1 = output.result.interpsolution.phase(1).control;
    u2 = output.result.interpsolution.phase(2).control;
    u = [u1; u2];
    p1 = output.result.interpsolution.phase(1).costate;
    p2 = output.result.interpsolution.phase(2).costate;
    p = [p1; p2];
else
    t = output.result.interpsolution.phase.time;
    s = output.result.interpsolution.phase.state;
    u = output.result.interpsolution.phase.control;
    p = output.result.interpsolution.phase.costate;
end%if

if (Min_Time && (Sched_Pol || Hybrid_Opt)) || (Min_Fuel && Decoup_Bang) % You have separate results for X/Y vs. Z dynamics
    
    % Fuel use in normalized units (integral of control magnitudes):
    norm_fuel_use = sum([trapz(t_xy,abs(u_xy)) trapz(t_z,abs(u_z))]);
    
    % For Min_Time problem with Sched_Pol, ensure solution vectors are same
    % length by interpolating
    if Min_Time && Sched_Pol
        if t_xy(end) < t_z(end)
            x_interp = interp1(t_xy,x,t_z);
            uFx_interp = interp1(t_xy,uFx,t_z);
            y_interp = interp1(t_xy,y,t_z);
            uFy_interp = interp1(t_xy,uFy,t_z);
            for i = 1:numel(y_interp)
                if isfinite(x_interp(i)) == 0
                    x_interp(i) = 0;
                end%if
                if isfinite(uFx_interp(i)) == 0
                    uFx_interp(i) = 0;
                end%if
                if isfinite(y_interp(i)) == 0
                    y_interp(i) = 0;
                end%if
                if isfinite(uFy_interp(i)) == 0
                    uFy_interp(i) = 0;
                end%if
            end%for
            x = x_interp;
            y = y_interp;
            uFx = uFx_interp;
            uFy = uFy_interp;
            t = t_z;
        elseif t_z(end) < t_xy(end)
            z_interp = interp1(t_z,z,t_xy);
            uFz_interp = interp1(t_z,uFz,t_xy);
            for i = 1:numel(z_interp)
                if isfinite(z_interp(i)) == 0
                    z_interp(i) = 0;
                end%if
                if isfinite(uFz_interp(i)) == 0
                    uFz_interp(i) = 0;
                end%if
            end%for
            z = z_interp;
            uFz = uFz_interp;
            t = t_xy;
        end%if
    end%if
    
    % Plot translation state histories (x,y,z, xd,yd,zd)
    figure('Position',[40 400 1600 400]);
    ax1 = subplot(1,4,1);
    plot(t_xy,s_xy(:,1:4)); hold on;
    plot(t_z,s_z(:,1:2));
    grid on;
    legend('X','Y', 'X dot', 'Y dot', 'Z','Z dot');
    title('Translation State Histories')
    xlabel('Time (s)');
    ylabel('Translation State (m, m/s)');

    % Plot thruster control history (uFx, uFy, uFz)
    ax2 = subplot(1,4,2);
    grid on; hold on;
    if Min_Time && Sched_Pol
        line1 = plot(t_max,  Fmax*ones(size(t_max(:,1))), 'k', 'LineWidth',2);
        line2 = plot(t_max, -Fmax*ones(size(t_max(:,1))), 'k', 'LineWidth',2);
        force_x = plot(t,Fmax.*uFx); 
        force_y = plot(t,Fmax.*uFy); 
        force_z = plot(t,Fmax.*uFz); 
    elseif Min_Fuel
        line1 = plot(t_xy,  Fmax*ones(size(t_xy(:,1))), 'k', 'LineWidth',2);
        line2 = plot(t_z, -Fmax*ones(size(t_z(:,1))), 'k', 'LineWidth',2);
        force_x = plot(t_xy, Fmax.*uFx);
        force_y = plot(t_xy, Fmax.*uFy);
        force_z = plot(t_z,  Fmax.*uFz);
    end%if
    legend([line1 force_x force_y force_z], 'Max/Min Thrust Limits','Radial Thrust','Transverse Thrust','Cross-Track Thrust');
    title('Thruster Control History')
    ylim([-1.1*Fmax 1.1*Fmax])
    xlim([0 t_xy(end)]);
    xlabel('Time (s)');
    ylabel('Control Thrust (N)');

%     % Plot Hamilton w.r.t. time
    ax3 = subplot(1,4,3);
    plot(t_xy, H_XY);
    grid on; hold on;
    title('Hamiltonian for XY dynamics')
    legend('H_{xy}')
    xlabel('Time');
    ylabel('Hamiltonian');

    ax4 = subplot(1,4,4);
    grid on; hold on;
    plot(t_z, H_Z);
    title('Hamiltonian for Z Dynamics')
    legend('H_z')
    xlabel('Time');
    ylabel('Hamiltonian');

    return % Don't use plotting tools below

elseif Incl_Att_Dynamics
            x = s(:,1);
            y = s(:,2);
            z = s(:,3);
            xdot = s(:,4);
            ydot = s(:,5);
            zdot = s(:,6);
            q1 = s(:,7);
            q2 = s(:,8);
            q3 = s(:,9);
            q4 = s(:,10);
            wx = s(:,11);
            wy = s(:,12);
            wz = s(:,13);
        uF  = u(:,1);
        uTx = u(:,2);
        uTy = u(:,3);
        uTz = u(:,4);

    % Calculate the Hamiltonian at each time step
    H = [];
    for i = 1:numel(t)

        q = [q1(i); q2(i); q3(i); q4(i)];       % Assemble quaternions into a vector
        w = [wx(i); wy(i); wz(i)];               % Assemble body rates into a vector
        R1 = q(4).^2 + q(1).^2 - q(2).^2 - q(3).^2;
        R2 = 2.*(q(1).*q(2) - q(3).*q(4));
        R3 = 2.*(q(1).*q(3) + q(2).*q(4));
        R4 = 2.*(q(1).*q(2) + q(3).*q(4));
        R5 = q(4).^2 - q(1).^2 + q(2).^2 - q(3).^2;
        R6 = 2.*(q(2).*q(3) - q(1).*q(4));
        R7 = 2.*(q(1).*q(3) - q(2).*q(4));
        R8 = 2.*(q(2).*q(3) + q(1).*q(4));
        R9 = q(4).^2 - q(1).^2 - q(2).^2 + q(3).^2;
        R = [R1 R2 R3; ...
             R4 R5 R6; ...
             R7 R8 R9];            % Calculate DCM from body frame to LVLH frame
        qdot = 0.5 * (alt_quatmultiply(q', [w' 0]) + alt_quatmultiply([0 0 n 0],q'));

        xd                 = xdot(i);
        yd                 = ydot(i);
        zd                 = zdot(i);
        % Clohessy-Wiltshire Equations of Motion (with forcing term)
        xddot              = 2*n.*ydot(i) + 3*n^2.*x(i) + b.*uF(i) * ([1 0 0]*R*Thrust_dir);
        yddot              = -2*n.*xdot(i)              + b.*uF(i) * ([0 1 0]*R*Thrust_dir); 
        zddot              = -n^2.*z(i)                 + b.*uF(i) * ([0 0 1]*R*Thrust_dir);
%         xddot              = 2*n.*ydot(i) + 3*n^2.*x(i) + b.*([1 0 0]*R*[uFx(i) uFy(i) uFz(i)]');
%         yddot              = -2*n.*xdot(i)              + b.*([0 1 0]*R*[uFx(i) uFy(i) uFz(i)]'); 
%         zddot              = -n^2.*z(i)                 + b.*([0 0 1]*R*[uFx(i) uFy(i) uFz(i)]');

        % Euler's Moment Equations (with forcing term)
        % (NOTE: I(1) = Ixx, I(5) == Iyy, I(9) == Izz)
        wxdot              = (I(5) - I(9))/I(1) .* wy(i) .* wz(i) + (Tmax/I(1) .* uTx(i));
        wydot              = (I(9) - I(1))/I(5) .* wx(i) .* wz(i) + (Tmax/I(5) .* uTy(i));
        wzdot              = (I(1) - I(5))/I(9) .* wx(i) .* wy(i) + (Tmax/I(9) .* uTz(i));

        s_dot = [xd; yd; zd; xddot; yddot; zddot; ...
                 qdot'; ...
                 wxdot; wydot; wzdot]; 

        H(i) = 1 + p(i,:)*s_dot;
    end%for  % End Hamiltonian calculation

    % Plot translation state histories (x,y,z, xd,yd,zd)
    figure('Position',[40 400 1600 400]);
    ax1 = subplot(2,3,1);
    plot(t,s(:,1:6)); hold on;
    grid on;
    legend('X','Y','Z', 'X dot', 'Y dot', 'Z dot');
    title('Translation State Histories')
    xlabel('Time (s)');
    ylabel('Translation State (m, m/s)');

    % Plot quaternion state histories (q1,q2,q3,q4)
    ax2 = subplot(2,3,2);
    plot(t,s(:,7:10)); hold on;
    grid on;
    if Period_Match
        ph_mk = plot([t1(end) t1(end)],[-1 1],'--k');
%         ph_mk = plot([t2(end) t2(end)],[-1 1],'--k');
        legend('q1', 'q2', 'q3', 'q4','Phase Boundary');
    else
        legend('q1', 'q2', 'q3', 'q4');
    end%if
    ylim([-1 1])
    title('Quaternion Histories')
    xlabel('Time (s)');
    ylabel('Quaternion Magnitude');

    % Plot rotation rate histories (wx,wy,wz)
    ax3 = subplot(2,3,3);
    plot(t,s(:,11:13)); hold on;
    grid on;
    if Period_Match
        ph_mk = plot([t1(end) t1(end)],[-1 1],'--k');
%         ph_mk = plot([t2(end) t2(end)],[-1 1],'--k');
        legend('wx','wy','wz','Phase Boundary');
    else
        legend('wx','wy','wz');
    end%if
    title('Rotation Rate Histories')
    xlabel('Time (s)');
    ylabel('Rotation Rates (rad/s)');

    % Plot thruster control history (uF)
    ax4 = subplot(2,3,4);
    grid on; hold on;
    line1 = plot(t, Fmax*ones(size(u(:,1))), 'k', 'LineWidth',2);
    line2 = plot(t, zeros(size(u(:,1))), 'k', 'LineWidth',2);
    force = plot(t,Fmax.*uF);
    if Period_Match
        ph_mk = plot([t1(end) t1(end)],[-Fmax Fmax],'--k');
%         ph_mk = plot([t2(end) t2(end)],[-Fmax Fmax],'--k');
        legend([line1 force ph_mk], 'Max/Min Thrust Limits','Thrust','Phase Boundary');
    else
        legend([line1 force], 'Max/Min Thrust Limits','Thrust');
    end%if
    title('Thruster Control History')
    ylim([-0.1*Fmax 1.1*Fmax])
    xlabel('Time (s)');
    ylabel('Control Thrust (N)');

%     ax4 = subplot(2,3,4);
%     grid on; hold on;
%     line1 = plot(t, Fmax*ones(size(u(:,1))), 'k', 'LineWidth',2);
%     line2 = plot(t, -Fmax*ones(size(u(:,1))), 'k', 'LineWidth',2);
%     forcex = plot(t,Fmax.*uFx);
%     forcey = plot(t,Fmax.*uFy);
%     forcez = plot(t,Fmax.*uFz);
%     legend([line1 forcex forcey forcez], 'Max/Min Thrust Limits','Radial Thrust','Transverse Thrust','Cross-Track Thrust');
%     title('Thruster Control History')
%     ylim([-1.1*Fmax 1.1*Fmax])
%     xlabel('Time (s)');
%     ylabel('Control Thrust (N)');



    % Plot reaction wheels torque histories (uTx,uTy,uTz)
    ax5 = subplot(2,3,5);
    grid on; hold on;
    line1 = plot(t, -Tmax*ones(size(u(:,1))), 'k', 'LineWidth',2);
    line2 = plot(t, Tmax*ones(size(u(:,1))), 'k', 'LineWidth',2);
    torque_x = plot(t,Tmax.*uTx);
    torque_y = plot(t,Tmax.*uTy);
    torque_z = plot(t,Tmax.*uTz);
    if Period_Match
        ph_mk = plot([t1(end) t1(end)],[-Tmax Tmax],'--k');
%         ph_mk = plot([t2(end) t2(end)],[-Tmax Tmax],'--k');
        legend([line1 torque_x torque_y torque_z ph_mk],'Max/Min Torque Limits','Tx','Ty','Tz','Phase Boundary');
    else
        legend([line1 torque_x torque_y torque_z],'Max/Min Torque Limits','Tx','Ty','Tz');
    end%if
    title('Reaction Wheel Torque History')
    ylim([-1.1*Tmax 1.1*Tmax])
    xlabel('Time (s)');
    ylabel('Control Torque (N-m)');

    % Plot Hamilton w.r.t. time
    ax6 = subplot(2,3,6);
    plot(t, H);
    grid on;
    title('Hamiltonian')
    xlabel('Time');
    ylabel('Hamiltonian');

    % Make it easier to zoom in synchronously on plots 
    linkaxes([ax1, ax2, ax3, ax4, ax5, ax6],'x');
    %%
    % Call the function that draws/animates trajectory
%     quiverPlotSat(MakeAnimation);


else % If not including attitude dynamics:

        x = s(:,1);
        y = s(:,2);
        z = s(:,3);
        xdot = s(:,4);
        ydot = s(:,5);
        zdot = s(:,6);
    uFx = u(:,1);
    uFy = u(:,2);
    uFz = u(:,3);

    if (Target_Ellipse && coast) || Period_Match
        norm_fuel_use = sum([trapz(t1,abs(u1(:,1))) trapz(t1,abs(u1(:,2))) trapz(t1,abs(u1(:,3)))]);
    else
        % Fuel use in normalized units (integral of control magnitudes):
        norm_fuel_use = sum([trapz(t,abs(uFx)) trapz(t,abs(uFy)) trapz(t,abs(uFz))]);
    end%if
    
    % Calculate the Hamiltonian at each time step
    H = [];
    for i = 1:numel(t)

        xd                 = xdot(i);
        yd                 = ydot(i);
        zd                 = zdot(i);
        if Nonlinear
            if e ~= 0 % Time-varying dynamics - lookup true anomaly 
                nu_interp = interp1qr(nu_table(:,1),nu_table(:,2),mod(t(i),2*pi/n));
                r_tgt = (a*(1-e^2))./(1+e.*cos(nu_interp));
                h = sqrt(mu .* a*(1-e^2));
                nu_dot = h./(r_tgt.^2);
                nu_ddot = -2.*(h.^2.*e.*sin(nu_interp))./(r_tgt.^4.*(1+e.*cos(nu_interp)));
                rd_tgt = (r_tgt.*nu_dot.*e.*sin(nu_interp))./(1+e.*cos(nu_interp));
            else % Time-invariant dynamics, true anomaly rate of change is constant.
                r_tgt = a;   
                rd_tgt = 0;
                h = sqrt(mu * a);
                nu_dot = n;
                nu_ddot = 0;
                rd_tgt = 0;
            end%if
            r_ch  = sqrt((r_tgt + x(i)).^2 + y(i).^2 + z(i).^2);

            % nonlinear Keplerian Equations of relative motion w/ forcing 
            xddot = 2.*nu_dot.*ydot(i) + nu_ddot.*y(i) + nu_dot.^2.*x(i) + mu./(r_tgt.^2) - mu./(r_ch.^3).*(r_tgt + x(i)) + b.*uFx(i);
            yddot = -2.*nu_dot.*xdot(i) - nu_ddot.*x(i) + nu_dot.^2.*y(i) - mu./(r_ch.^3).*y(i) + b.*uFy(i);
            zddot = -mu./(r_ch.^3).*z(i) + b.*uFz(i);
        else
            % Clohessy-Wiltshire Equations of Motion (with forcing term)
            xddot              = 2*n.*ydot(i) + 3*n^2.*x(i) + b.*uFx(i);
            yddot              = -2*n.*xdot(i)              + b.*uFy(i); 
            zddot              = -n^2.*z(i)                 + b.*uFz(i);
        end%if
        s_dot = [xd; yd; zd; xddot; yddot; zddot];
        
        if Min_Time
            H(i) = 1 + p(i,:)*s_dot;
        elseif Min_Energy
            if Nonlinear && enforce_CW_ellipse
                x_ellipse = CW_alpha/2.*sin(n.*t(i));
                x_error = x(i) - x_ellipse;
                y_ellipse = CW_alpha.*cos(n.*t(i));
                y_error = y(i) - y_ellipse;
                z_ellipse = CW_alpha/2*tand(CW_zeta).*sin(n.*t(i));
                z_error = z(i) - z_ellipse;
                H(i) = 0.5*1e-6*(x_error.^2 + y_error.^2 + z_error.^2) + ...
                        0.5.*(uFx(i).^2 + uFy(i).^2 + uFz(i).^2) + ...
                        p(i,:)*s_dot;
            else
                H(i) = 0.5*uFx(i)^2 + 0.5*uFy(i)^2 + 0.5*uFz(i)^2 + p(i,:)*s_dot;
            end%if
        elseif Min_Fuel
            H(i) = uFx(i) + uFy(i) + uFz(i) + p(i,:)*s_dot;
        elseif Quadratic_Cost
            H(i) = 0.1*Q_R_ratio*(x(i)^2 + y(i)^2 + z(i)^2 + ...
                                  xdot(i)^2 + ydot(i)^2 + zdot(i)^2) + ...
                   0.1*(uFx(i)^2 + uFy(i)^2 + uFz(i)^2) + ...
                   p(i,:)*s_dot;
        end%if
    end%for

   % Plot translation state histories (x,y,z, xd,yd,zd)
    figure('Position',[40 400 1600 400]);
    ax1 = subplot(1,3,1);
    plot(t,s(:,1:6)); hold on;
    if num_phases > 1 && Min_Fuel
        plot([t1(end), t1(end)], [-0.5*max(max(s)),0.5*max(max(s))],'--k')
        plot([t2(end), t2(end)], [-0.5*max(max(s)),0.5*max(max(s))],'--k')
        legend('X','Y','Z', 'X dot', 'Y dot', 'Z dot');
    elseif Target_Ellipse || Period_Match
        plot([t1(end), t1(end)], [-0.5*max(max(s)),0.5*max(max(s))],'--k')
        legend('X','Y','Z', 'X dot', 'Y dot', 'Z dot','Phase Boundary');
    else
        legend('X','Y','Z', 'X dot', 'Y dot', 'Z dot');
    end%if
    grid on;
    xlim([0 t(end)]);
    title('Translation State Histories')
    xlabel('Time (s)');
    ylabel('Translation State (m, m/s)');

    % Plot thruster control history (uFx, uFy, uFz)
    ax2 = subplot(1,3,2);
    grid on; hold on;
    line1 = plot(t,  Fmax*ones(size(u(:,1))), 'k', 'LineWidth',2);
    line2 = plot(t, -Fmax*ones(size(u(:,1))), 'k', 'LineWidth',2);
    if num_phases > 1 && Min_Fuel
        if Decoup_Bang
            force_x = plot(t1_xy, Fmax.*u1_xy(:,1),'b');
            force_y = plot(t1_xy, Fmax.*u1_xy(:,2),'r');
            force_z = plot(t1_z, Fmax.*u1_z(:,1),'y');
            plot(t2_xy, Fmax.*u2_xy(:,1),'b')
            plot(t2_xy, Fmax.*u2_xy(:,2),'r')
            plot(t2_z, Fmax.*u2_z(:,1),'y')
            plot(t3_xy, Fmax.*u3_xy(:,1),'b')
            plot(t3_xy, Fmax.*u3_xy(:,2),'r')
            plot(t3_z, Fmax.*u3_z(:,1),'y')
            ph_mk_xy = plot([t1_xy(end), t1_xy(end)], [-0.5*Fmax,0.5*Fmax],'--o');
            ph_mk_xy = plot([t2_xy(end), t2_xy(end)], [-0.5*Fmax,0.5*Fmax],'--o');
            ph_mk_z = plot([t1_z(end), t1_z(end)], [-0.5*Fmax,0.5*Fmax],'--y');
            ph_mk_z = plot([t2_z(end), t2_z(end)], [-0.5*Fmax,0.5*Fmax],'--y');
            legend([line1 force_x force_y force_z ph_mk_xy ph_mk_z], 'Max/Min Thrust Limits','Radial Thrust','Transverse Thrust','Cross-Track Thrust', 'XY Phase Separation Marker', 'Z Phase Separation Marker');
        else
            force_x = plot(t1, Fmax.*u1(:,1),'b');
            force_y = plot(t1, Fmax.*u1(:,2),'r');
            force_z = plot(t1, Fmax.*u1(:,3),'y');
            plot(t2, Fmax.*u2(:,1),'b')
            plot(t2, Fmax.*u2(:,2),'r')
            plot(t2, Fmax.*u2(:,3),'y')
            plot(t3, Fmax.*u3(:,1),'b')
            plot(t3, Fmax.*u3(:,2),'r')
            plot(t3, Fmax.*u3(:,3),'y')
            ph_mk = plot([t1(end), t1(end)], [-0.5*Fmax,0.5*Fmax],'--k');
            ph_mk = plot([t2(end), t2(end)], [-0.5*Fmax,0.5*Fmax],'--k');
            legend([line1 force_x force_y force_z ph_mk], 'Max/Min Thrust Limits','Radial Thrust','Transverse Thrust','Cross-Track Thrust', 'Phase Separation Marker');
        end%if
    elseif Target_Ellipse || Period_Match
        force_x = plot(t, Fmax.*uFx);
        force_y = plot(t, Fmax.*uFy);
        force_z = plot(t, Fmax.*uFz);
        ph_mk = plot([t1(end), t1(end)], [-0.5*Fmax,0.5*Fmax],'--k');
        legend([line1 force_x force_y force_z ph_mk], 'Max/Min Thrust Limits','Radial Thrust','Transverse Thrust','Cross-Track Thrust', 'Phase Separation Marker');
    else % Single-phase
        force_x = plot(t, Fmax.*uFx);
        force_y = plot(t, Fmax.*uFy);
        force_z = plot(t, Fmax.*uFz);
        legend([line1 force_x force_y force_z], 'Max/Min Thrust Limits','Radial Thrust','Transverse Thrust','Cross-Track Thrust');
    end%if
    title('Thruster Control History')
    ylim([-1.1*Fmax 1.1*Fmax])
    xlim([0 t(end)]);
    xlabel('Time (s)');
    ylabel('Control Thrust (N)');

    % Plot Hamilton w.r.t. time
    ax3 = subplot(1,3,3);
    plot(t, H);
    grid on;
    title('Hamiltonian')
    xlabel('Time');
    ylabel('Hamiltonian');
    xlim([0 t(end)]);
    % Make it easier to zoom in synchronously on plots 
    linkaxes([ax1, ax2, ax3],'x');
    %%
    % Call the function that draws/animates trajectory
    % quiverPlotSat(MakeAnimation);

end%if Pos+Attitude vs. Position Only formulation
       
%%
if Plot3D_Trajectory
    
    figure();
    plot3(x,y,z,'c','LineWidth',2)
    grid on
    hold on
    if num_phases >1
    plot3(s1(end,1),s1(end,2),s1(end,3),'ok');
    end%if

    plot3(0,0,0,'ok','LineWidth',1,'MarkerFaceColor','k')
    axis equal
    xlabel('Radial [m]')
    ylabel('In-Track [m]')
    zlabel('Cross-Track [m]')
    
    scale = max(max(abs(s(:,1:3))));
    [cylX, cylY, cylZ] = cylinder([0 scale],100);
    cylZ = cylZ*scale*sqrt(cotd(FOV/2)^2);
    cone = surf(-cylZ,cylX,cylY); 
    set(cone,'FaceColor',[0.9 0.9 0.9],'EdgeColor',[0.9 0.9 0.9],'FaceLighting','GOURAUD','FaceAlpha',0.5,'EdgeAlpha',0.5);
   
    if Avoid_FOV
        set(cone,'FaceColor',[0.8 0.2 0.3],'EdgeColor',[0.8 0.2 0.3],'FaceLighting','GOURAUD','FaceAlpha',1);
        distance_to_fov = s(:,1) + sqrt(cotd(FOV/2)^2.*(s(:,2).^2 + s(:,3).^2));
        figure();
        plot(t,distance_to_fov);
        xlabel('Time(s)');
        ylabel('Radial Displacement from FOV Cone (m)');
        grid on;
    end
end%if
    
    
    
    
    
    