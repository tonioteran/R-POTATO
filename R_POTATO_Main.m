%---------------------------------------------------%
% Zachary K. Funke                                 %
%---------------------------------------------------%
% The problem solved here is given as follows:      %
%   Minimize J = int( u' * R * u ) from 0 to T      %
% subject to the simplified relative dynamics of    %
% nearby orbiting objects given by the Clohessy -   %
% Wiltshire Equations of motion                     %
%    x'' - 2*n*y' = a_x                             %
%    y'' + 2*n*x' - 3*n^2*y = a_y                   %
%    z'' + n^2*z = a_z                              %
% and the boundary conditions                       %
%    x,y,z(0) = s_0                                 %
%    x,y,z(t_f) = 0, t_f = T = FREE                 %
%---------------------------------------------------%
% clear all; clc

switch InPlane_Only + OutOfPlane_Only
    case 0
        if Incl_Att_Dynamics
            state_dim = 13;
            ctrl_dim  = 4; 
            %ctrl_dim = 6;
        else
            state_dim = 6;
            ctrl_dim = 3;
        end%if
    case 1
        if InPlane_Only
            state_dim = 4;
            ctrl_dim = 2;
        elseif OutOfPlane_Only
            state_dim = 2;
            ctrl_dim = 1;
        end%if
    case 2
        error('Both InPlane_Only and OutOfPlane_Only are set.')
end%switchcase

auxdata.n           = n;
auxdata.mu          = mu;
auxdata.m           = m;
auxdata.dim         = state_dim;
auxdata.decoup      = Decoup_Bang;
auxdata.mintime     = Min_Time;
auxdata.minfuel     = Min_Fuel;
auxdata.minenergy   = Min_Energy;
auxdata.quadcost    = Quadratic_Cost;
auxdata.qr_ratio    = Q_R_ratio;
auxdata.Fmax        = Fmax;
auxdata.Thrust_dir  = Thrust_dir;
auxdata.Tmax        = Tmax;
auxdata.I           = I;
auxdata.num_phases  = num_phases;
auxdata.ellipse     = Target_Ellipse;
auxdata.alpha       = alpha;
auxdata.zeta        = zeta;
auxdata.gamma       = gamma;
auxdata.coast       = coast;
auxdata.range_path  = range_path;
auxdata.range       = range;
auxdata.period_path = period_path;
auxdata.center_path = center_path;
auxdata.center      = center;
auxdata.vbar_offset = Vbar_Offset;
auxdata.Nonlinear   = Nonlinear;
auxdata.e           = e;
auxdata.avoid       = Avoid_FOV;
auxdata.FOV         = FOV;
auxdata.FOV_buff    = FOV_buff;
auxdata.P_Match     = Period_Match;
auxdata.T_Point     = Target_Point;
auxdata.de_ON       = de_ON;
auxdata.x_size      = x_size;
auxdata.dw_ON       = dw_ON;
auxdata.dw          = dw;
auxdata.di_ON       = di_ON;
auxdata.di          = di;
auxdata.dG_ON       = dG_ON;
auxdata.dG          = dG;
auxdata.enforceCW   = (Nonlinear && enforce_CW_ellipse);    
auxdata.CW_alpha    = CW_alpha;
auxdata.CW_zeta     = CW_zeta;
auxdata.mnvr_epoch  = Mnvr_Epoch;

paths = [range_path, ...
         period_path, ...
         center_path];
     
if e ~=0 
    auxdata.nu_table    = nu_table;
end%if

avoidance_constr_lower = 0;
avoidance_constr_upper = 10e5;


DAY      = 86400;
SID_DAY  = 86164;

t0 = 0; 
tfmin = 1; 
if (Min_Energy || Min_Fuel)
    tfmax = time_constr;
else
    tfmax = Max_Time;
end%if


if Incl_Att_Dynamics %-------------------% 6-DOF, 13 state vars
    s_0   = [ps_0; rs_0]'; 

    s_min = [-MaxRange -MaxRange -MaxRange ...
              -MaxVel  -MaxVel -MaxVel    ...
             -1 -1 -1 -1       ...
             -1 -1 -1       ]; 

    s_max = [ MaxRange  MaxRange  MaxRange ...
               MaxVel   MaxVel   MaxVel       ...
             1 1 1 1           ...
             1 1 1          ];

    s_f   = [ 0 0 0   ...   % Rendezvous'd
              0 0 0   ...   % no transl velocity
              0 0 0 1 ...   % end with axes aligned 
              0 0 0      ]; % no rotational velocity
    umin = [0 -1 -1 -1]; 
    umax = [1  1  1  1];
elseif InPlane_Only %-------------------% 2-DOF, 4 state vars
    s_0   = [ps_0(1:2); ps_0(4:5)]';
    s_min = [-MaxRange -MaxRange ...
              -MaxVel  -MaxVel  ];
    s_max = [ MaxRange  MaxRange ...
               MaxVel   MaxVel  ];      
    s_f   = [ 0 0    ...
              0 0   ];
    umin = [-1 -1 ]; 
    umax = [ 1  1 ];
elseif OutOfPlane_Only || Decoup_Bang %-----% 1-DOF, 2 state vars
    s_0   = [ps_0(3); ps_0(6)]';
    s_min = [-MaxRange ...
              -MaxVel  ];
    s_max = [ MaxRange ...
               MaxVel  ];      
    s_f   = [ 0   ...
              0  ];
    umin = [-1 ]; 
    umax = [ 1 ];
else %--------------------------- % 3-DOF, 6 state vars
    s_0   = [ps_0]';
    s_min = [-MaxRange -MaxRange -MaxRange ...
              -MaxVel  -MaxVel -MaxVel  ];
    s_max = [ MaxRange  MaxRange  MaxRange ...
               MaxVel   MaxVel   MaxVel ];      
    s_f   = [ 0 0 0    ...
              0 0 0   ];
    umin = [-1 -1 -1 ]; 
    umax = [1  1  1 ];
end%if

%-------------------------------------------------------------------------%
%----------------------- Setup for Problem Bounds ------------------------%
%-------------------------------------------------------------------------%
for p = 1:num_phases
bounds.phase(p).initialtime.lower = t0; 
bounds.phase(p).initialtime.upper = t0;

if (Min_Time || Quadratic_Cost) && ~time_constr
    bounds.phase(p).finaltime.lower = tfmin;   % final time FREE (objective) 
    bounds.phase(p).finaltime.upper = tfmax;
elseif (Min_Fuel || Min_Energy) && ~time_constr
    bounds.phase(p).finaltime.lower = tfmax;   % final time FIXED (constraint)
    bounds.phase(p).finaltime.upper = tfmax;
elseif time_constr
    bounds.phase(p).finaltime.lower = tfmax;   % final time FIXED to match limiting trajectory
    bounds.phase(p).finaltime.upper = tfmax;
end%if

    if Min_Fuel  && num_phases > 1 % Need three phases for bang-off-bang.
        if p == 1 
            bounds.phase(p).initialtime.lower = t0; 
            bounds.phase(p).initialtime.upper = t0;
            bounds.phase(p).finaltime.lower = 0;   
            bounds.phase(p).finaltime.upper = Bang_Width;
            bounds.phase(p).initialstate.lower = [s_0]; 
            bounds.phase(p).initialstate.upper = [s_0]; 
            bounds.phase(p).state.lower = [s_min]; 
            bounds.phase(p).state.upper = [s_max]; 
            bounds.phase(p).finalstate.lower = [s_min]; 
            bounds.phase(p).finalstate.upper = [s_max]; 
            if Target_Point && ~Incl_Att_Dynamics % 3-DOF State regulation problem
                for j = 1:ctrl_dim
                    if ps_0(j) > 0 
                        bounds.phase(p).control.lower(j) = umax(j); 
                        bounds.phase(p).control.upper(j) = umax(j);
                    elseif ps_0(j) < 0
                        bounds.phase(p).control.lower(j) = umin(j); 
                        bounds.phase(p).control.upper(j) = umin(j);
                    elseif ps_0(j) == 0
                        bounds.phase(p).control.lower(j) = 0; 
                        bounds.phase(p).control.upper(j) = 0;
                    end%if 
                end%for
            end%if
            bounds.phase(p).integral.lower = 0;
            bounds.phase(p).integral.upper = 7e5;
        elseif p == 2
            bounds.phase(p).initialtime.lower = 0; 
            bounds.phase(p).initialtime.upper = Bang_Width;
            bounds.phase(p).finaltime.lower = tfmax-Bang_Width;   
            bounds.phase(p).finaltime.upper = tfmax;
            bounds.phase(p).initialstate.lower = [s_min]; 
            bounds.phase(p).initialstate.upper = [s_max]; 
            bounds.phase(p).state.lower = [s_min]; 
            bounds.phase(p).state.upper = [s_max]; 
            bounds.phase(p).finalstate.lower = [s_min]; 
            bounds.phase(p).finalstate.upper = [s_max]; 
            if Incl_Att_Dynamics % Still allow attitude control during coast phase
                bounds.phase(p).control.lower = [0 umin(2:end)]; 
                bounds.phase(p).control.upper = [0 umax(2:end)];
            else
                bounds.phase(p).control.lower = zeros(size(umin)); 
                bounds.phase(p).control.upper = zeros(size(umax));
            end%if
            bounds.phase(p).integral.lower = 0;
            bounds.phase(p).integral.upper = 7e5;
        elseif p == 3
            bounds.phase(p).initialtime.lower = tfmax-Bang_Width; 
            bounds.phase(p).initialtime.upper = tfmax;
            bounds.phase(p).finaltime.lower = tfmax;   
            bounds.phase(p).finaltime.upper = tfmax;
            bounds.phase(p).initialstate.lower = [s_min]; 
            bounds.phase(p).initialstate.upper = [s_max]; 
            bounds.phase(p).state.lower = [s_min]; 
            bounds.phase(p).state.upper = [s_max]; 
            if Incl_Att_Dynamics
                bounds.phase(p).finalstate.lower = [s_f(1:6), s_min(7:10), s_f(11:13)]; % final state FIXED! (except for quats)
                bounds.phase(p).finalstate.upper = [s_f(1:6), s_max(7:10), s_f(11:13)]; % final state FIXED! (except for quats)
            else
                bounds.phase(p).finalstate.lower = [s_f]; 
                bounds.phase(p).finalstate.upper = [s_f]; 
            end%if
            if Target_Point && ~Incl_Att_Dynamics % 3-DOF State regulation problem
                for j = 1:ctrl_dim
                    if ps_0(j) > 0 
                        bounds.phase(p).control.lower(j) = umin(j); 
                        bounds.phase(p).control.upper(j) = umin(j);
                    elseif ps_0(j) < 0
                        bounds.phase(p).control.lower(j) = umax(j); 
                        bounds.phase(p).control.upper(j) = umax(j);
                    elseif ps_0(j) == 0
                        bounds.phase(p).control.lower(j) = 0; 
                        bounds.phase(p).control.upper(j) = 0;
                    end%if 
                end%for
            end%if
            bounds.phase(p).integral.lower = 0;
            bounds.phase(p).integral.upper = 7e5;
        end%if
    elseif Min_Fuel && num_phases == 1
            bounds.phase(p).initialtime.lower = t0; 
            bounds.phase(p).initialtime.upper = t0;
            bounds.phase(p).finaltime.lower = tfmax;   
            bounds.phase(p).finaltime.upper = tfmax;
            bounds.phase(p).initialstate.lower = [s_0]; 
            bounds.phase(p).initialstate.upper = [s_0]; 
            bounds.phase(p).state.lower = [s_min]; 
            bounds.phase(p).state.upper = [s_max]; 
            bounds.phase(p).finalstate.lower = [s_f]; 
            bounds.phase(p).finalstate.upper = [s_f]; 
            bounds.phase(p).control.lower = [umin]; 
            bounds.phase(p).control.upper = [umax];
%             bounds.phase(p).path.lower = zeros(1, ctrl_dim);
%             bounds.phase(p).path.upper = zeros(1, ctrl_dim);
            bounds.phase(p).integral.lower = 0;
            bounds.phase(p).integral.upper = 7e5;
    elseif ~Incl_Att_Dynamics && (Target_Ellipse || Period_Match) % Include second phase for ellipse
        if p == 1   % Approaching Ellipse
            bounds.phase(p).initialstate.lower = [s_0]; 
            bounds.phase(p).initialstate.upper = [s_0]; 
            bounds.phase(p).state.lower = [s_min]; 
            bounds.phase(p).state.upper = [s_max]; 
            bounds.phase(p).finalstate.lower = s_min; % Don't constrain first phase final state, this is handled by Endpoint
            bounds.phase(p).finalstate.upper = s_max; 
            bounds.phase(p).control.lower = umin; 
            bounds.phase(p).control.upper = umax;
            bounds.phase(p).integral.lower = 0;
            bounds.phase(p).integral.upper = 5e7;
        elseif p == 2   % On Ellipse
            if Min_Energy
                bounds.phase(p).initialtime.lower = tfmax; 
                bounds.phase(p).initialtime.upper = tfmax;
                bounds.phase(p).finaltime.lower = tfmax+ 1.1*SID_DAY;   
                bounds.phase(p).finaltime.upper = tfmax+ 1.1*SID_DAY;
            elseif Min_Time
                bounds.phase(p).initialtime.lower = tfmin; 
                bounds.phase(p).initialtime.upper = tfmax;
                bounds.phase(p).finaltime.lower = tfmin+ 2*DAY;   
                bounds.phase(p).finaltime.upper = tfmax+ 2*DAY;
            end%if
            bounds.phase(p).initialstate.lower = [s_min]; 
            bounds.phase(p).initialstate.upper = [s_max]; 
            bounds.phase(p).state.lower = [s_min]; 
            bounds.phase(p).state.upper = [s_max]; 
            bounds.phase(p).finalstate.lower = [s_min]; % final state 
            bounds.phase(p).finalstate.upper = [s_max]; % final state 
            if coast || Period_Match % Disallow thrusting on ellipse
                bounds.phase(p).control.lower = zeros(size(umin)); 
                bounds.phase(p).control.upper = zeros(size(umax));
            elseif ~coast % Allow thrusting on ellipse (and follow path constraints)
                bounds.phase(p).control.lower = umin; 
                bounds.phase(p).control.upper = umax;
                if Nonlinear % Follow path constraint for NLEM
                    bounds.phase(2).path.lower = [0  ];
                    bounds.phase(2).path.upper = [0  ];
                else % Follow path constraints for CW
                    bounds.phase(2).path.lower = zeros(size(paths));
                    bounds.phase(2).path.upper = zeros(size(paths));
                end%if
            end%if
            bounds.phase(p).integral.lower = 0;
            bounds.phase(p).integral.upper = 5e14;
        end%if
    elseif Incl_Att_Dynamics && Period_Match
        thrust_time_max = ThrustTime;
        if p == 1
            bounds.phase(p).finaltime.lower = 0;   
            bounds.phase(p).finaltime.upper = thrust_time_max;
            bounds.phase(p).initialstate.lower = s_0; 
            bounds.phase(p).initialstate.upper = s_0; 
            bounds.phase(p).state.lower = [s_min(1:13)]; 
            bounds.phase(p).state.upper = [s_max(1:13)]; 
            bounds.phase(p).finalstate.lower = [s_min(1:10) 0 0 0]; 
            bounds.phase(p).finalstate.upper = [s_max(1:10) 0 0 0]; 
            bounds.phase(p).control.lower = umin; 
            bounds.phase(p).control.upper = umax;
            bounds.phase(p).integral.lower = 0;
            bounds.phase(p).integral.upper = 5e7;
        elseif p == 2
            bounds.phase(p).initialtime.lower = 0; 
            bounds.phase(p).initialtime.upper = thrust_time_max;
            bounds.phase(p).finaltime.lower = SID_DAY;%*SID_DAY;   
            bounds.phase(p).finaltime.upper = SID_DAY;%*SID_DAY
            bounds.phase(p).initialstate.lower = s_min; 
            bounds.phase(p).initialstate.upper = s_max; 
            bounds.phase(p).state.lower = s_min; 
            bounds.phase(p).state.upper = s_max; 
            bounds.phase(p).finalstate.lower = s_min; 
            bounds.phase(p).finalstate.upper = s_max; 
            bounds.phase(p).control.lower = [0 0 0 0]; 
            bounds.phase(p).control.upper = [0 0 0 0];
            bounds.phase(p).integral.lower = 0;
            bounds.phase(p).integral.upper = 5e7;
        end%if
    else % Normal single-phase case.
        bounds.phase(p).initialstate.lower = [s_0]; 
        bounds.phase(p).initialstate.upper = [s_0]; 
        bounds.phase(p).state.lower = [s_min]; 
        bounds.phase(p).state.upper = [s_max]; 
        if Incl_Att_Dynamics
            bounds.phase(p).finalstate.lower = [s_f(1:6), s_min(7:10), s_f(11:13)]; % final state FIXED! (except for quats)
            bounds.phase(p).finalstate.upper = [s_f(1:6), s_max(7:10), s_f(11:13)]; % final state FIXED! (except for quats)
        elseif Nonlinear && enforce_CW_ellipse
            bounds.phase(p).finalstate.lower = s_min; % final state FIXED! (except for quats)
            bounds.phase(p).finalstate.upper = s_max;  % final state FIXED! (except for quats)
        else
            bounds.phase(p).finalstate.lower = [s_f]; % final state FIXED! (except for quats)
            bounds.phase(p).finalstate.upper = [s_f]; % final state FIXED! (except for quats)
        end%if
        if Avoid_FOV && Target_Point
            bounds.phase(p).path.lower = avoidance_constr_lower;
            bounds.phase(p).path.upper = avoidance_constr_upper;
        end%if
            bounds.phase(p).control.lower = umin; 
            bounds.phase(p).control.upper = umax;
        bounds.phase(p).integral.lower = 0;
        bounds.phase(p).integral.upper = 5e14;
    end%if
end%for

%-------------------------------------------------------------------------%
%---------------------- Provide Guess of Solution ------------------------%
%-------------------------------------------------------------------------%
for p = 1:num_phases
    
    if Min_Fuel % Guesses for 3 phases: bang-off-bang
        if p == 1
            guess.phase(p).time    = [t0; t0+1]; 
            guess.phase(p).state   = [s_0; s_0];
            guess.phase(p).control = [umax; zeros(size(umin))];
            guess.phase(p).integral = 5;
        elseif p == 2
            guess.phase(p).time    = [Bang_Width; tfmax-Bang_Width]; 
            guess.phase(p).state   = [s_0; s_0];
            guess.phase(p).control = [zeros(size(umin)); zeros(size(umax))];
            guess.phase(p).integral = 25;
        elseif p == 3
            guess.phase(p).time    = [tfmax-1; tfmax]; 
            guess.phase(p).state   = [s_0; s_min];
            guess.phase(p).control = [zeros(size(umin)); umax];
            guess.phase(p).integral = 50;
        end%if        
    elseif Target_Ellipse || Period_Match % Guesses for 2 phases: approach, ellipse
        if p == 1
            guess.phase(p).time    = [t0; 100]; 
            guess.phase(p).state   = [s_0; ones(size(s_0))];
            guess.phase(p).control = [umax; zeros(size(umin))];
            guess.phase(p).integral = 5;
        elseif p == 2
            guess.phase(p).time    = [100; tfmax]; 
            if Incl_Att_Dynamics
                guess.phase(p).state   = [[ones(1,6),zeros(1,3),ones(1,4)]; [ones(1,6),zeros(1,3),ones(1,4)]];       
            else
                guess.phase(p).state   = [ones(size(s_0)); ones(size(s_0))];
            end%if
            guess.phase(p).control = [zeros(size(umin)); zeros(size(umax))];
            guess.phase(p).integral = 25;
        end%if
    elseif Incl_Att_Dynamics && Target_Point 
        guess.phase(p).time    = [0; tfmax]; 
        guess.phase(p).state   = [s_0; s_f];
        guess.phase(p).control = [umin; umax];
        guess.phase(p).integral = 300;
    else % Normal single-phase guesses
        guess.phase(p).time    = [0; tfmax]; 
        guess.phase(p).state   = [s_0; ones(size(s_0))];
        guess.phase(p).control = [umin; umax];
        guess.phase(p).integral = 300;
    end%if
    
end%for

%-------------------------------------------------------------------------%
%------------- Set up Event Constraints That Link Phases -----------------%
%-------------------------------------------------------------------------%

if Min_Fuel && num_phases > 1 % Have three phases with only linkage events.
    bounds.eventgroup(1).lower = zeros(1,state_dim+1); % "Link states and times from phase 1 -> phase 2"
    bounds.eventgroup(1).upper = zeros(1,state_dim+1);
    bounds.eventgroup(2).lower = zeros(1,state_dim+1); % "Link states and times from phase 2 -> phase 3"
    bounds.eventgroup(2).upper = zeros(1,state_dim+1);
    bounds.eventgroup(3).lower = 0.1*ones(1,num_phases); % "Phase final times follow phase initial times"
    bounds.eventgroup(3).upper = 10*tfmax*ones(1,num_phases);
elseif Target_Ellipse % Have two phases with an added dynamical event constraint at approach terminus  
    bounds.eventgroup(1).lower = zeros(1,state_dim+1); % "Link states and times from phase 1 -> phase 2"
    bounds.eventgroup(1).upper = zeros(1,state_dim+1);
    if Nonlinear && coast  % Keplerian rel. dynamics
        bounds.eventgroup(2).lower = [0 0 0 0 0]; % "Satisfy dynamics to begin ellipse at approach terminus"
        bounds.eventgroup(2).upper = [0 0 0 0 0];
        bounds.eventgroup(3).lower = 0.1*ones(1,num_phases); % "Phase final times follow phase initial times"
        bounds.eventgroup(3).upper = 100000*ones(1,num_phases);
    elseif ~Nonlinear && coast % C.W. rel. dynamics
        bounds.eventgroup(2).lower = [0 0 0 0 0]; % "Satisfy dynamics to begin ellipse at approach terminus"
        bounds.eventgroup(2).upper = [0 0 0 0 0];
        bounds.eventgroup(3).lower = 0.1*ones(1,num_phases); % "Phase final times follow phase initial times"
        bounds.eventgroup(3).upper = 100000*ones(1,num_phases);
    elseif ~coast
        bounds.eventgroup(2).lower = 0.1*ones(1,num_phases); % "Phase final times follow phase initial times"
        bounds.eventgroup(2).upper = 100000*ones(1,num_phases);
    end%if
elseif Period_Match 
        bounds.eventgroup(1).lower = zeros(1,state_dim+1); % "Link states and times from phase 1 -> phase 2"
        bounds.eventgroup(1).upper = zeros(1,state_dim+1);
        bounds.eventgroup(2).lower = 0; % "Satisfy dynamics to null secular dynamics"
        bounds.eventgroup(2).upper = 0;
        bounds.eventgroup(3).lower = 0.1*ones(1,num_phases); % "Phase final times follow phase initial times"
        bounds.eventgroup(3).upper = 100000*ones(1,num_phases);
end%if

%-------------------------------------------------------------------------%
%----------Provide Mesh Refinement Method and Initial Mesh ---------------%
%-------------------------------------------------------------------------%
mesh.method       = 'hp-PattersonRao';
mesh.tolerance    = MeshTolerance;
mesh.maxiterations = MaxMeshIters;
mesh.colpointsmin = ColPointsMin;
mesh.colpointsmax = ColPointsMax;

%-------------------------------------------------------------------------%
%------------- Assemble Information into Problem Structure ---------------%        
%-------------------------------------------------------------------------%
setup.name                               = 'R_POTATO_Main';
setup.functions.continuous               = @R_POTATO_Continuous;
setup.functions.endpoint                 = @R_POTATO_Endpoint;
setup.auxdata                            = auxdata;
setup.bounds                             = bounds;
setup.guess                              = guess;
setup.mesh                               = mesh; 
if AutomaticBounds
    setup.scales.method                      = 'automatic-bounds';
else
    setup.scales.method                      = 'none';
end%if
setup.nlp.solver                         = 'ipopt';
setup.nlp.ipoptoptions.maxiterations     = MaxIPOPTiters;
setup.derivatives.supplier               = 'sparseCD';
setup.derivatives.derivativelevel        = 'second';
setup.method                             = 'RPM-Differentiation';
setup.displaylevel                       = 2;
%-------------------------------------------------------------------------%
%------------------------- Solve Problem Using GPOP2 ---------------------%
%-------------------------------------------------------------------------%
output   = gpops2(setup);
solution = output.result.solution;
