%-------------------------------------------%
% BEGIN: function R_POTATO_Endpoint.m    %
%-------------------------------------------%
function output = R_POTATO_Endpoint(input)
 
% Obtain auxiliary data
mu           = input.auxdata.mu;
num_phases   = input.auxdata.num_phases;
dim          = input.auxdata.dim;
ellipse      = input.auxdata.ellipse;
n            = input.auxdata.n;
Fmax         = input.auxdata.Fmax;
m            = input.auxdata.m;
coast        = input.auxdata.coast;
vbar_offset  = input.auxdata.vbar_offset;
alpha        = input.auxdata.alpha;
zeta         = input.auxdata.zeta;
gamma        = input.auxdata.gamma;
Nonlinear    = input.auxdata.Nonlinear;
e            = input.auxdata.e;  
P_Match      = input.auxdata.P_Match;
T_Point      = input.auxdata.T_Point;
de_ON        = input.auxdata.de_ON;
x_size       = input.auxdata.x_size;
dw_ON        = input.auxdata.dw_ON;
dw           = input.auxdata.dw;
di_ON        = input.auxdata.di_ON;
di           = input.auxdata.di;
dG_ON        = input.auxdata.dG_ON;
dG           = input.auxdata.dG;
enforceCW    = input.auxdata.enforceCW;
mnvr_t0      = input.auxdata.mnvr_epoch;
if e~=0
    nu_table = input.auxdata.nu_table;
end%if
b = Fmax/m;



    if num_phases > 1 && T_Point % Min-Fuel Rendezvous: Bang-off-bang
        
        % Variables at Start and Terminus of Phase 1
        t0{1} = input.phase(1).initialtime;
        tf{1} = input.phase(1).finaltime;
        s0{1} = input.phase(1).initialstate;
        sf{1} = input.phase(1).finalstate;
        
        % Variables at Start and Terminus of Phase 2
        t0{2} = input.phase(2).initialtime;
        tf{2} = input.phase(2).finaltime;
        s0{2} = input.phase(2).initialstate;
        sf{2} = input.phase(2).finalstate;
        
        % Variables at Start and Terminus of Phase 3
        t0{3} = input.phase(3).initialtime;
        tf{3} = input.phase(3).finaltime;
        s0{3} = input.phase(3).initialstate;
        sf{3} = input.phase(3).finalstate;        

        % Event Group 1: Linkage Constraints Between Phases 1 and 2
        output.eventgroup(1).event = [s0{2}(1:dim)-sf{1}(1:dim),t0{2}-tf{1}];
        
        % Event Group 2: Linkage Constraints Between Phases 2 and 3
        output.eventgroup(2).event = [s0{3}(1:dim)-sf{2}(1:dim),t0{3}-tf{2}];
        
        % Event Group 3:  Final Time of Each Phase Larger Than Initial Time of Phase
        output.eventgroup(3).event = [tf{1}-t0{1}, tf{2}-t0{2}, tf{3}-t0{3}];

        
        integral = input.phase(1).integral + ...
                   input.phase(2).integral + ...
                   input.phase(3).integral; 
        
    elseif ~Nonlinear && ~T_Point % C.W. Ellipse
        if ellipse % Design periodic trajectory
            switch coast
                case 1 % Entrance to ellipse is an endpoint constraint
                    xf                  = input.phase(1).finalstate(:,1);
                    yf                  = input.phase(1).finalstate(:,2);
                    zf                  = input.phase(1).finalstate(:,3);
                    xdotf               = input.phase(1).finalstate(:,4);
                    ydotf               = input.phase(1).finalstate(:,5);
                    zdotf               = input.phase(1).finalstate(:,6);

                    % Determine the coefficients of combined sinusoid
                    % solution for periodic trajectory in each coordinate
                    ax = -3.*xf - 2/n.*ydotf;
                    bx = 1/n.*xdotf;
                    Rx = sqrt(ax.^2 + bx.^2);

                    ay = 2/n.*xdotf;
                    by = 6.*xf + 4/n.*ydotf;
                    Ry = sqrt(ay.^2 + by.^2);

                    az = zf;
                    bz = 1/n.*zdotf;
                    Rz = sqrt(az.^2 + bz.^2);
                    
                    if zeta ~= 0
                        % Try to simplify trigonometry for 'easy' cases
                        if gamma == 0 
                            tau_y = atan2(by,ay);
                            tau_z = atan2(bz,az); % Z and Y dynamics in phase;
                            R_below_apex = Ry;
                            phase_error = tau_y - tau_z;
                        elseif gamma == 90 
                            tau_x = atan2(bx,ax);
                            tau_z = atan2(bz,az); % Z and X in phase;
                            R_below_apex = Rx;
                            phase_error = tau_x - tau_z;
                        elseif gamma == -90 
                            tau_x = atan2(bx,ax);
                            tau_z = atan2(-bz,-az); % Z and X are 180 deg out of phase;
                            R_below_apex = Rx;
                            phase_error = tau_x - tau_z;
                        elseif abs(gamma) == 180
                            tau_y = atan2(by,ay);
                            tau_z = atan2(-bz,-az); % Z and Y are 180 deg out of phase;
                            R_below_apex = Ry;
                            phase_error = tau_y - tau_z;
                        else % No 'easy' case specified for azimuth calcs
                            r = sqrt((alpha^2) / (5 - cosd(gamma)^2 - 4*sind(gamma)^2));
                            % What are the X and Y components of the local ascending node vector?
                            x_asc_node = r*cosd(180 - gamma);
                            y_asc_node = r*sind(180 - gamma); 
                            for i = 1:numel(Ry)
                                if y_asc_node/Ry(i) > 1
                                    y_phase(i) = 0;
                                elseif y_asc_node/Ry(i) < -1
                                    y_phase(i) = pi;
                                else
                                    y_phase(i) = acosd(y_asc_node/Ry);
                                end%if
                            end%for
                            % What is the magnitude of the orbit-plane projection of the range vector at max cross-track distance?
                            R_below_apex = sqrt((Rx.*cosd(y_phase))^2 + (Ry.*cosd(90+y_phase))^2);
                            tau_y = atan2(by,ay); 
                            tau_z = atan2(bz,az) - deg2rad(90 - y_phase);
                            % Specify azimuth as the phase difference between Y and Z dynamics
                            phase_error = tau_y - tau_z;
                        end%if
                        % Inclination is the interior angle between the range vector at max cross-track and its projection on orbit plane
                        incl_error = R_below_apex.*tand(zeta) - Rz;
                    elseif zeta == 0 % Try to simplify for in-plane case
                        phase_error = zeros(size(xf));
                        incl_error = abs(zf) + abs(zdotf);
                    end%if

                        % Variables at Start and Terminus of Phase 1
                        t0{1} = input.phase(1).initialtime;
                        tf{1} = input.phase(1).finaltime;
                        s0{1} = input.phase(1).initialstate;
                        sf{1} = input.phase(1).finalstate;

                        % Variables at Start and Terminus of Phase 2
                        t0{2} = input.phase(2).initialtime;
                        tf{2} = input.phase(2).finaltime;
                        s0{2} = input.phase(2).initialstate;
                        sf{2} = input.phase(2).finalstate;

                        % Event Group 1: Linkage Constraints Between Phases 1 and 2
                        output.eventgroup(1).event = [s0{2}(1:dim)-sf{1}(1:dim),t0{2}-tf{1}];

                        % Event Group 2: Dynamics that Specify Free Orbit Ellipse         
                        output.eventgroup(2).event = [ydotf + 2*n.*xf, ...                            % Null secular dynamics
                                                      ydotf.*(yf-vbar_offset)./(xdotf.*xf) + 4, ...   % Location
                                                      4.*xf.^2 + (yf-vbar_offset).^2 - alpha^2, ...   % Size              
                                                      incl_error, ...                                 % Inclination   
                                                      phase_error];                                   % Azimuth
                                                                                 

                        % Event Group 3:  Final Time of Each Phase Larger Than Initial Time of Phase
                        output.eventgroup(3).event = [tf{1}-t0{1}, tf{2}-t0{2}];

%                         % Pack objective values for all phases
%                         integral = input.phase(1).integral + ...
%                                    input.phase(2).integral;
                        % If coasting, there's nothing you can do to
                        % minimize cost function on second phase, so don't
                        % let GPOPS-II work with it
                        integral = input.phase(1).integral;

                       
                case 0 % Ellipse is enforced as a path constraint (Just need phase linking here)
                    % Variables at Start and Terminus of Phase 1
                    t0{1} = input.phase(1).initialtime;
                    tf{1} = input.phase(1).finaltime;
                    s0{1} = input.phase(1).initialstate;
                    sf{1} = input.phase(1).finalstate;

                    % Variables at Start and Terminus of Phase 2
                    t0{2} = input.phase(2).initialtime;
                    tf{2} = input.phase(2).finaltime;
                    s0{2} = input.phase(2).initialstate;
                    sf{2} = input.phase(2).finalstate;
                    
                    % Event Group 1: Linkage Constraints Between Phases 1 and 2
                    output.eventgroup(1).event = [s0{2}(1:dim)-sf{1}(1:dim),t0{2}-tf{1}];
                    % Event Group 2:  Final Time of Each Phase Larger Than Initial Time of Phase
                    output.eventgroup(2).event = [tf{1}-t0{1}, tf{2}-t0{2}];

                    integral = input.phase(1).integral + ...
                               input.phase(2).integral;
            end%switchcase
        elseif P_Match 
            xf                  = input.phase(2).finalstate(:,1);
            ydotf               = input.phase(2).finalstate(:,5);
            
            % Variables at Start and Terminus of Phase 1
            t0{1} = input.phase(1).initialtime;
            tf{1} = input.phase(1).finaltime;
            s0{1} = input.phase(1).initialstate;
            sf{1} = input.phase(1).finalstate;

            % Variables at Start and Terminus of Phase 2
            t0{2} = input.phase(2).initialtime;
            tf{2} = input.phase(2).finaltime;
            s0{2} = input.phase(2).initialstate;
            sf{2} = input.phase(2).finalstate;

            % Event Group 1: Linkage Constraints Between Phases 1 and 2
            output.eventgroup(1).event = [s0{2}(1:dim)-sf{1}(1:dim),t0{2}-tf{1}];

            % Event Group 2: Condition for Entering Periodic Trajectory
            output.eventgroup(2).event = ydotf + 2*n.*xf;

            % Event Group 3:  Final Time of Each Phase Larger Than Initial Time of Phase
            output.eventgroup(3).event = [tf{1}-t0{1}, tf{2}-t0{2}];

%             integral = input.phase(1).integral + ...
%                        input.phase(2).integral;
            % If coasting, there's nothing you can do to
            % minimize cost function on second phase, so don't
            % let GPOPS-II work with it
            integral = input.phase(1).integral;

        end%if
    elseif Nonlinear && ~T_Point % Keplerian relative dynamics
        if ellipse % Design NLEM periodic trajectory
            switch coast
                case 1 % Entrance to ellipse is an endpoint constraint
                    xf                  = input.phase(1).finalstate(:,1);
                    yf                  = input.phase(1).finalstate(:,2);
                    zf                  = input.phase(1).finalstate(:,3);
                    xdotf               = input.phase(1).finalstate(:,4);
                    ydotf               = input.phase(1).finalstate(:,5);
                    zdotf               = input.phase(1).finalstate(:,6);

                    % Variables at Start and Terminus of Phase 1
                    t0{1} = input.phase(1).initialtime;
                    tf{1} = input.phase(1).finaltime;
                    s0{1} = input.phase(1).initialstate;
                    sf{1} = input.phase(1).finalstate;

                    % Variables at Start and Terminus of Phase 2
                    t0{2} = input.phase(2).initialtime;
                    tf{2} = input.phase(2).finaltime;
                    s0{2} = input.phase(2).initialstate;
                    sf{2} = input.phase(2).finalstate;
                    
                    
                    a = (mu/(n^2))^(1/3);
                    if e > 0 % Time-varying dynamics - lookup true anomaly 
                        nu_interp = interp1qr(nu_table(:,1),nu_table(:,2),mod(tf{1}+mnvr_t0,2*pi/n));
                        r_tgt = (a*(1-e^2))./(1+e.*cos(nu_interp));
                        h = sqrt(mu .* a*(1-e^2));
                        nu_dot = h./(r_tgt.^2);
                        rd_tgt = (r_tgt.*nu_dot.*e.*sin(nu_interp))./(1+e.*cos(nu_interp));
                    elseif e == 0 % Time-invariant dynamics, true anomaly rate of change is constant. 
                        r_tgt = a;   
                        rd_tgt = 0;
                        h = sqrt(mu * a);
                        nu_dot = n;
                        rd_tgt = 0;
                    end%if

                    % Determine R and V for chaser in LVLH frame.
                    r_ch_x = r_tgt + xf;
                    r_ch_y = yf;
                    r_ch_z = zf;
                    r_ch  = sqrt(r_ch_x.^2 + r_ch_y.^2 + r_ch_z.^2);
                    v_ch_x = rd_tgt + xdotf - nu_dot.*yf;
                    v_ch_y = ydotf + nu_dot.*(xf+r_tgt);
                    v_ch_z = zdotf;
                    v_ch = sqrt( v_ch_x.^2 + v_ch_y.^2 + v_ch_z.^2 );
                    
                    % Determine specific mechanical energies based on
                    % semi-major axis (target) and vis-viva equation (chaser).
                    E_tgt = -mu/(2*a);
                    E_ch  = 0.5.*v_ch.^2 - mu./r_ch;

                    % Obtain the dot product of R and V for chaser 
                    dot_rv = (r_tgt+xf).*(xdotf-nu_dot.*yf+rd_tgt) + yf.*(ydotf+nu_dot.*(xf+r_tgt)) + zf.*zdotf;

                    % Obtain the target's perigee vector in LVLH frame
                    e_tgt_x = 1/mu*(nu_dot.^2*r_tgt.^3 - mu);
                    e_tgt_y = 1/mu*(-nu_dot.*rd_tgt.*r_tgt.^2);
                    e_tgt_z = 0;
                    
                    % Determine the chaser's eccentricity (perigee) vector in LVLH frame
                    e_ch_x = 1./mu.*((v_ch.^2 - mu./r_ch).*(r_tgt+xf) - dot_rv.*(xdotf-nu_dot.*yf+rd_tgt));
                    e_ch_y = 1./mu.*((v_ch.^2 - mu./r_ch).*yf         - dot_rv.*(ydotf + nu_dot.*(xf+r_tgt)));
                    e_ch_z = 1./mu.*((v_ch.^2 - mu./r_ch).*zf         - dot_rv.*zdotf);
                    % Eccentricity from magnitude of eccentricity vector
                    e_ch = sqrt(e_ch_x.^2 + e_ch_y.^2 + e_ch_z.^2);
                    
                    % Determine the differential eccentricity vector and
                    % magnitude
                    de_x = e_ch_x - e_tgt_x;
                    de_y = e_ch_y - e_tgt_y;
                    de_z = e_ch_z - e_tgt_z;
                    de = e_ch - e;
                    
                    % Determine the cosine of the angle between the
                    % eccentricity vectors
                    if e ~= 0 
                        cos_arg_peri = (e_ch_x.*e_tgt_x + e_ch_y.*e_tgt_y)./(sqrt(e_ch_x.^2+e_ch_y.^2).*sqrt(e_tgt_x.^2+e_tgt_y.^2));
%                        cos_arg_peri = ((nu_dot.^2.*r_tgt.^3-mu).*e_ch_x + (-nu_dot.*rd_tgt.*r_tgt.^2).*e_ch_y)./(mu.^2.*e.*e_ch);
                    elseif e == 0 % Perigee is undefined for target, so find the angle between true anomalies instead
                        a_ch = -mu./(2*E_ch);
                        nu_ch = acos((a_ch.*(1-e_ch^2) - 1)./e_ch);
                        cos_arg_peri = abs(nu_ch - n.*mod(tf{1},2*pi/n));
                        %             cos_arg_peri = ((nu_dot.^2.*r_tgt.^3-mu).*e_ch_x + (-nu_dot.*rd_tgt.*r_tgt.^2).*e_ch_y)./(mu.^2.*e.*e_ch);
                    end%if
                        
                    % Compute the chaser's angular momentum vector in the
                    % LVLH frame
                    h_ch_x = yf*zdotf - zf*(ydotf+nu_dot*(xdotf+r_tgt));
                    h_ch_y = -zdotf*(r_tgt+xf) + zf*(xdotf+rd_tgt-nu_dot*(yf));
                    h_ch_z = (r_tgt+xf)*(ydotf+nu_dot*(xf+r_tgt)) - yf*(xdotf+rd_tgt-nu_dot*yf);
                    h_ch = sqrt(mu * a*(1-e_ch^2));
                    
                    % Determine the cosine of the angle between the local angular
                    % momentum vectors (differential inclination)
                    cos_inc = h_ch_z / h_ch;

                    % Determine the X and Y components of the cross product
                    % between the target's and chaser's angular momentum
                    % vectors - this is the local ascending node
                    n_loc_x = -h*h_ch_y;
                    n_loc_y = sign(de)*h*h_ch_x;
%                     cos_RAAN = h_ch_x / h_ch; % sqrt(h_ch_x^2 + h_ch_y^2);

                    % Determine the cosine of the angle between the local
                    % ascending node vector and the differential
                    % eccentricity vector
                    cos_Gamma = (n_loc_x*de_x + n_loc_y*de_y)./(sqrt(n_loc_x^2 + n_loc_y^2)*sqrt(de_x^2+de_y^2));

                    % At the very least, chaser & target must have matching
                    % energy during the periodic phase
                    energy_match      = E_tgt - E_ch;
                    
                    % Compute the desired differential eccentricity based
                    % on the specified size
                    delta_ecc_reqd    = x_size/a;
                    
                    % Determine how parallel the local vectors are
                    parallel_arg_peri = cos_arg_peri - 1;
                    if di == 0 % Try to simplify the case where an in-plane ellipse is enforced
                        parallel_ang_mom  = abs(zf) + abs(zdotf); 
                    else 
                        parallel_ang_mom  = cos_inc - 1;
                    end%if
                    parallel_LAN     = cos_Gamma - 1;
                    
                    % Event Group 1: Linkage Constraints Between Phases 1 and 2
                    output.eventgroup(1).event = [s0{2}(1:dim)-sf{1}(1:dim),t0{2}-tf{1}];

                    % Event Group 2: Condition for Entering Ellipse Phase
                    % Note: scale all conditions by 1e6 to accentuate
                    % minima -- differential COEs are very "flat" at close
                    % range. Also, only enforce the constraints that are
                    % switched "ON" in the Instantiation Shell. 
                    output.eventgroup(2).event = 1e4.*[energy_match, ...                        % Period-matching condition
                                                       1e6*de_ON.*(de - delta_ecc_reqd), ...        % Size 
                                                       1e2.*dw_ON.*(parallel_arg_peri - dw), ...     % Location 
                                                       di_ON.*(parallel_ang_mom - di), ...      % Inclination 
                                                       dG_ON.*(parallel_LAN - dG)];             % Azimuth

                    % Event Group 3:  Final Time of Each Phase Larger Than Initial Time of Phase
                    output.eventgroup(3).event = [tf{1}-t0{1}, tf{2}-t0{2}];

%                     % Pack the objective values from all phases
%                     integral = input.phase(1).integral + ...
%                                input.phase(2).integral;
                    % If coasting, there's nothing you can do to
                    % minimize cost function on second phase, so don't
                    % let GPOPS-II work with it
                    integral = input.phase(1).integral;
                           
                case 0 % Free-orbit motion is enforced as a path constraint
                    xf                  = input.phase(1).finalstate(:,1);
                    yf                  = input.phase(1).finalstate(:,2);
                    zf                  = input.phase(1).finalstate(:,3);
                    xdotf               = input.phase(1).finalstate(:,4);
                    ydotf               = input.phase(1).finalstate(:,5);
                    zdotf               = input.phase(1).finalstate(:,6);

                    % Variables at Start and Terminus of Phase 1
                    t0{1} = input.phase(1).initialtime;
                    tf{1} = input.phase(1).finaltime;
                    s0{1} = input.phase(1).initialstate;
                    sf{1} = input.phase(1).finalstate;

                    % Variables at Start and Terminus of Phase 2
                    t0{2} = input.phase(2).initialtime;
                    tf{2} = input.phase(2).finaltime;
                    s0{2} = input.phase(2).initialstate;
                    sf{2} = input.phase(2).finalstate;
                    % Event Group 1: Linkage Constraints Between Phases 1 and 2
                    output.eventgroup(1).event = [s0{2}(1:dim)-sf{1}(1:dim),t0{2}-tf{1}];

        %             output.eventgroup(2).event = [xf-100, yf, xdotf];
        %             % Event Group 2:  Final Time of Each Phase Larger Than Initial Time of Phase
        %             output.eventgroup(3).event = [tf{1}-t0{1}, tf{2}-t0{2}];

                    output.eventgroup(2).event = [tf{1}-t0{1}, tf{2}-t0{2}];

                    integral = input.phase(1).integral + ...
                                input.phase(2).integral;
            end%switchcase
        elseif P_Match
            xf                  = input.phase(1).finalstate(:,1);
            yf                  = input.phase(1).finalstate(:,2);
            zf                  = input.phase(1).finalstate(:,3);
            xdotf               = input.phase(1).finalstate(:,4);
            ydotf               = input.phase(1).finalstate(:,5);
            zdotf               = input.phase(1).finalstate(:,6);

            % Variables at Start and Terminus of Phase 1
            t0{1} = input.phase(1).initialtime;
            tf{1} = input.phase(1).finaltime;
            s0{1} = input.phase(1).initialstate;
            sf{1} = input.phase(1).finalstate;

            % Variables at Start and Terminus of Phase 2
            t0{2} = input.phase(2).initialtime;
            tf{2} = input.phase(2).finaltime;
            s0{2} = input.phase(2).initialstate;
            sf{2} = input.phase(2).finalstate;
            
            a = (mu/(n^2))^(1/3);
            if e ~= 0 
                nu_interp = interp1qr(nu_table(:,1),nu_table(:,2),mod(tf{1}+mnvr_t0,2*pi/n));
                r_tgt = (a*(1-e^2))./(1+e.*cos(nu_interp));
                h = sqrt(mu .* a*(1-e^2));
                nu_dot = h./(r_tgt.^2);
                rd_tgt = (r_tgt.*nu_dot.*e.*sin(nu_interp))./(1+e.*cos(nu_interp));
            else
                r_tgt = a;   
                rd_tgt = 0;
                h = sqrt(mu * a);
                nu_dot = n;
                rd_tgt = 0;
            end%if
            
            r_ch  = sqrt((r_tgt + xf).^2 + yf.^2 + zf.^2);
            v_ch = sqrt( (rd_tgt + xdotf - nu_dot.*yf).^2 + (ydotf + nu_dot.*(xf+r_tgt)).^2 + zdotf.^2 );
            E_tgt = -mu/(2*a);
            E_ch  = 0.5.*v_ch.^2 - mu./r_ch;
            
            % Event Group 1: Linkage Constraints Between Phases 1 and 2
            output.eventgroup(1).event = [s0{2}(1:dim)-sf{1}(1:dim),t0{2}-tf{1}];
            
            % Event Group 1: Condition for end of phase 1
            output.eventgroup(2).event = 1e6*(E_tgt - E_ch); % Period-matching condition only
            
            % Event Group 3:  Final Time of Each Phase Larger Than Initial Time of Phase
            output.eventgroup(3).event = [tf{1}-t0{1}, tf{2}-t0{2}];
            
%             integral = input.phase(1).integral + ...
%                        input.phase(2).integral;

            % If coasting, there's nothing you can do to
            % minimize cost function on second phase, so don't
            % let GPOPS-II work with it
            integral = input.phase(1).integral;
        else
            integral = input.phase(1).integral; 
        end%if
        
    else
        integral = input.phase(1).integral; 
        
    end%if


output.objective = integral;



%-------------------------------------------%
% END: function R_POTATO_Endpoint.m      %
%-------------------------------------------%

