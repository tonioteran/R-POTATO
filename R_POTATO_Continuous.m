%---------------------------------------------%
% BEGIN: function R_POTATO_Continuous.m    %
%---------------------------------------------%
function phaseout = R_POTATO_Continuous(input)


n                  = input.auxdata.n;
mu                 = input.auxdata.mu;              
m                  = input.auxdata.m;
dim                = input.auxdata.dim;
mintime            = input.auxdata.mintime;
minfuel            = input.auxdata.minfuel;
minenergy          = input.auxdata.minenergy;
quadcost           = input.auxdata.quadcost;
qr_ratio           = input.auxdata.qr_ratio;
Fmax               = input.auxdata.Fmax;
Tmax               = input.auxdata.Tmax;
I                  = input.auxdata.I;
Thrust_dir         = input.auxdata.Thrust_dir;
num_phases         = input.auxdata.num_phases;
Nonlinear          = input.auxdata.Nonlinear;
ellipse            = input.auxdata.ellipse;
alpha              = input.auxdata.alpha;
zeta               = input.auxdata.zeta;
gamma              = input.auxdata.gamma;
e                  = input.auxdata.e;
avoid              = input.auxdata.avoid;
FOV                = input.auxdata.FOV;
buffer             = input.auxdata.FOV_buff;
if e~= 0 
    nu_table           = input.auxdata.nu_table;
end%if
vbar_offset        = input.auxdata.vbar_offset;
coast              = input.auxdata.coast;
P_Match            = input.auxdata.P_Match;
activate_range_constr = input.auxdata.range_path;
range              = input.auxdata.range;
activate_per_constr = input.auxdata.period_path;
activate_center_constr = input.auxdata.center_path;
center             = input.auxdata.center;
enforceCW          = input.auxdata.enforceCW;
CW_alpha           = input.auxdata.CW_alpha;
CW_zeta            = input.auxdata.CW_zeta;
t0                 = input.auxdata.mnvr_epoch;

b                  = Fmax/m;
I_inv              = inv(I);


for p = 1:num_phases
    if dim == 13
        s                  = input.phase(p).state;
            x                  = input.phase(p).state(:,1);
            y                  = input.phase(p).state(:,2);
            z                  = input.phase(p).state(:,3);
            xdot               = input.phase(p).state(:,4);
            ydot               = input.phase(p).state(:,5);
            zdot               = input.phase(p).state(:,6);
            q1                 = input.phase(p).state(:,7);
            q2                 = input.phase(p).state(:,8);
            q3                 = input.phase(p).state(:,9);
            q4                 = input.phase(p).state(:,10);
            wx                 = input.phase(p).state(:,11);
            wy                 = input.phase(p).state(:,12);
            wz                 = input.phase(p).state(:,13);
        uF                = input.phase(p).control(:,1);
        uTx               = input.phase(p).control(:,2);
        uTy               = input.phase(p).control(:,3);
        uTz               = input.phase(p).control(:,4);
%         uFx                = input.phase(p).control(:,1);
%         uFy                = input.phase(p).control(:,2);
%         uFz                = input.phase(p).control(:,3);
%         uTx               = input.phase(p).control(:,4);
%         uTy               = input.phase(p).control(:,5);
%         uTz               = input.phase(p).control(:,6);
    
% -------------------- INTERMEDIATE CALCULATIONS ------------------------- 
% --------------------  FOR ATTITUDE FORMULATION -------------------------


        % Rotation Matrix (DCM) that resolves body vector into LVLH frame
        R1 = q4.^2 + q1.^2 - q2.^2 - q3.^2;
        R2 = 2.*(q1.*q2 - q3.*q4);
        R3 = 2.*(q1.*q3 + q2.*q4);
        R4 = 2.*(q1.*q2 + q3.*q4);
        R5 = q4.^2 - q1.^2 + q2.^2 - q3.^2;
        R6 = 2.*(q2.*q3 - q1.*q4);
        R7 = 2.*(q1.*q3 - q2.*q4);
        R8 = 2.*(q2.*q3 + q1.*q4);
        R9 = q4.^2 - q1.^2 - q2.^2 + q3.^2;
        
        
        % Resolve the thrust direction vector into the LVLH frame
        
        Nr = R1*Thrust_dir(1) + R2*Thrust_dir(2) + R3*Thrust_dir(3);
        Ns = R4*Thrust_dir(1) + R5*Thrust_dir(2) + R6*Thrust_dir(3);
        Nw = R7*Thrust_dir(1) + R8*Thrust_dir(2) + R9*Thrust_dir(3);
        
%         Nx_r = R1;
%         Nx_s = R4;
%         Nx_w = R7;
%        
%         Ny_r = R2;
%         Ny_s = R5;
%         Ny_w = R8;
%         
%         Nz_r = R3;
%         Nz_s = R6;
%         Nz_w = R9;
        
        
        % Convert chaser body frame rates to LVLH frame (Boyarko et al.)
%         Wr = R1.*wx + R2.*wy + R3.*wz;
%         Ws = R4.*wx + R5.*wy + R6.*wz;
%         Ww = R7.*wx + R8.*wy + R9.*wz - n;
        
        % Quaternion Rate Determined with Quat-Multiplied Frame
        % Compositions
        Qdot1 = 0.5 * ( q2.*wz - q3.*wy + q4.*wx         + n.*q2);   
        Qdot2 = 0.5 * (-q1.*wz + q3.*wx + q4.*wy         - n.*q1);
        Qdot3 = 0.5 * ( q1.*wy - q2.*wx + q4.*wz         - n.*q4);
        Qdot4 = 0.5 * (-q1.*wx - q2.*wy - q3.*wz         + n.*q3);
    
    elseif dim == 6
        t                  = input.phase(p).time;
        s                  = input.phase(p).state;
            x                  = input.phase(p).state(:,1);
            y                  = input.phase(p).state(:,2);
            z                  = input.phase(p).state(:,3);
            xdot               = input.phase(p).state(:,4);
            ydot               = input.phase(p).state(:,5);
            zdot               = input.phase(p).state(:,6);
        uFx                = input.phase(p).control(:,1);
        uFy                = input.phase(p).control(:,2);
        uFz                = input.phase(p).control(:,3);
    elseif dim == 4
        t                  = input.phase(p).time;
        s                  = input.phase(p).state;
            x                  = input.phase(p).state(:,1);
            y                  = input.phase(p).state(:,2);
            xdot               = input.phase(p).state(:,3);
            ydot               = input.phase(p).state(:,4);
        uFx                = input.phase(p).control(:,1);
        uFy                = input.phase(p).control(:,2);
    elseif dim == 2
        t                  = input.phase(p).time;
        s                  = input.phase(p).state;
            z                  = input.phase(p).state(:,1);
            zdot               = input.phase(p).state(:,2);
        uFz                = input.phase(p).control(:,1);
    end%if



%----------------------------- DYNAMICS ---------------------------------
    
    if dim == 13 % If including attitude dynamics:
        
            % Clohessy-Wiltshire Equations of Motion (with forcing term)
            xddot              = 2*n.*ydot + 3*n^2.*x   + b*Nr.*uF; 
            yddot              = -2*n.*xdot             + b*Ns.*uF; 
            zddot              = -n^2.*z                + b*Nw.*uF;
%             xddot              = 2*n.*ydot + 3*n^2.*x   + b.*(Nx_r.*uFx + Ny_r.*uFy + Nz_r.*uFz); 
%             yddot              = -2*n.*xdot             + b.*(Nx_s.*uFx + Ny_s.*uFy + Nz_s.*uFz); 
%             zddot              = -n^2.*z                + b.*(Nx_w.*uFx + Ny_w.*uFy + Nz_w.*uFz);



        %     % Quaternion Kinematic Differential Equation 
    %         q1dot              = Q1 .*q1 + Q2 .*q2 + Q3 .*q3 + Q4 .*q4;
    %         q2dot              = Q5 .*q1 + Q6 .*q2 + Q7 .*q3 + Q8 .*q4;
    %         q3dot              = Q9 .*q1 + Q10.*q2 + Q11.*q3 + Q12.*q4;
    %         q4dot              = Q13.*q1 + Q14.*q2 + Q15.*q3 + Q16.*q4;
            % Quaternion Rate Determined with Quat-Multiplied Frame
            % Compositions
            q1dot              = Qdot1;
            q2dot              = Qdot2;
            q3dot              = Qdot3;
            q4dot              = Qdot4;
            % Euler's Moment Equations (with forcing term)
            % (NOTE: I(1) = Ixx, I(5) == Iyy, I(9) == Izz)
            wxdot              = (I(5) - I(9))/I(1) .* wy .* wz + (Tmax/I(1) .* uTx);
            wydot              = (I(9) - I(1))/I(5) .* wx .* wz + (Tmax/I(5) .* uTy);
            wzdot              = (I(1) - I(5))/I(9) .* wx .* wy + (Tmax/I(9) .* uTz);
       
        phaseout(p).dynamics  = [xdot, ydot, zdot, ...
                         xddot, yddot, zddot, ...
                         q1dot, q2dot, q3dot, q4dot, ...
                         wxdot, wydot, wzdot];
                     
    elseif (dim == 6 || dim == 4 || dim == 2) % Else, if not including attitude dynamics:
        
        if Nonlinear % Use Nonlinear Equations of Relative Motion
            a = (mu/(n^2))^(1/3);
            if e ~= 0 % Time-varying dynamics - lookup true anomaly 
                nu_interp = interp1qr(nu_table(:,1),nu_table(:,2),mod(t+t0,2*pi/n));
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
            
            r_ch  = sqrt((r_tgt + x).^2 + y.^2 + z.^2);
%             v_tgt = sqrt(2*mu./(r_tgt) - mu/a);
%             cos_phi = h./(r_tgt.*v_tgt);
%             for i = 1:numel(cos_phi)
%                 if cos_phi(i) > 1.0
%                     cos_phi(i) = 1.0;
%                 end%if
%             end%for
%             phi_tgt = acos(cos_phi);
%             v_ch = sqrt( (rd_tgt + xdot - nu_dot.*y).^2 + (ydot + nu_dot.*(x+r_tgt)).^2 + zdot.^2 );
%             E_tgt = -mu/(2*a);
%             E_ch  = 0.5.*v_ch.^2 - mu./r_ch;
            
%             dot_rv = (r_tgt+x).*(xdot-nu_dot.*y+rd_tgt) + y.*(ydot+nu_dot.*(x+r_tgt)) + z.*zdot;
% 
%             e_ch_x = (v_ch.^2 - mu./r_ch).*(r_tgt+x) - dot_rv.*(xdot-nu_dot.*y+rd_tgt);
%             e_ch_y = (v_ch.^2 - mu./r_ch).*y        - dot_rv.*(ydot + nu_dot.*(x+r_tgt));
%             e_ch_z = (v_ch.^2 - mu./r_ch).*z        - dot_rv.*zdot;

%             e_ch = 1./mu.*sqrt(e_ch_x.^2 + e_ch_y.^2 + e_ch_z.^2);

%             cos_arg_peri = ((nu_dot.^2.*r_tgt.^3-mu).*e_ch_x + (-nu_dot.*rd_tgt.*r_tgt.^2).*e_ch_y)./(mu.^2.*e.*e_ch);
            
            xddot = 2.*nu_dot.*ydot + nu_ddot.*y + nu_dot.^2.*x + mu./(r_tgt.^2) - mu./(r_ch.^3).*(r_tgt + x) + b.*uFx;
            yddot = -2.*nu_dot.*xdot - nu_ddot.*x + nu_dot.^2.*y - mu./(r_ch.^3).*y + b.*uFy;
            zddot = -mu./(r_ch.^3).*z + b.*uFz;
            phaseout(p).dynamics  = [xdot, ydot, zdot, ...
                                     xddot, yddot, zddot ];

            if ~coast
                if enforceCW
                    x_ellipse = CW_alpha/2.*sin(n.*(t+t0));
                    x_error = x - x_ellipse;
                    y_ellipse = CW_alpha.*cos(n.*(t+t0));
                    y_error = y - y_ellipse;
                    z_ellipse = CW_alpha/2*tand(CW_zeta).*sin(n.*(t+t0));
                    z_error = z - z_ellipse;
                end%if
            end%if
            
            
            
        else
            % Clohessy-Wiltshire Equations of Motion (with forcing term)
            switch dim % Clohessy Wiltshire dynamics are separable between in-plane and out-of-plane motion
                case 6
                    xddot              = 2*n.*ydot + 3*n^2.*x   + b.*uFx; 
                    yddot              = -2*n.*xdot             + b.*uFy; 
                    zddot              = -n^2.*z                + b.*uFz;
                    phaseout(p).dynamics  = [xdot, ydot, zdot, ...
                                             xddot, yddot, zddot ];
                case 4
                    xddot              = 2*n.*ydot + 3*n^2.*x   + b.*uFx; 
                    yddot              = -2*n.*xdot             + b.*uFy; 
                    phaseout(p).dynamics  = [xdot, ydot, ...
                                             xddot, yddot ];
                case 2
                    zddot              = -n^2.*z                + b.*uFz;
                    phaseout(p).dynamics  = [zdot, ...
                                             zddot ];
            end%switchcase
            
            
            if avoid
                avoidance_constr = x + sqrt(cotd((FOV+buffer)/2)^2.*(y.^2 + z.^2));
                phaseout(1).path = avoidance_constr;
            else
                avoidance_constr = [zeros(size(t))];
            end%if
                     
            if ellipse && ~coast
                  range_constr = activate_range_constr.*((sqrt(x.^2 + y.^2 + z.^2) - range))*1e4;
                  period_match_constr =  activate_per_constr.*(ydot+2*n.*x)*1e4 ;
                  center_constr = activate_center_constr.*((x.*xdot)./((y + center).*ydot) + 0.25).*1e4;
                  phaseout(2).path = [range_constr, period_match_constr, center_constr]; 
            end%if

            
        end%if
    end%if
                         
%     phaseout(p).path = ones(size(input.phase(p).time)) - quatnorm([q1 q2 q3 q4]);

% ------------------------------ COST FUNCTIONS ---------------------------                         
                         
    if dim == 13 % For coupled position + attitude dynamics (6-DOF problem)
        % For min thruster + torque ctrl energy:
        if minenergy
            phaseout(p).integrand = 0.5*uF.^2 + 0.5*(uTx).^2 + 0.5*(uTy).^2 + 0.5*(uTz).^2; % + ones(size(input.phase.time));
        end%if

        % For min fuel:
        if minfuel
            phaseout(p).integrand =  abs(uF); 
        end%if

        % For min time:
        if mintime
            if P_Match
                phaseout(p).integrand = ones(size(input.phase(p).time)) + 0.5.*(uTx.^2 + uTy.^2 + uTz.^2);
            else
                phaseout(p).integrand = ones(size(input.phase(p).time));
            end%if
        end%if
    else % For 3-DOF problem
        % For min thruster energy
        if minenergy
                switch dim
                    case 6
                        phaseout(p).integrand = 0.5.*(uFx.^2 + uFy.^2 + uFz.^2); 
                    case 4
                        phaseout(p).integrand = 0.5.*(uFx.^2 + uFy.^2); 
                    case 2
                        phaseout(p).integrand = 0.5.*(uFz.^2); 
                end%switchcase
            end%if
        end%if

        % For min fuel:
        if minfuel
            switch dim
                case 6
%                     phaseout(p).integrand = uFx + uFy + uFz; 
                    phaseout(p).integrand = abs(uFx) + abs(uFy) + abs(uFz); 
                case 4
%                     phaseout(p).integrand = uFx + uFy; 
                    phaseout(p).integrand = abs(uFx) + abs(uFy); 
                case 2
%                     phaseout(p).integrand = uFz; 
                    phaseout(p).integrand = abs(uFz); 
            end%switchcase
        end%if

        % For min time:
        if mintime
            phaseout(p).integrand = ones(size(input.phase(p).time));
        end%if
        
        % For LQR cost: 
        if quadcost
            switch dim
                case 6
                    phaseout(p).integrand = 0.1.*qr_ratio.*(x.^2 + y.^2 + z.^2 + xdot.^2 + ydot.^2 + zdot.^2) + ...
                                            0.1.*(uFx.^2 + uFy.^2 + uFz.^2); 
                case 4
                    phaseout(p).integrand = 0.1.*qr_ratio.*(x.^2 + y.^2 + xdot.^2 + ydot.^2) + ...
                                            0.1.*(uFx.^2 + uFy.^2); 
                case 2
                    phaseout(p).integrand = 0.1.*qr_ratio.*(z.^2 + zdot.^2) + ...
                                            0.1.*(uFz.^2); 
            end%switchcase
        end%if
        
        if Nonlinear && enforceCW % This overrides any other cost function specification. Enforce CW ellipse with LQR.
            phaseout(p).integrand = 0.5.*1e-6*(x_error.^2 + y_error.^2 + z_error.^2) + 0.5.*(uFx.^2 + uFy.^2 + uFz.^2); 
        end%if

    end%if
              
        

end%for

%---------------------------------------------%
% END: function R_POTATO_Continuous.m      %
%---------------------------------------------%
