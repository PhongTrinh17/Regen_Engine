% Test script for 1d regen engine

clc; clear; close all;
set(groot, 'defaultFigureColor', 'w');
set(groot, 'defaultAxesColor', 'w');
set(groot, 'defaultAxesXColor', 'k');
set(groot, 'defaultAxesYColor', 'k');
set(groot, 'defaultTextColor', 'k');

%% Basic Parameters

F = 1480 * 4.44822; % Target thrust in N
Pamb = 9.94; % psia At 10,500 ft altitude
T_amb = 298; % K (from feed calc sheet)
mfrac_eth = 0.75;
mfrac_h2o = 1 - mfrac_eth;
oxidizer = 'LOX';
% Assumed efficiencies
cstar_eff = 0.85;
cf_eff = 0.98;

% CEA parameters where Pamb = 9.94 psia (interpolate for diff values across nozzle after)
% Taken from throat (A/At = 1.00)
o_f = 1.3;
Pc_us = 300; % psia, target
Pc = convpres(Pc_us, 'psi', 'Pa');
card_str = sprintf(['fuel C2H5OH(L)   C 2 H 6 O 1\n', ...
    'h,cal=-66370.0      t(k)=298.00      wt%%=75.00\n', ...
    'fuel water H 2.0 O 1.0  wt%%=25.00\n', ...
    'h,cal=-68308.  t(k)=298.00 rho,g/cc = 0.9998']);

py.rocketcea.cea_obj.add_new_fuel('ETHANOL_WATER_75_25(L)', card_str);
fuel = 'ETHANOL_WATER_75_25(L)';
c = py.rocketcea.cea_obj.CEA_Obj(pyargs('oxName', oxidizer,'fuelName', fuel));
eps = c.get_eps_at_PcOvPe(pyargs('Pc', Pc_us, 'MR', o_f, 'PcOvPe',(Pc_us / Pamb)));
% transport properties [Cp, mu, k, Prandtl]
transport_chamber = c.get_Chamber_Transport(pyargs('Pc', Pc_us, 'MR', o_f, 'eps', eps));
transport_throat = c.get_Throat_Transport(pyargs('Pc', Pc_us, 'MR', o_f, 'eps', eps));
transport_exit = c.get_Exit_Transport(pyargs('Pc', Pc_us, 'MR', o_f, 'eps', eps));
cp_g_ref = [double(transport_chamber{1}), double(transport_throat{1}), double(transport_exit{1})] .* 4184; % J/kg*K
mu_g_ref = [double(transport_chamber{2}), double(transport_throat{2}), double(transport_exit{2})] .* 0.0001; % Pa*s viscosity
k_g_ref = [double(transport_chamber{3}), double(transport_throat{3}), double(transport_exit{3})] .* 0.4184; % W/m*K thermal conductivity
prandtl_ref = [double(transport_chamber{4}), double(transport_throat{4}), double(transport_exit{4})];

gamma_throat = c.get_Throat_MolWt_gamma(pyargs('Pc', Pc_us, 'MR', o_f, 'eps', eps));
gamma = double(gamma_throat{2}); % specific heat ratio, constant for isentropic relations

temps_cea = c.get_Temperatures(pyargs('Pc', Pc_us, 'MR', o_f, 'eps', eps)); % R
T_stag = double(temps_cea{1}) * 5/9; % K, Chamber (adiabatic Wall -> constant stagnation temp)

cstar_theo = double(c.get_Cstar(pyargs('Pc', Pc_us, 'MR', o_f))) * 0.3048; % m/s
cf_cea = c.get_PambCf(pyargs('Pamb', Pamb, 'Pc', Pc_us, 'MR', o_f, 'eps', eps));
cf_theo = double(cf_cea{1});
MW = 23.446; % g/mol

% Calculating target injector manifold pressure
fuel_stiffness = 0.20; % standard
dP_inj = fuel_stiffness * Pc; % Pa
mdot_total = F / (cstar_theo * cstar_eff * cf_theo * cf_eff); % kg/s
cstar_act = cstar_theo .* cstar_eff; % m/s
cstar_act_us = cstar_act * 39.3701; %in/s
mdot_f = mdot_total ./ (1 + o_f); % kg/s
rho_f_inj = 789 * mfrac_eth + 1000 * mfrac_h2o; % kg/m^3 (20C)
CdA_f_inj = mdot_f / sqrt(2 * rho_f_inj * dP_inj); % m^2 (Heritage)
P_target = Pc + dP_inj; % Pa manifold pressure (target value)

%% Throat Geometry (Rao)
At = (cstar_act * mdot_total)/Pc; % m^2
Rt = sqrt(At/pi); % m
R_curve = 1.5 * Rt;
R_curve_us = R_curve * 39.3701; % m - in

%% Diverging Geometry
%Me = sqrt((2/(gamma-1))*((Pc_us/Pamb)^((gamma - 1)/gamma) - 1));
%Ae = At * (1/Me)*((2/(gamma + 1))*(1 + ((gamma - 1)/2) * Me^2))^((gamma + 1)/(2*(gamma - 1)));
Ae = At * eps; % m^2
Re = sqrt(Ae/pi); % m
percent_len = 0.8; % input (0.8 is optimal fractional length for most cases)
len_div = percent_len * ((Re - Rt)/tan(deg2rad(15))); % m diverging length from radii
theta_e = deg2rad(13); % deg, from HH fig 4.16
theta_n = deg2rad(23); % deg

%% Chamber and Converging Geometry
conv_angle = deg2rad(45); % deg, standard
L_star = 40; % in, optimal for ethanol/lox
id_chamber = 4.75; % in, heritage
V_total = convlength(L_star, 'in', 'm') * At; % m^3
Rc = convlength(id_chamber / 2, 'in', 'm'); % m, chamber radius
x_conv_tangent = -(R_curve) * sin(conv_angle);
y_conv_tangent = Rt + (R_curve) * (1 - cos(conv_angle));
% 1.5*Rt converging fillet
R_fillet = Rt;
y_conv_fillet = Rc - R_fillet; % center height of fillet circle
y_f_tangent = y_conv_fillet + R_fillet*cos(conv_angle);
x_f_tangent = x_conv_tangent - (y_f_tangent - y_conv_tangent) / tan(conv_angle);
x_conv_fillet = x_f_tangent - R_fillet * sin(conv_angle); % x center of fillet circle
x_chamber_end = x_conv_fillet;

len_conv = -x_chamber_end; % m
V_conv = (1/3 * pi * len_conv) * (Rc^2 + Rt^2 + Rc*Rt); % m^3
V_chamber = V_total - V_conv;
len_chamber = V_chamber/(pi*Rc^2);
len_total = len_chamber + len_conv + len_div; % m

%% Encoding Geometry

dx = 0.001; % 1 mm step size

% Recentering coordinates to around throat
x_t = 0;
x_chamber_start = -(len_chamber + len_conv);
x_exit = len_div;
% Position arrrays
pos_i = x_chamber_start:dx:x_exit; % axial position
pos_j = zeros(size(pos_i)); % radii
% Diverging arc
xn = (0.382 * Rt) * sin(theta_n); % diverging tangent point where parabola starts (theta_n)
Rn = Rt + (0.382 * Rt) * (1 - cos(theta_n)); % y coordinate of xn
% Rao parabola (y = ax^2 + bx + c)
% Need to derive parabola that hits (xn, Rn), (xe, Re) w/ starting slope tan(theta_n)
matrix_A = [xn^2, xn, 1; x_exit^2, x_exit, 1; 2*xn, 1, 0];
matrix_B = [Rn; Re; tan(theta_n)];
coeffs = matrix_A \ matrix_B;
a_Rao = coeffs(1); b_Rao = coeffs(2); c_Rao = coeffs(3);

% Encoder loop
for k = 1:length(pos_i)
    x = pos_i(k);
    if x < x_chamber_end % chamber
        pos_j(k) = Rc;
    elseif x >= x_chamber_end && x < x_f_tangent
        pos_j(k) = y_conv_fillet + sqrt(R_fillet^2 - (x - x_conv_fillet)^2);
    elseif x >= x_f_tangent && x < x_conv_tangent % straight converging cone
        pos_j(k) = y_f_tangent - tan(conv_angle) * (x - x_f_tangent);
    elseif x >= x_conv_tangent && x < 0 % converging throat arc
        pos_j(k) = (Rt + R_curve) - sqrt((R_curve)^2 - x^2); % Rt + arc length - height of curve at point
    elseif x >= 0 && x < xn % diverging throat arc
        pos_j(k) = (Rt + 0.382*Rt) - sqrt((0.382*Rt)^2 - x^2);
    else % Rao parabola approx
        pos_j(k) = a_Rao * x^2 + b_Rao * x + c_Rao;
    end
end

%% Cooling Channel and Outer Jacket Geometry
% Fixed channel height and wall thickness
min_tol = 0.001; % m 3d printer tolerance
D_gas = 2 * pos_j; % Array of gas-side diameter at every node
D_t = 2 * Rt;
D_t_us = D_t * 39.3701; % m - in
h_channel = min_tol; % channel height (radial)
wall_thickness = min_tol; % HW/CW
w_rib = min_tol; % fixed, rib width
D_channel_base = D_gas + 2 * wall_thickness; % Engine diameters with added wall thickness
D_t_base = D_t + 2 * wall_thickness; % Throat diameter with added wall thickness
% Number of channels determined at throat
circ_t_base = pi * (D_t_base); % Gas side diameter + wall thickness
num_channel = floor(circ_t_base / (w_rib + min_tol));
% Local circumferences across engine
circ_local_base = pi * D_channel_base;
w_channel = (circ_local_base - (num_channel * w_rib)) ./ num_channel; % variable, channel widths
D_h = (4 .* w_channel .* h_channel) ./ (2 * w_channel + 2 * h_channel); % hydraulic diameter of rectangular channels
Per_heated = w_channel + 2 * h_channel; % heated perimeter

%% HT Areas
A_gas = pi .* D_gas .* dx; % Gas-wall convection SA
A_w = pi .* ((D_gas + D_channel_base)./2) .* dx; % wall-wall conduction SA (average diameter)
A_co = Per_heated .* num_channel .* dx; % coolant side surface area (heated)
A_wc = w_channel .* dx; % cool wall area per increment

%% Visualization Plot

figure('Name', '1D Engine Geometry', 'Color', 'w');
hold on; grid on;
plot(pos_i, pos_j, 'k', 'LineWidth', 2, 'DisplayName', 'Hot Wall');
plot(pos_i, -pos_j, 'k', 'LineWidth', 2, 'HandleVisibility','off');

r_channel_base = D_channel_base ./ 2;
plot(pos_i, r_channel_base, 'b--', 'LineWidth', 1.5, 'DisplayName', 'Cold Wall');
plot(pos_i, -r_channel_base, 'b--', 'LineWidth', 1.5, 'HandleVisibility', 'off');

r_outer_jacket = r_channel_base + h_channel;
plot(pos_i, r_outer_jacket, 'b', 'LineWidth', 2, 'DisplayName', 'Outer Jacket');
plot(pos_i, -r_outer_jacket, 'b', 'LineWidth', 2, 'HandleVisibility', 'off');

title('Regen 1D Profile')
xlabel('Axial Position x (m)');
ylabel('Radial Position y (m)');
axis equal;
xline(0, 'r--', 'Throat', 'LabelVerticalAlignment', 'bottom', 'HandleVisibility', 'off');
exportgraphics(gcf, 'geometry.pdf', 'ContentType','vector');
hold off;


%% Gas Properties
% 1D interpolation from 3 CEA points for Cp, gamma, k, mu for Bartz
M_local = zeros(size(pos_i));
A_local = D_gas.^2 .* (pi/4);
AR_local = A_local ./ At;
for k = 1:length(pos_i)
    if pos_i(k) < 0 % Chamber and Converging
        M_local(k) = flowisentropic(gamma, AR_local(k), 'sub');
    elseif pos_i(k) == 0
        M_local(k) = 1.0;
    else
        M_local(k) = flowisentropic(gamma, AR_local(k), 'sup');
    end
end
M_chamber =  flowisentropic(gamma, Rc^2/Rt^2, 'sub');
M_ref = [M_chamber, 1.0, M_local(end)];

cp_g_local = interp1(M_ref, cp_g_ref, M_local, 'linear', 'extrap');
mu_g_local = interp1(M_ref, mu_g_ref, M_local, 'linear', 'extrap');
k_g_local = interp1(M_ref, k_g_ref, M_local, 'linear', 'extrap');
prandtl_g_local = interp1(M_ref, prandtl_ref, M_local, 'linear', 'extrap');

Taw = T_stag * ((1 + prandtl_g_local.^(1/3).*((gamma - 1) / 2) .* M_local.^2) ...
    ./ (1 + ((gamma - 1) / 2) .* M_local.^2));
% Old code, doesn't work because interp1 need x values to be monotonically increasing, switch to using M_local
%{ 
AR_ref = [Rc^2/Rt^2, 1, eps];
A_local = D_gas.^2 .* (pi/4);
AR_local = A_local ./ At;
cp_g_local = interp1(AR_ref, cp_g_ref, AR_local, 'linear', 'extrap');
mu_g_local = interp1(AR_ref, mu_g_ref, AR_local, 'linear', 'extrap');
k_g_local = interp1(AR_ref, k_g_ref, AR_local, 'linear', 'extrap');
prandtl_g_local = interp1(AR_ref, prandtl_ref, AR_local, 'linear', 'extrap');
% Local Mach from isentropic area-mach relation
M_local = zeros(size(pos_i));
for k = 1:length(pos_i)
    if pos_i(k) < 0 % Chamber and Converging
        M_local(k) = flowisentropic(gamma, AR_local(k), 'sub');
    elseif pos_i(k) == 0
        M_local(k) = 1.0;
    else
        M_local(k) = flowisentropic(gamma, AR_local(k), 'sup');
    end
end
Taw = T_stag * ((1 + prandtl_g_local.^(1/3).*((gamma - 1) / 2) .* M_local.^2) ...
    ./ (1 + ((gamma - 1) / 2) .* M_local.^2));
%}

%% Coolant Properties

% Required properties update for equilibrium: Cp, k, mu, rho
% Define Table bounds (Tmax < T_sat at Pmin or coolprop will crash)
table_filename = 'coolprop_tables.mat';
if isfile(table_filename)
    load(table_filename);
else % create tables first time (delete file when changing parameters)
    res = 50; % 50x50 data grid
    P_min = convpres(300, 'psi', 'Pa'); % Pa, at chamber
    P_max = convpres(900, 'psi', 'Pa'); % Pa, at fuel tank
    T_min = 290; % K, standard inlet temp
    T_max = 450; % K, below 300 psi boiling point of ethanol
    % 1D vectors for P and T (fast solve for Coolprop, easy to get)
    P_vec = linspace(P_min, P_max, res);
    T_vec = linspace(T_min, T_max, res);
    % use ndgrid for higher dimensionality use
    [P_grid, T_grid] = ndgrid(P_vec, T_vec);
    % 2D grids for bulk properties
    rho_grid = zeros(res,res);
    cp_grid = zeros(res,res);
    mu_grid = zeros(res,res);
    k_grid = zeros(res,res);
    % 1D grids
    T_sat_grid = zeros(1,res);
    rho_l_grid = zeros(1,res); % sat. liquid density
    rho_v_grid = zeros(1,res); % sat. vapor density
    surften_grid = zeros(1,res); % surface tension
    h_fg_grid = zeros(1,res); % latent heat of vaporization

    % Generate values
    for i = 1:res
        P_val = P_vec(i);
        if P_val < 6140000 % Ethanol critical pressure in Pa
            T_sat_eth = py.CoolProp.CoolProp.PropsSI('T','P',P_val,'Q',0, 'ethanol');
            T_sat_h2o = py.CoolProp.CoolProp.PropsSI('T','P',P_val,'Q',0, 'water');
            T_sat_grid(i) = mfrac_eth * T_sat_eth + mfrac_h2o * T_sat_h2o;
            % Liquid Density
            rho_l_eth = py.CoolProp.CoolProp.PropsSI('D','P',P_val,'Q',0, 'ethanol');
            rho_l_h2o = py.CoolProp.CoolProp.PropsSI('D','P',P_val,'Q',0, 'water');
            rho_l_grid(i) = 1 / ((mfrac_eth / rho_l_eth) + (mfrac_h2o / rho_l_h2o));
            % Vapor Density
            rho_v_eth = py.CoolProp.CoolProp.PropsSI('D','P',P_val,'Q',1, 'ethanol');
            rho_v_h2o = py.CoolProp.CoolProp.PropsSI('D','P',P_val,'Q',1, 'water');
            rho_v_grid(i) = 1 / ((mfrac_eth / rho_v_eth) + (mfrac_h2o / rho_v_h2o));
            % Surface Tension (Interfacial)
            surften_eth = py.CoolProp.CoolProp.PropsSI('I','P',P_val,'Q',0, 'ethanol');
            surften_h2o = py.CoolProp.CoolProp.PropsSI('I','P',P_val,'Q',0, 'water');
            surften_grid(i) = mfrac_eth * surften_eth + mfrac_h2o * surften_h2o;
            % Latent Heat of Vaporization (H_vap - H_liq)
            h_v_eth = py.CoolProp.CoolProp.PropsSI('H','P',P_val,'Q',1, 'ethanol');
            h_l_eth = py.CoolProp.CoolProp.PropsSI('H','P',P_val,'Q',0, 'ethanol');
            h_fg_eth = h_v_eth - h_l_eth;
            h_v_h2o = py.CoolProp.CoolProp.PropsSI('H','P',P_val,'Q',1, 'water');
            h_l_h2o = py.CoolProp.CoolProp.PropsSI('H','P',P_val,'Q',0, 'water');
            h_fg_h2o = h_v_h2o - h_l_h2o;
            h_fg_grid(i) = mfrac_eth * h_fg_eth + mfrac_h2o * h_fg_h2o;
        else % undefined for superheated vapor
            T_sat_grid(i) = NaN;
            rho_l_grid(i) = NaN;
            rho_v_grid(i) = NaN;
            surften_grid(i) = NaN;
            h_fg_grid(i) = NaN;
        end
        for j = 1:res
            P_val = P_grid(i,j);
            T_val = T_grid(i,j);
            
            rho_eth = py.CoolProp.CoolProp.PropsSI('D','T',T_val,'P',P_val,'ethanol');
            rho_h2o = py.CoolProp.CoolProp.PropsSI('D','T',T_val,'P',P_val,'water');
            cp_eth = py.CoolProp.CoolProp.PropsSI('C','T',T_val,'P',P_val,'ethanol');
            cp_h2o = py.CoolProp.CoolProp.PropsSI('C','T',T_val,'P',P_val,'water');
            mu_eth = py.CoolProp.CoolProp.PropsSI('V','T',T_val,'P',P_val,'ethanol');
            mu_h2o = py.CoolProp.CoolProp.PropsSI('V','T',T_val,'P',P_val,'water');
            k_eth = py.CoolProp.CoolProp.PropsSI('L','T',T_val,'P',P_val,'ethanol');
            k_h2o = py.CoolProp.CoolProp.PropsSI('L','T',T_val,'P',P_val,'water');
            
            rho_grid(i,j) =  1 / ((mfrac_eth / rho_eth) + (mfrac_h2o / rho_h2o));
            cp_grid(i,j) = mfrac_eth * cp_eth + mfrac_h2o * cp_h2o;
            mu_grid(i,j) = mfrac_eth * mu_eth + mfrac_h2o * mu_h2o;
            k_grid(i,j) = mfrac_eth * k_eth + mfrac_h2o * k_h2o;
        end
    end
    % Interpolation objects
    get_rho = griddedInterpolant(P_grid, T_grid, rho_grid, 'linear', 'nearest');
    get_cp = griddedInterpolant(P_grid, T_grid, cp_grid, 'linear', 'nearest');
    get_mu = griddedInterpolant(P_grid, T_grid, mu_grid, 'linear', 'nearest');
    get_k = griddedInterpolant(P_grid, T_grid, k_grid, 'linear', 'nearest');
    get_T_sat = griddedInterpolant(P_vec, T_sat_grid, 'linear', 'nearest');
    get_rho_l = griddedInterpolant(P_vec, rho_l_grid, 'linear', 'nearest');
    get_rho_v = griddedInterpolant(P_vec, rho_v_grid, 'linear', 'nearest');
    get_surften = griddedInterpolant(P_vec, surften_grid, 'linear', 'nearest');
    get_h_fg = griddedInterpolant(P_vec, h_fg_grid, 'linear', 'nearest');
    % P_sat grid for nucleate boiling
    psat_index = ~isnan(T_sat_grid);
    T_sat_clean = T_sat_grid(psat_index);
    P_clean = P_vec(psat_index);
    get_P_sat = griddedInterpolant(T_sat_clean, P_clean, 'linear', 'nearest');

    save(table_filename,'get_rho', 'get_cp', 'get_mu', 'get_k', 'get_T_sat', ...
        'get_rho_l', 'get_rho_v', 'get_surften', 'get_h_fg', 'get_P_sat', 'P_vec', 'T_vec')
end

%% Material Properties (316 Stainless)
T_mp = 1643.15; % K
% Wall thermal conductivity k_w = f(T)
k_w_ref_temps = [-0.15, 19.85, 26.85, 76.85, 126.85, 226.85, 326.85, 426.85, 526.85, 626.85, ...
    726.85, 826.85, 926.85, 1026.9, 1126.9, 1226.9, 1326.9, 1370.9, 1398.9, 1426.9] + 273.15; % K
k_w_ref = [12.97, 13.31, 13.44, 14.32, 15.16, 16.8, 18.36, 19.87, 21.39, 22.79, 24.06, 25.46, 26.74, ...
    28.02, 29.32, 30.61, 31.86, 32.41, 26.9, 27.24];
r = 1e-4; % m, surface roughness
%get_k_w = interp1(k_w_ref_temps, k_w_ref, T_w_local, 'linear'); % callout

%% Output Arrays
T_hw_array = zeros(size(pos_i));
T_cw_array = zeros(size(pos_i));
T_bulk_array = zeros(size(pos_i));
P_array = zeros(size(pos_i));
P_loss_array = zeros(size(pos_i));
q_flux_array = zeros(size(pos_i));
h_g_array = zeros(size(pos_i));
h_c_array = zeros(size(pos_i));
CHF_array = zeros(size(pos_i));

%% Main Loop
P_guess = convpres(500, 'psi', 'Pa');
P_prev_guess = convpres(450, 'psi', 'Pa');
P_error = realmax;
P_prev_error = 0;
tol_P = 100; % Pa
iter_P = 0;
while abs(P_error) > tol_P % Pressure guess loop
    iter_P = iter_P + 1;
    P_loc = P_guess;
    T_bulk = T_amb; 

    for d = length(pos_i):-1:1 % Axial marching loop
        % Local Geometry, add channel height array earlier
        A_g_loc = A_gas(d);
        A_w_loc = A_w(d);
        D_g_loc = D_gas(d);
        cw = w_channel(d);
        ch = h_channel;
        D_h_loc = (2*cw*ch)/(cw+ch);
        A_c_cs = cw*ch; % Cross sectional area of channel
        A_c_cs_tot = A_c_cs * num_channel;
        if (d ~= 1) % If not at final chamber station
            A_c_cs_next = w_channel(d-1)*h_channel;
            D_h_loc_next = (2*w_channel(d-1)*h_channel/(w_channel(d-1)+h_channel));
        else % At final station
            A_c_cs_next = A_c_cs;
            D_h_loc_next = D_h_loc;
        end

        % Coolant properties
        rho_c = get_rho(P_loc, T_bulk);
        cp_c = get_cp(P_loc, T_bulk);
        mu_c = get_mu(P_loc, T_bulk);
        k_c = get_k(P_loc, T_bulk);
        T_sat = get_T_sat(P_loc);
        rho_c_l = get_rho_l(P_loc);
        rho_c_v = get_rho_v(P_loc);
        surften = get_surften(P_loc);
        h_fg = get_h_fg(P_loc);

        mu_g_local_us = mu_g_local(d) * 0.0559974; % Pa*s - (lb/in)*s
        cp_g_local_us = cp_g_local(d) * 0.00023885; % J/(kg*K) - Btu/(lb*deg F)
       
        vel_c = mdot_f / (A_c_cs_tot * rho_c); % coolant cs area instead
        Re = (rho_c * D_h_loc * vel_c) / mu_c;

        if (Re>4000) % Turbulent
            f_guess = 64/Re; % Initial guess
            f_error = realmax;
            tol_f = 0.0001;
            iter_F = 0;
            while abs(f_error) > tol_f
                iter_F = iter_F + 1;
                f_calc = 1/(-2*log10(2.51/(Re*sqrt(f_guess))+r/(D_h_loc*3.72)))^2;
                f_error = f_calc - f_guess;
                f_guess = 0.5*f_guess + 0.5*f_calc; % avoid overshooting
                if iter_F > 100
                    break;
                end
            end
        elseif (Re<2300) % Laminar
            f_calc = 64/Re;
        else
            f_calc = 0.02783 + (7.17*10^-6)*(Re - 2300);
        end

        Pr = (cp_c*mu_c)/k_c;
        Nu = ((f_calc/8)*(Re-1000)*Pr)/...
            (1+12.7*(f_calc/8)^0.5*(Pr^(2/3)-1)); % Gnielinski 
        h_c = Nu * k_c / D_h_loc; % Coolant convection htc hydraulic diameter
        h_c_array(d) = h_c;
        
        % Gas properties
        Taw_loc = Taw(d); % Local adiabatic wall temp

        % HW Temp Iteration
        T_hw_guess = 1315; % Initial Guess (1.25 FOS applied to material melting point)
        T_hw_prev_guess = 650; % Cooler lower bound (Luca's value)
        tol_q = 0.1; % Watts
        q_error = realmax; % placeholder
        q_prev_error = 0;
        iter_T = 0;
        while abs(q_error) > tol_q
            iter_T = iter_T + 1;
            % Fin Efficiency
            k_w_loc = interp1(k_w_ref_temps, k_w_ref, T_hw_guess, 'linear', 'extrap');
            fin_m = sqrt((2*h_c*(dx + w_rib))/(k_w_loc * dx * w_rib));
            fin_eff = tanh(fin_m * h_channel) / (fin_m * h_channel);
            % Gas convection HTC with Bartz
            sigma = 1 / ...
                ((0.5 * T_hw_guess/T_stag * (1 + (gamma-1)/2 * M_local(d)^2) + 0.5)^0.68 *...
                (1 + (gamma-1)/2 * M_local(d)^2)^0.12);
            h_g_us = ((0.026/D_t_us^0.2)*...
                (((mu_g_local_us^0.2)*cp_g_local_us)/prandtl_g_local(d)^0.6)*...
                (Pc_us/cstar_act_us)^0.8)*...
                (D_t_us/R_curve_us)^0.1*...
                (At/A_local(d))^0.9*...
                sigma;
            h_g = h_g_us * 2943611.72;
            
            q_1 = h_g*A_g_loc*(Taw_loc - T_hw_guess);

            T_cw = T_hw_guess - (q_1)*...
                (num_channel*log(1+2*wall_thickness/D_g_loc))/(2*pi*dx*k_w_loc); % Cold wall temp derived from guess, change dx to real dl

            % For wall-coolant HT: t_i = wall_thickness, Nch = num_channel, Di = D_g_loc, Li = dx
            if (T_cw <= T_sat)
                R_th = (1/(h_c*dx*(2*fin_eff*ch+cw)))+...
                    ((num_channel*log(1+2*wall_thickness/D_g_loc))/(2*pi*dx*k_w_loc));
                q_3 = num_channel * (T_hw_guess-T_bulk) / R_th;
            else % Nucleate
                P_sat_T_cw = get_P_sat(T_cw);
                h_nb = 0.00122*...
                    (((k_c^0.79)*(cp_c^0.45)*(rho_c_l^0.49))/...
                    ((surften^0.5)*(mu_c^0.29)*(h_fg^0.24)*(rho_c_v^0.24)))*...
                    ((T_cw-T_sat)^0.24)*...
                    (max(0,P_sat_T_cw - P_loc))^0.75;
                S = 1/(1+(2.53*(10^-6))*(Re^1.17));
                q_3 = A_co(d)*(h_c*(T_cw-T_bulk)+S*h_nb*(T_cw-T_sat));
            end
            
            current_q_error = q_1 - q_3;
            % Secant
            if iter_T == 1
                q_prev_error = current_q_error;
                temp_T = T_hw_guess;
                T_hw_guess = T_hw_prev_guess;
                T_hw_prev_guess = temp_T;
            else
                T_next = T_hw_guess - current_q_error * (T_hw_guess - T_hw_prev_guess) / (current_q_error - q_prev_error + 1e-10); % small buffer to avoid divide by 0
                T_next = max(T_bulk + 1, min(T_next, Taw_loc - 1));
                T_hw_prev_guess = T_hw_guess;
                q_prev_error = current_q_error;
                T_hw_guess = T_next;
            end
            q_error = current_q_error;

            % Break
            if iter_T > 100
                disp('T_hw failed to converge at station d =' + d);
                break;
            end
        end

        T_hw_array(d) = T_hw_guess;
        T_cw_array(d) = T_cw;
        h_g_array(d) = h_g;

        % Check CHF
        q_flux = q_1/(pi*(D_g_loc+wall_thickness*2)*dx);
        q_flux_array(d) = q_flux;
        F_p = 1.17-8.56*(10^(-4))*convpres(P_loc, 'Pa', 'psi');
        CHF_base = 0.1003+0.05264*sqrt(max(0,convvel(vel_c, 'm/s', 'ft/s')*...
            (convtemp(T_sat, 'K', 'F')-convtemp(T_bulk, 'K', 'F'))));
        CHF_array(d) = CHF_base*F_p*1635000;
        %{
        if (q_flux >= CHF) %prob change this to smth else
            error('CHF exceeded at station ' + d);
        end
        %}

        % Prepare for next station
        T_bulk_array(d) = T_bulk; % Store the updated bulk temperature
        T_bulk = T_bulk + q_1/(mdot_f*cp_c); % K
        
        % Calculate pressure losses 
        P_loss_viscous = (f_calc*rho_c*vel_c^2*dx)/(2*D_h_loc);

        if (d ~= length(pos_i))
            if (A_c_cs < A_c_cs_next)
                K = ((A_c_cs/A_c_cs_next)^2-1)^2;
            elseif (A_c_cs > A_c_cs_next)
                K = 0.5-0.167*(A_cs_next/A_c_cs)-...
                    0.125*(A_cs_next/A_c_cs)^2-...
                    0.208*(A_cs_next/A_c_cs)^3;
            else
                K = 0;
            end
            P_loss_area = 0.5*K*rho_c*vel_c^2;
        else
            P_loss_area = 0;
        end

        P_loss_mom = mdot_f^2*... % Assume den diff is negligible, unless can find a way to get next station den 
            (2/(A_c_cs*num_channel+A_c_cs_next*num_channel))*...
            (1/(rho_c*A_c_cs*num_channel) - 1/(rho_c*A_c_cs_next*num_channel));

        P_loss_tot = P_loss_mom + P_loss_area + P_loss_viscous;
        P_loss_array(d) = P_loss_tot;
        P_loc = P_loc - P_loss_tot; 
        P_array(d) = convpres(P_loc, 'Pa', 'psi');
    end

    current_P_error = P_loc - P_target;
    % Secant
    if iter_P == 1
        P_prev_error = current_P_error;
        temp_P = P_guess;
        P_guess = P_prev_guess;
        P_prev_guess = temp_P;
    else
        P_next = P_guess - current_P_error * (P_guess - P_prev_guess) / (current_P_error - P_prev_error + 1e-10);
        % Prevent negative or physically impossible pressure guesses
        P_next = max(Pamb * 6894.75, P_next);
        P_prev_guess = P_guess;
        P_prev_error = current_P_error;
        P_guess = P_next;
    end
    P_error = current_P_error;
    if iter_P > 50
        disp('Pressure loop failed to converge');
        break;
    end
end

%% Plots

% Temps
figure('Name', 'Temperature', 'Color', 'w');
hold on; grid on;
plot(pos_i, T_hw_array, 'r', 'LineWidth', 2);
plot(pos_i, T_cw_array, 'b', 'LineWidth', 2);
plot(pos_i, T_bulk_array, 'c', 'LineWidth', 2);
legend('Hot Wall', 'Cold Wall', 'Bulk Coolant')
title('Temperatures')
xlabel('Axial Position x (m)');
ylabel('Temperature (K)');
xline(0, 'k--', 'Throat', 'LabelVerticalAlignment', 'bottom', 'HandleVisibility', 'off');
exportgraphics(gcf, 'temperatures.pdf', 'ContentType','vector');
hold off;

% Pressure
figure('Name', 'Pressure', 'Color', 'w');
hold on; grid on;
plot(pos_i, P_array, 'b', 'LineWidth', 2);
title('Coolant Static Pressure')
xlabel('Axial Position x (m)');
ylabel('Pressure (psi)');
xline(0, 'k--', 'Throat', 'LabelVerticalAlignment', 'bottom', 'HandleVisibility', 'off');
exportgraphics(gcf, 'pressure.pdf', 'ContentType','vector');

hold off;

% Heat Flux
figure('Name', 'HeatFlux', 'Color', 'w');
hold on; grid on;
plot(pos_i, q_flux_array, 'r', 'LineWidth', 2);
plot(pos_i, CHF_array, 'k--', 'LineWidth', 1.5, 'DisplayName', 'CHF Limit');
title('Coolant Heat Flux')
xlabel('Axial Position x (m)');
ylabel('Heat Flux (W/m^2)');
xline(0, 'k--', 'Throat', 'LabelVerticalAlignment', 'bottom', 'HandleVisibility', 'off');
exportgraphics(gcf, 'heatflux.pdf', 'ContentType','vector');

hold off;

% HTC
figure('Name', 'GasHTC', 'Color', 'w');
hold on; grid on;
plot(pos_i, h_g_array, 'r', 'LineWidth', 2);
title('Gas Heat Transfer Coefficient')
xlabel('Axial Position x (m)');
ylabel('Heat Transfer Coefficient (W/m^2*K)');
xline(0, 'k--', 'Throat', 'LabelVerticalAlignment', 'bottom', 'HandleVisibility', 'off');
exportgraphics(gcf, 'gashtc.pdf', 'ContentType','vector');
hold off;

figure('Name', 'CoolantHTC', 'Color', 'w');
hold on; grid on;
plot(pos_i, h_c_array, 'b', 'LineWidth', 2);
title('Coolant Heat Transfer Coefficient')
xlabel('Axial Position x (m)');
ylabel('Heat Transfer Coefficient (W/m^2*K)');
xline(0, 'k--', 'Throat', 'LabelVerticalAlignment', 'bottom', 'HandleVisibility', 'off');
exportgraphics(gcf, 'coolhtc.pdf', 'ContentType','vector');
hold off;