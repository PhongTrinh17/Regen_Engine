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
exp_ratio = double(c.get_eps_at_PcOvPe(pyargs('Pc', Pc_us, 'MR', o_f, 'PcOvPe',(Pc_us / Pamb))));
% transport properties [Cp, mu, k, Prandtl]
transport_chamber = c.get_Chamber_Transport(pyargs('Pc', Pc_us, 'MR', o_f, 'eps', exp_ratio));
transport_throat = c.get_Throat_Transport(pyargs('Pc', Pc_us, 'MR', o_f, 'eps', exp_ratio));
transport_exit = c.get_Exit_Transport(pyargs('Pc', Pc_us, 'MR', o_f, 'eps', exp_ratio));
cp_g_ref = [double(transport_chamber{1}), double(transport_throat{1}), double(transport_exit{1})] .* 4184; % J/kg*K
mu_g_ref = [double(transport_chamber{2}), double(transport_throat{2}), double(transport_exit{2})] .* 0.0001; % Pa*s viscosity
k_g_ref = [double(transport_chamber{3}), double(transport_throat{3}), double(transport_exit{3})] .* 0.4184; % W/m*K thermal conductivity
prandtl_ref = [double(transport_chamber{4}), double(transport_throat{4}), double(transport_exit{4})];

gamma_throat = c.get_Throat_MolWt_gamma(pyargs('Pc', Pc_us, 'MR', o_f, 'eps', exp_ratio));
gamma = double(gamma_throat{2}); % specific heat ratio, constant for isentropic relations

temps_cea = c.get_Temperatures(pyargs('Pc', Pc_us, 'MR', o_f, 'eps', exp_ratio)); % R
% derated by cstar_eff^2 (c* ~ sqrt(Tc)); adiabatic wall -> constant stagnation temp
T_stag = double(temps_cea{1}) * 5/9 * cstar_eff^2; % K

cstar_theo = double(c.get_Cstar(pyargs('Pc', Pc_us, 'MR', o_f))) * 0.3048; % m/s
cf_cea = c.get_PambCf(pyargs('Pamb', Pamb, 'Pc', Pc_us, 'MR', o_f, 'eps', exp_ratio));
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
Ae = At * exp_ratio; % m^2
R_exit = sqrt(Ae/pi); % m, nozzle exit radius
percent_len = 0.8; % input (0.8 is optimal fractional length for most cases)
len_div = percent_len * ((R_exit - Rt)/tan(deg2rad(15))); % m diverging length from radii
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
% Need to derive parabola that hits (xn, Rn), (xe, R_exit) w/ starting slope tan(theta_n)
matrix_A = [xn^2, xn, 1; x_exit^2, x_exit, 1; 2*xn, 1, 0];
matrix_B = [Rn; R_exit; tan(theta_n)];
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

dpos_dx = gradient(pos_j, dx);
dl = dx .* sqrt(1 + dpos_dx.^2);

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
num_channel = 50;
% Local circumferences across engine
circ_local_base = pi * D_channel_base;
w_channel = (circ_local_base - (num_channel * w_rib)) ./ num_channel; % variable, channel widths
assert(all(w_channel > 0), ...
    'Channel width non-positive (min = %.3g mm): num_channel/w_rib too large for throat circumference.', ...
    min(w_channel)*1000);
D_h = (4 .* w_channel .* h_channel) ./ (2 * w_channel + 2 * h_channel); % hydraulic diameter of rectangular channels
Per_heated = w_channel + 2 * h_channel; % heated perimeter

%% HT Areas
A_gas = pi .* D_gas .* dl; % Gas-wall convection SA
A_w = pi .* ((D_gas + D_channel_base)./2) .* dl; % wall-wall conduction SA (average diameter)
A_co = Per_heated .* num_channel .* dl; % coolant side surface area (heated)
A_wc = w_channel .* dl; % cool wall area per increment

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
table_filename = 'coolprop_tables.mat';
table_version_req = 4; % bump to force cached tables to rebuild
tables_loaded = false;
if isfile(table_filename)
    tbl = load(table_filename);
    if isfield(tbl, 'table_version') && tbl.table_version == table_version_req
        get_rho = tbl.get_rho; get_cp = tbl.get_cp; get_mu = tbl.get_mu; get_k = tbl.get_k;
        get_T_sat = tbl.get_T_sat; get_rho_l = tbl.get_rho_l; get_rho_v = tbl.get_rho_v;
        get_surften = tbl.get_surften; get_h_fg = tbl.get_h_fg; get_P_sat = tbl.get_P_sat;
        P_vec = tbl.P_vec; T_vec = tbl.T_vec;
        tables_loaded = true;
    end
end
if ~tables_loaded % create tables first time (or rebuild when table_version_req is bumped)
    res = 50; % 50x50 data grid
    P_min = convpres(300, 'psi', 'Pa'); % Pa, at chamber
    P_max = convpres(900, 'psi', 'Pa'); % Pa, at fuel tank
    T_min = 290; % K, standard inlet temp
    T_max = 500; % K, above mixture T_sat; bulk lookups below clamp to the liquid side
    % 1D vectors for P and T (fast solve for Coolprop, easy to get)
    P_vec = linspace(P_min, P_max, res);
    T_vec = linspace(T_min, T_max, res);
    x_eth = (mfrac_eth/46.068) / (mfrac_eth/46.068 + mfrac_h2o/18.015);
    x_h2o = 1 - x_eth;
    A12_vl = 1.6798; A21_vl = 0.9227;
    gam_eth = exp(A12_vl*(A21_vl*x_h2o/(A12_vl*x_eth + A21_vl*x_h2o))^2);
    gam_h2o = exp(A21_vl*(A12_vl*x_eth/(A12_vl*x_eth + A21_vl*x_h2o))^2);
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
        T_sat_h2o = py.CoolProp.CoolProp.PropsSI('T','P',P_val,'Q',0, 'water');
        if P_val < 6140000 % Ethanol critical pressure in Pa
            T_sat_eth = py.CoolProp.CoolProp.PropsSI('T','P',P_val,'Q',0, 'ethanol');
            % bubble point by bisection on modified Raoult; a mass-weighted
            % average of the pure Tsats is not a bubble point
            T_lo = T_sat_eth - 15;
            T_hi = min(T_sat_h2o, 513.5); % capped below ethanol critical temp
            for it_bp = 1:40
                T_mid = 0.5*(T_lo + T_hi);
                P_bub = x_eth*gam_eth*py.CoolProp.CoolProp.PropsSI('P','T',T_mid,'Q',0,'ethanol') + ...
                    x_h2o*gam_h2o*py.CoolProp.CoolProp.PropsSI('P','T',T_mid,'Q',0,'water');
                if P_bub < P_val
                    T_lo = T_mid;
                else
                    T_hi = T_mid;
                end
            end
            T_sat_grid(i) = 0.5*(T_lo + T_hi);
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
            % clamp lookups to the liquid side of each component's dome so the
            % grid never mixes liquid and vapor states
            if P_val < 6140000
                T_eth_q = min(T_val, T_sat_eth - 1);
            else % supercritical: single phase at any T
                T_eth_q = T_val;
            end
            T_h2o_q = min(T_val, T_sat_h2o - 1);

            rho_eth = py.CoolProp.CoolProp.PropsSI('D','T',T_eth_q,'P',P_val,'ethanol');
            rho_h2o = py.CoolProp.CoolProp.PropsSI('D','T',T_h2o_q,'P',P_val,'water');
            cp_eth = py.CoolProp.CoolProp.PropsSI('C','T',T_eth_q,'P',P_val,'ethanol');
            cp_h2o = py.CoolProp.CoolProp.PropsSI('C','T',T_h2o_q,'P',P_val,'water');
            mu_eth = py.CoolProp.CoolProp.PropsSI('V','T',T_eth_q,'P',P_val,'ethanol');
            mu_h2o = py.CoolProp.CoolProp.PropsSI('V','T',T_h2o_q,'P',P_val,'water');
            k_eth = py.CoolProp.CoolProp.PropsSI('L','T',T_eth_q,'P',P_val,'ethanol');
            k_h2o = py.CoolProp.CoolProp.PropsSI('L','T',T_h2o_q,'P',P_val,'water');
            
            rho_grid(i,j) =  1 / ((mfrac_eth / rho_eth) + (mfrac_h2o / rho_h2o));
            cp_grid(i,j) = mfrac_eth * cp_eth + mfrac_h2o * cp_h2o;
            % Grunberg-Nissan; mass-weighted mu is ~2x low for this blend
            G12 = 840 / T_val;
            mu_grid(i,j) = exp(x_eth*log(mu_eth) + x_h2o*log(mu_h2o) + x_eth*x_h2o*G12);
            k_grid(i,j) = mfrac_eth * k_eth + mfrac_h2o * k_h2o;
        end
    end
    % Interpolation objects
    get_rho = griddedInterpolant(P_grid, T_grid, rho_grid, 'linear', 'nearest');
    get_cp = griddedInterpolant(P_grid, T_grid, cp_grid, 'linear', 'nearest');
    get_mu = griddedInterpolant(P_grid, T_grid, mu_grid, 'linear', 'nearest');
    get_k = griddedInterpolant(P_grid, T_grid, k_grid, 'linear', 'nearest');
    % sub-critical points only; a NaN inside the grid bleeds into nearby queries
    sat_idx = ~isnan(T_sat_grid);
    get_T_sat = griddedInterpolant(P_vec(sat_idx), T_sat_grid(sat_idx), 'linear', 'nearest');
    get_rho_l = griddedInterpolant(P_vec(sat_idx), rho_l_grid(sat_idx), 'linear', 'nearest');
    get_rho_v = griddedInterpolant(P_vec(sat_idx), rho_v_grid(sat_idx), 'linear', 'nearest');
    get_surften = griddedInterpolant(P_vec(sat_idx), surften_grid(sat_idx), 'linear', 'nearest');
    get_h_fg = griddedInterpolant(P_vec(sat_idx), h_fg_grid(sat_idx), 'linear', 'nearest');
    % P_sat grid for nucleate boiling
    get_P_sat = griddedInterpolant(T_sat_grid(sat_idx), P_vec(sat_idx), 'linear', 'nearest');

    table_version = table_version_req;
    save(table_filename,'get_rho', 'get_cp', 'get_mu', 'get_k', 'get_T_sat', ...
        'get_rho_l', 'get_rho_v', 'get_surften', 'get_h_fg', 'get_P_sat', ...
        'P_vec', 'T_vec', 'table_version')
end
P_table_max = max(P_vec); % Pa, upper clamp for pressure guesses

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

% Test Arrays (remove when done)
f_array = zeros(size(pos_i));
Nu_array = zeros(size(pos_i));
T_sat_array = zeros(size(pos_i));
h_c_f_array = zeros(size(pos_i));
Re_array = zeros(size(pos_i));
vel_c_array = zeros(size(pos_i));

%% Main Loop
P_guess = convpres(500, 'psi', 'Pa');
P_prev_guess = convpres(450, 'psi', 'Pa');
P_prev_error = NaN;
tol_P = 100; % Pa
iter_P = 0;
max_iter_P = 50;
P_converged = false;
while ~P_converged % Pressure guess loop
    iter_P = iter_P + 1;
    P_loc = P_guess;
    T_bulk = T_amb; 

    for d = length(pos_i):-1:1 % Axial marching loop
        % Local Geometry, add channel height array earlier
        A_g_loc = A_gas(d);
        A_w_loc = A_w(d);
        D_g_loc = D_gas(d);
        dl_loc = dl(d);
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
        % a NaN here makes the solver loops below exit silently as if converged
        if ~all(isfinite([rho_c, cp_c, mu_c, k_c, T_sat, rho_c_l, rho_c_v, surften, h_fg]))
            error('Non-finite coolant property at station %d (P = %.4g Pa, T_bulk = %.1f K)', ...
                d, P_loc, T_bulk);
        end

        mu_g_local_us = mu_g_local(d) * 0.0559974; % Pa*s - (lb/in)*s
        cp_g_local_us = cp_g_local(d) * 0.00023885; % J/(kg*K) - Btu/(lb*deg F)
       
        vel_c = mdot_f / (A_c_cs_tot * rho_c); % coolant cs area instead
        Re = (rho_c * D_h_loc * vel_c) / mu_c;
        Re_array(d) = Re;
        vel_c_array(d) = vel_c;

        % rough f for pressure drop, smooth f for Gnielinski (a smooth-tube correlation)
        if (Re>4000) % Turbulent
            f_calc = colebrook_f(Re, r/D_h_loc);
            f_smooth = colebrook_f(Re, 0);
        elseif (Re<2300) % Laminar
            f_calc = 64/Re;
            f_smooth = f_calc;
        else % Transition
            wt_f = (Re - 2300)/(4000 - 2300);
            f_calc = (1-wt_f)*(64/2300) + wt_f*colebrook_f(4000, r/D_h_loc);
            f_smooth = (1-wt_f)*(64/2300) + wt_f*colebrook_f(4000, 0);
        end

        f_array(d) = f_calc;

        Pr = (cp_c*mu_c)/k_c;
        Nu_lam = 4.36; % laminar Nu, uniform heat flux
        if (Re > 4000)
            Nu = ((f_smooth/8)*(Re-1000)*Pr)/...
                (1+12.7*(f_smooth/8)^0.5*(Pr^(2/3)-1));
        elseif (Re < 2300) % laminar (Gnielinski would go negative for Re < 1000)
            Nu = Nu_lam;
        else
            Nu_turb = ((f_smooth/8)*(Re-1000)*Pr)/...
                (1+12.7*(f_smooth/8)^0.5*(Pr^(2/3)-1));
            Nu = Nu_lam + (Nu_turb - Nu_lam)*(Re - 2300)/(4000 - 2300);
        end
        h_c = Nu * k_c / D_h_loc; % Coolant convection htc hydraulic diameter
        h_c_array(d) = h_c;
        Nu_array(d) = Nu;
        
        % Gas properties
        Taw_loc = Taw(d); % Local adiabatic wall temp

        % HW Temp Iteration
        T_hw_guess = 1315; % Initial Guess (1.25 FOS applied to material melting point)
        T_hw_prev_guess = 650; % Cooler lower bound (Luca's value)
        tol_q = 0.1; % Watts
        q_prev_error = NaN;
        iter_T = 0;
        max_iter_T = 100;
        while true
            iter_T = iter_T + 1;
            % Fin Efficiency
            k_w_loc = interp1(k_w_ref_temps, k_w_ref, T_hw_guess, 'linear', 'extrap');
            % continuous rib, so no dx term (that would make fin_eff mesh-dependent)
            fin_m = sqrt(2*h_c/(k_w_loc*w_rib));
            fin_eff = tanh(fin_m * h_channel) / (fin_m * h_channel);
            h_c_f = h_c*(cw + 2*fin_eff*ch)/(cw + w_rib); % diagnostic only, not used in q_3
            h_c_f_array(d) = h_c_f;
            % Gas convection HTC with Bartz
            sigma = 1 / ...
                ((0.5 * (T_hw_guess/T_stag) * (1 + (gamma-1)/2 * M_local(d)^2) + 0.5)^0.68 *...
                (1 + (gamma-1)/2 * M_local(d)^2)^0.12);
            h_g = ((0.026/D_t^0.2)*...
                ((mu_g_local(d)^0.2*cp_g_local(d))/prandtl_g_local(d)^0.6)*...
                (Pc/cstar_act)^0.8)*...
                (D_t/R_curve)^0.1*...
                (At/A_local(d))^0.9*...
                sigma;
            
            q_1 = h_g*A_g_loc*(Taw_loc - T_hw_guess);

            T_cw = T_hw_guess - (q_1)*...
                log(1+2*wall_thickness/D_g_loc)/(2*pi*dl_loc*k_w_loc); % Cold wall temp derived from guess

            % both regimes share this convection term so q_3 stays continuous
            % across T_cw = T_sat; fin efficiency is applied once, here
            A_conv_eff = num_channel * (cw + 2*fin_eff*ch) * dl_loc;
            q_conv = h_c * A_conv_eff * (T_cw - T_bulk);
            if (T_cw <= T_sat)
                q_3 = q_conv;
            else % nucleate boiling
                P_sat_T_cw = get_P_sat(T_cw);
                dP_boil = max(P_sat_T_cw - P_loc, 0); % clamp: negative would make q_3 complex
                h_nb = 0.00122*...
                    (((k_c^0.79)*(cp_c^0.45)*(rho_c_l^0.49))/...
                    ((surften^0.5)*(mu_c^0.29)*(h_fg^0.24)*(rho_c_v^0.24)))*...
                    ((T_cw-T_sat)^0.24)*...
                    dP_boil^0.75;
                S = 1/(1+(2.53*(10^-6))*(Re^1.17));
                q_3 = q_conv + S*h_nb*A_conv_eff*(T_cw-T_sat);
            end

            current_q_error = q_1 - q_3;
            % NaN or complex would otherwise exit the loop silently and look converged
            if ~isfinite(current_q_error) || ~isreal(current_q_error)
                error('Non-finite/complex heat-flux residual at station %d (T_hw = %.1f K, P = %.4g Pa)', ...
                    d, T_hw_guess, P_loc);
            end
            if abs(current_q_error) <= tol_q
                break; % q_1, T_cw, h_g all correspond to this final T_hw_guess
            end
            if iter_T >= max_iter_T
                warning('T_hw failed to converge at station d = %d (residual = %.3g W)', ...
                    d, current_q_error);
                break;
            end
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
        end

        T_hw_array(d) = T_hw_guess;
        T_cw_array(d) = T_cw;
        h_g_array(d) = h_g;

        % Check CHF
        q_flux = q_1/A_conv_eff; % same reference area as q_3
        q_flux_array(d) = q_flux;
        dT_sub = T_sat - T_bulk; % K, bulk subcooling
        if dT_sub > 0
            CHF_base = 0.1003+0.05264*sqrt(convvel(vel_c, 'm/s', 'ft/s')*(dT_sub*9/5));
            F_p = 1.17-8.56*(10^(-4))*convpres(P_loc, 'Pa', 'psi');
            CHF_array(d) = CHF_base*F_p*1635000;
        else
            CHF_array(d) = NaN; % correlation invalid at zero subcooling
        end
        T_sat_array(d) = T_sat;
        %{
        if (q_flux >= CHF_array(d)) %prob change this to smth else
            error('CHF exceeded at station %d', d);
        end
        %}

        % Prepare for next station
        T_bulk_array(d) = T_bulk; % Store the updated bulk temperature
        T_bulk = T_bulk + q_1/(mdot_f*cp_c); % K
        
        % Calculate pressure losses 
        P_loss_viscous = (f_calc*rho_c*vel_c^2*dl_loc)/(2*D_h_loc);

        if (d ~= length(pos_i))
            if (A_c_cs < A_c_cs_next)
                K = ((A_c_cs/A_c_cs_next)^2-1)^2;
            elseif (A_c_cs > A_c_cs_next)
                K = 0.5-0.167*(A_c_cs_next/A_c_cs)-...
                    0.125*(A_c_cs_next/A_c_cs)^2-...
                    0.208*(A_c_cs_next/A_c_cs)^3;
            else
                K = 0;
            end
            P_loss_area = 0.5*K*rho_c*vel_c^2;
        else
            P_loss_area = 0;
        end

        % dP = G_avg*(V_next - V_here); relies on T_bulk already being advanced
        rho_c_next = get_rho(P_loc, T_bulk);
        P_loss_mom = mdot_f^2*...
            (2/(A_c_cs*num_channel+A_c_cs_next*num_channel))*...
            (1/(rho_c_next*A_c_cs_next*num_channel) - 1/(rho_c*A_c_cs*num_channel));

        P_loss_tot = P_loss_mom + P_loss_area + P_loss_viscous;
        P_loss_array(d) = P_loss_tot;
        P_loc = P_loc - P_loss_tot; 
        P_array(d) = convpres(P_loc, 'Pa', 'psi');
    end

    current_P_error = P_loc - P_target;
    % NaN or complex would otherwise exit the loop silently and look converged
    if ~isfinite(current_P_error) || ~isreal(current_P_error)
        error('Non-finite/complex pressure residual on iteration %d (P_guess = %.4g Pa)', ...
            iter_P, P_guess);
    end
    if abs(current_P_error) <= tol_P
        P_converged = true;
    elseif iter_P >= max_iter_P
        warning('Pressure loop failed to converge after %d iterations (error = %.3g Pa)', ...
            iter_P, current_P_error);
        break;
    elseif iter_P == 1
        P_prev_error = current_P_error;
        temp_P = P_guess;
        P_guess = P_prev_guess;
        P_prev_guess = temp_P;
    else
        P_next = P_guess - current_P_error * (P_guess - P_prev_guess) / (current_P_error - P_prev_error + 1e-10);
        P_next = max(P_target, min(P_next, P_table_max));
        P_prev_guess = P_guess;
        P_prev_error = current_P_error;
        P_guess = P_next;
    end
end

% the tables clamp outside their bounds with no warning
n_boil = nnz(T_bulk_array >= T_sat_array);
if n_boil > 0
    warning('T_bulk reached T_sat at %d stations: bulk boiling is not modeled', n_boil);
end
T_sat_ceiling = get_T_sat(P_table_max);
n_Tb_clamped = nnz(T_bulk_array > max(T_vec));
if n_Tb_clamped > 0
    warning('T_bulk exceeds the property-table ceiling (%.0f K) at %d stations (max %.1f K): bulk properties were clamped there', ...
        max(T_vec), n_Tb_clamped, max(T_bulk_array));
end
n_Tcw_clamped = nnz(T_cw_array > T_sat_ceiling);
if n_Tcw_clamped > 0
    warning('T_cw exceeds the saturation-table ceiling (%.1f K) at %d stations: h_nb used a clamped P_sat there', ...
        T_sat_ceiling, n_Tcw_clamped);
end
if convpres(min(P_array), 'psi', 'Pa') < min(P_vec) - 1
    warning('coolant pressure fell below the table floor (%.0f psi): properties were clamped there', ...
        convpres(min(P_vec), 'Pa', 'psi'));
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
plot(pos_i, h_c_f_array, 'g', 'LineWidth', 2);
title('Coolant Heat Transfer Coefficient')
xlabel('Axial Position x (m)');
ylabel('Heat Transfer Coefficient (W/m^2*K)');
xline(0, 'k--', 'Throat', 'LabelVerticalAlignment', 'bottom', 'HandleVisibility', 'off');
exportgraphics(gcf, 'coolhtc.pdf', 'ContentType','vector');
hold off;

% Friction Factor
figure('Name', 'Darcy Friction Factor', 'Color', 'w');
hold on; grid on;
plot(pos_i, f_array, 'b', 'LineWidth', 2);
title('Friction Factor')
xlabel('Axial Position x (m)');
ylabel('Dimensionless');
xline(0, 'k--', 'Throat', 'LabelVerticalAlignment', 'bottom', 'HandleVisibility', 'off');

% Nu
figure('Name', 'Nusselt Number', 'Color', 'w');
hold on; grid on;
plot(pos_i, Nu_array, 'b', 'LineWidth', 2);
title('Nusselt Number')
xlabel('Axial Position x (m)');
ylabel('Dimensionless');
xline(0, 'k--', 'Throat', 'LabelVerticalAlignment', 'bottom', 'HandleVisibility', 'off');

% Saturation Temp
figure('Name', 'Saturation Temp', 'Color', 'w');
hold on; grid on;
plot(pos_i, T_sat_array, 'b', 'LineWidth', 2);
title('Saturation Temp')
xlabel('Axial Position x (m)');
ylabel('K');
xline(0, 'k--', 'Throat', 'LabelVerticalAlignment', 'bottom', 'HandleVisibility', 'off');

% Re
figure('Name', 'Re', 'Color', 'w');
hold on; grid on;
plot(pos_i, Re_array, 'b', 'LineWidth', 2);
title('Re')
xlabel('Axial Position x (m)');
ylabel('Dimensionless');
xline(0, 'k--', 'Throat', 'LabelVerticalAlignment', 'bottom', 'HandleVisibility', 'off');

% Coolant Velocity
figure('Name', 'Coolant Velocity', 'Color', 'w');
hold on; grid on;
plot(pos_i, vel_c_array, 'b', 'LineWidth', 2);
title('Coolant Velocity')
xlabel('Axial Position x (m)');
ylabel('m/s');
xline(0, 'k--', 'Throat', 'LabelVerticalAlignment', 'bottom', 'HandleVisibility', 'off');

%% Local Functions
function f = colebrook_f(Re, rel_rough)
% rel_rough = r/D_h; pass 0 for smooth
f = 64/Re;
for it_f = 1:100
    f_new = 1/(-2*log10(rel_rough/3.7065 + 2.5226/(Re*sqrt(f))))^2;
    if abs(f_new - f) < 1e-4
        f = f_new;
        return;
    end
    f = 0.5*f + 0.5*f_new; % damped to avoid overshooting
end
end