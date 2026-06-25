% Test script for 1d regen engine

clc; clear; close all;

%% Basic Parameters (change to arrays after)

F = 1480 * 4.44822; % Target thrust in N
Pamb = 9.94; % psia At 10,500 ft altitude
Tamb = 298; % K (from feed calc sheet)
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
Pc = Pc_us * 6894.7573; % Pa
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
mdot_f = mdot_total ./ (1 + o_f); % kg/s
rho_f_inj = 789 * mfrac_eth + 1000 * mfrac_h2o; % kg/m^3 (20C)
CdA_f_inj = mdot_f / sqrt(2 * rho_f_inj * dP_inj); % m^2 (Heritage)
P_target = Pc + dP_inj; % Pa manifold pressure (target value)

%% Throat Geometry (Rao)
At = (cstar_act * mdot_total)/Pc; % m^2
Rt = sqrt(At/pi); % m

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
L_star = 35; % in, optimal for ethanol/lox
id_chamber = 4.75; % in, heritage
V_total = convlength(L_star, 'in', 'm') * At; % m^3
Rc = convlength(id_chamber / 2, 'in', 'm'); % m
x_conv_tangent = -(1.5*Rt) * sin(conv_angle);
y_conv_tangent = Rt + (1.5 * Rt) * (1 - cos(conv_angle));
x_conv_start = x_conv_tangent - (Rc - y_conv_tangent) / tan(conv_angle);
len_conv = -x_conv_start; % m
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
    if x < x_conv_start % chamber
        pos_j(k) = Rc;
    elseif x >= x_conv_start && x < x_conv_tangent % straight converging cone
        pos_j(k) = Rc - tan(conv_angle) * (x - x_conv_start);
    elseif x >= x_conv_tangent && x < 0 % converging throat arc
        pos_j(k) = (Rt + 1.5*Rt) - sqrt((1.5*Rt)^2 - x^2); % Rt + arc length - height of curve at point
    elseif x >= 0 && x < xn % diverging throat arc
        pos_j(k) = (Rt + 0.382*Rt) - sqrt((0.382*Rt)^2 - x^2);
    else % Rao parabola approx
        pos_j(k) = a_Rao * x^2 + b_Rao * x + c_Rao;
    end
end

%% Cooling Channel and Outer Jacket Geometry
% Fixed channel height and wall thickness
min_tol = 0.001; % m 3d printer tolerance
D_engine = 2 * pos_j; % Array of diameter at every node
D_t = 2 * Rt;
h_channel = min_tol; % channel height (radial)
wall_thickness = min_tol;
w_rib = min_tol; % fixed, rib width
D_channel_base = D_engine + 2 * wall_thickness; % Engine diameters wtih added wall thickness
D_t_base = D_t + 2 * wall_thickness; % Throat diameter with added wall thickness
% Number of channels determined at throat
circ_t_base = pi * (D_t_base); % Gas side diameter + wall thickness
num_channel = floor(circ_t_base / (w_rib + min_tol));
% Local circumferences across engine
circ_local_base = pi * D_channel_base;
w_channel = (circ_local_base - (num_channel * w_rib)) ./ num_channel; % variable, channel widths
D_h = (4 .* w_channel .* h_channel) ./ (2 * w_channel + 2 * h_channel); % hydraulic diameter of rectangular channels
Per_heated = w_channel + 2 * h_channel; % heated perimeter 

%Areas
%Stot = Per_heated .* num_channel .* dx; % coolant side surface area (heated)
%Sg = pi .* D_engine .* dx; % gas side surface area
Aw_cool = w_channel .* dx; % cool wall area per increment
Afin = h_channel .* dx; % fin area per increment

%% Visualization Plot
%{
figure('Name', '1D Engine Geometry', 'Color', 'w');
hold on; grid on;
plot(pos_i, pos_j, 'k', 'LineWidth', 2);
plot(pos_i, -pos_j, 'k', 'LineWidth', 2);
title('Regen 1D Profile')
xlabel('Axial Position x (m)');
ylabel('Radial Position y (m)');
axis equal;
xline(0, 'r--', 'Throat', 'LabelVerticalAlignment', 'bottom');
hold off;
%} 

%% Gas Properties
% 1D interpolation from 3 CEA points for Cp, gamma, k, mu for Bartz
AR_ref = [Rc^2/Rt^2, 1, eps];
A_local = D_engine.^2 .* (pi/4);
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
Taw_local = T_stag * ((1 + prandtl_g_local.^(1/3).*((gamma - 1) / 2) .* M_local.^2) ...
            ./ (1 + ((gamma - 1) / 2) .* M_local.^2));

%% Coolant Properties

% Required properties update for equilibrium: Cp, k, mu, rho
% Define Table bounds (Tmax < Tsat at Pmin or coolprop will crash)
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
    Tsat_grid = zeros(1,res); % 1D
    % Generate values
    for i = 1:res
        P_val = P_vec(i);
        if P_val < 6400000 % Ethanol critical pressure in Pa
            Tsat_eth = py.CoolProp.CoolProp.PropsSI('T','P',P_val,'Q',0, 'ethanol');
            Tsat_h2o = py.CoolProp.CoolProp.PropsSI('T','P',P_val,'Q',0, 'water');
            Tsat_grid(i) = mfrac_eth * Tsat_eth + mfrac_h2o * Tsat_h2o;
        else % undefined for superheated vapor
            Tsat_grid(i) = NaN;
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
    % 2D interpolation objects
    get_rho = griddedInterpolant(P_grid, T_grid, rho_grid, 'linear', 'nearest');
    get_cp = griddedInterpolant(P_grid, T_grid, cp_grid, 'linear', 'nearest');
    get_mu = griddedInterpolant(P_grid, T_grid, mu_grid, 'linear', 'nearest');
    get_k = griddedInterpolant(P_grid, T_grid, k_grid, 'linear', 'nearest');
    get_Tsat = griddedInterpolant(P_vec, Tsat_grid, 'linear', 'nearest');
    save(table_filename,'get_rho', 'get_cp', 'get_mu', 'get_k', 'get_Tsat', 'P_vec', 'T_vec')
end

%% Material Properties (Copper)
T_mp = 1357.77; % K
% Wall thermal conductivity k_w = f(T)
mat_tstep = 20;
k_w_ref_temps = [273, 400, 600, 800, 1000, 1200];
k_w_ref = [401, 392, 383, 371, 357, 342];
get_k_w = griddedInterpolant(k_w_ref_temps, k_w_ref, 'linear'); % R^2 = 0.994

%% Main Loop
P_guess = convpres(500, 'psi', 'Pa');
for d = length(pos_i):1

end




