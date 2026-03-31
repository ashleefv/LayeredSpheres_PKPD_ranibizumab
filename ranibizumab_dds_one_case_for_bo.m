

%%%%%%%%%% Ranibizumab_DDS_All_Case.m %%%%%%%%%%

%% Time
tinitial = 0; tfinal = 2500; time_interval = 10;
t = linspace(tinitial, tfinal, tfinal*time_interval);
nt = numel(t);

%% Fixed parameters
DDS_geometry = string('Chitosan_PCL');
k_off = 0.864;
k_D   = 19000;

pv_ILM = 1.89E-7*86400; %cm/day
pr_ILM = 1.89E-7*86400; %cm/day
ph_ILM = 1.73E-7*86400; %cm/day
pc_ILM = 1.79E-7*86400; %cm/day
pv_RPE = 2.66E-7*86400; %cm/day
pr_RPE = 2.63E-7*86400; %cm/day
pc_RPE = 2.28E-7*86400; %cm/day
ph_RPE = 2.06E-7*86400; %cm/day

Rv_h = 2.39; %nm
Rr_h = 2.45; %nm
Rc_h = 3.29; %nm
Rh_h = 4.07; %nm

V_ret = 0.22; %cc
V_vit = 4.5; %cc
V_aq = 0.16; %cc
S_ret = 9.71; %cm^2
C_L = 3.6; %cc/day
V_in_REFERENCE = 18.5; %fmol/day
V_in = 18.5; %pmol/day
tr_half = 7.5; %days

%-------------------------

%Calculating k_el
tv_half = (tr_half/Rr_h)*Rv_h;  %Days
tr_half = (tr_half/Rr_h)*Rr_h;  %Days
tc_half = (tr_half/Rr_h)*Rc_h;  %Days
th_half = (tr_half/Rr_h)*Rh_h;  %Days

lambda_v = (log(2)/tv_half);
lambda_r = (log(2)/tr_half);
lambda_c = (log(2)/tc_half);
lambda_h = (log(2)/th_half);

kv_el = lambda_v - (S_ret/V_vit)*pv_ILM + (S_ret*pv_ILM)^2/(V_vit*S_ret*(pv_ILM + pv_RPE) - V_vit*V_ret*lambda_v);
kr_el = lambda_r - (S_ret/V_vit)*pr_ILM + (S_ret*pr_ILM)^2/(V_vit*S_ret*(pr_ILM + pr_RPE) - V_vit*V_ret*lambda_r);
kc_el = lambda_c - (S_ret/V_vit)*pc_ILM + (S_ret*pc_ILM)^2/(V_vit*S_ret*(pc_ILM + pc_RPE) - V_vit*V_ret*lambda_c);
kh_el = lambda_h - (S_ret/V_vit)*ph_ILM + (S_ret*ph_ILM)^2/(V_vit*S_ret*(ph_ILM + ph_RPE) - V_vit*V_ret*lambda_h);

%Calculate Initial VEGF Concentration
E_Q = (1/(V_vit*kv_el))*(V_in/((S_ret*pv_RPE)/(V_vit*kv_el) + 1 + pv_RPE/pv_ILM));
v_ret_Initial = E_Q*(1 + (kv_el/((S_ret/V_vit)*pv_ILM)));
v_vit_Initial = E_Q;
v_aq_Initial = E_Q*(V_vit/C_L)*kv_el;

%% Simulation (WITH DDS branch)
par.k_off = k_off;
par.k_D   = k_D;

[time_sec, drug_dose, initial_drug_dose, R1, R2, delR] = ...
    solve_FD_spheres_variable_diffusivity(dose_in, tfinal, DDS_geometry, radius_scale, thickness_scale);

RealTime_days = time_sec/(60*60*24);

Dose_profile = (drug_dose*1e-3)/((48.35*1000)*(4.5E-3)) * 1e12;           % pM (or pM/day depending on your meaning)
initDose_pM  = (initial_drug_dose*1e-3)/((48.35*1000)*(4.5E-3)) * 1e12;

rpar = [Dose_profile(:), RealTime_days(:)];  % [amount, time]
Ci = [v_ret_Initial,0,0,0,  v_vit_Initial,initDose_pM,0,0,  v_aq_Initial,0,0,0];
soln = ode45(@(tt,yy) ODEs(tt,yy,rpar,par), [tinitial tfinal], Ci);

%% ---- evaluate solution on your common t-grid ----
C_vret = deval(soln,t,1);
C_rret = deval(soln,t,2);
C_vvit = deval(soln,t,5);
C_rvit = deval(soln,t,6);
C_vaq  = deval(soln,t,9);
C_raq  = deval(soln,t,10);

%% % ---- suppression time calculation ----
[lowest_vret, Iret] = min(C_vret);
[lowest_vvit, Ivit] = min(C_vvit);
[lowest_vaq,  Iaq ] = min(C_vaq);

Index_min = 10;
if lowest_vret <= 0.5*v_ret_Initial || lowest_vvit <= 0.5*v_vit_Initial || lowest_vaq <= 0.5*v_aq_Initial
    Index_min = min([Iret, Ivit, Iaq]);
end

editedt = t(Index_min:end);
eRet = C_vret(Index_min:end);
eVit = C_vvit(Index_min:end);
eAq  = C_vaq(Index_min:end);

% helper safe pick (avoid "Index exceeds" if not found)
f = @(vec,thr) find(vec >= thr, 1, 'first');

ir10 = f(eRet, 0.1*v_ret_Initial);
iv10 = f(eVit, 0.1*v_vit_Initial);
ia10 = f(eAq,  0.1*v_aq_Initial);

ir50 = f(eRet, 0.5*v_ret_Initial);
iv50 = f(eVit, 0.5*v_vit_Initial);
ia50 = f(eAq,  0.5*v_aq_Initial);

Data_time_at_target_ret_10 = NaN; Data_time_at_target_ret_50 = NaN;
Data_time_at_target_vit_10 = NaN; Data_time_at_target_vit_50 = NaN;
Data_time_at_target_aq_10  = NaN; Data_time_at_target_aq_50  = NaN;

if ~isempty(ir10), Data_time_at_target_ret_10 = editedt(ir10); end
if ~isempty(iv10), Data_time_at_target_vit_10 = editedt(iv10); end
if ~isempty(ia10), Data_time_at_target_aq_10  = editedt(ia10); end
if ~isempty(ir50), Data_time_at_target_ret_50 = editedt(ir50); end
if ~isempty(iv50), Data_time_at_target_vit_50 = editedt(iv50); end
if ~isempty(ia50), Data_time_at_target_aq_50  = editedt(ia50); end


%%%%%%%%%% Table_Figure_together_Ranibizumab.m %%%%%%%%%%

%% Geometry metrics 
geo = lower(strtrim(DDS_geometry));
total_thickness = max(R2 - R1, 0);
if contains(geo,"chitosan_pcl") || (contains(geo,"chitosan") && contains(geo,"pcl"))
    delR_col = total_thickness;
elseif contains(geo,"pcl")
    delR_col = R2;
    thickness_scale = 0;
    radius_scale = 0;
elseif contains(geo,"chitosan")
    delR_col = 0;
    thickness_scale = 0;
else
    delR_col = 0;
end
fractionPCL = delR_col / R2;

%% Building new row and appending to CSV

existingT = readtable(csv_name, "VariableNamingRule", "preserve");
iteration = max(existingT.iteration) + 1;

colNames = {'10%_retina','50%_retina','10%_vitreous','50%_vitreous', ...
            '10%_aqueous','50%_aqueous','dose','k_d_value','k_off_value', ...
            'r_c=5.1_multiplier','del_r=1.25_multiplier','dds_geometry', ...
            'total_r1','total_r2','total_delr','fraction_pcl','dataset', 'iteration'};

newRow = table( ...
    Data_time_at_target_ret_10, Data_time_at_target_ret_50, ...
    Data_time_at_target_vit_10, Data_time_at_target_vit_50, ...
    Data_time_at_target_aq_10,  Data_time_at_target_aq_50, ...
    double(dose_in), k_D, k_off, ...
    radius_scale, thickness_scale, DDS_geometry, ...
    R1, R2, delR_col, fractionPCL, string("BO_suggestion"), iteration, ...
    'VariableNames', colNames);

updatedT  = [existingT; newRow];
writetable(updatedT, csv_name);

%fprintf("New row appended!\n");
%fprintf("10%% suppression times (days): Ret=%.3f Vit=%.3f Aq=%.3f\n", ...
%    Data_time_at_target_ret_10, Data_time_at_target_vit_10, Data_time_at_target_aq_10);
%fprintf("50%% suppression times (days): Ret=%.3f Vit=%.3f Aq=%.3f\n", ...
%    Data_time_at_target_ret_50, Data_time_at_target_vit_50, Data_time_at_target_aq_50);