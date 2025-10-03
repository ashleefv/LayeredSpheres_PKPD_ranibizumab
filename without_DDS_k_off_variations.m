clc;
clear;

% Removing existing PNG files from the current directory 
pngFiles = dir('*.png');
for i = 1:length(pngFiles)
    delete(pngFiles(i).name);
end

% DDS geometry
DDS_geometry = "Chitosan_PCL";
%DDS_geometry = "Chitosan";
%DDS_geometry = "PCL";


figure_count = 1;
Data_time_at_target_ret_10 = [];
Data_time_at_target_ret_50 = [];
Data_time_at_target_vit_10 = []; 
Data_time_at_target_vit_50 = [];
Data_time_at_target_aq_10 = []; 
Data_time_at_target_aq_50 = [];

dose_in = [2, 2, 2, 2];
radius_scale = [1, 1, 1, 1];
thickness_scale = [1, 1, 1, 1];

%Alter K_off values here
%k_off Values
k_off1 = 0.3456; 
k_off2 = 0.864;    %day^-1
%k_off2 = 1.296;    %day^-1
k_off3 = 1.728;    %day^-1
k_off4 = 3.456;    %day^-1

k_OffStore = [k_off1,k_off2,k_off3,k_off4];


% Time range
tinitial = 0; %days
tfinal=1000; %time (days)
time_interval = 10;
numberOfTimes = tfinal*time_interval;
t=linspace(tinitial, tfinal, numberOfTimes); %t is the time

%Parameters
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

%----------------------------------------------------------------

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


for i = 1:length(k_OffStore) %for loop for each dose

k_Off_value = k_OffStore(i);
dose_specific = dose_in(i);
Dose = (dose_specific*10^-3)/((48.35*1000)*(4.5E-3)) * 10^12;

Ci=[v_ret_Initial,0,0,0,v_vit_Initial,Dose,0,0,v_aq_Initial,0,0,0]; %[VEGF, ranibizumab, VEGF-ranibizumab, ranibizumab-VEGF-ranibizumab] repeated 3 times for each chamber (Retina, Vitreous, Aqueous)
odeFunc = @(t,Ci) ODEs(t,Ci, k_Off_value);
soln=ode45(odeFunc,t,Ci);

%Retina
C_vret=deval(soln,t,1);
C_rret=deval(soln,t,2);
C_cret=deval(soln,t,3);
C_hret=deval(soln,t,4);
%Vitreous
C_vvit=deval(soln,t,5);
C_rvit=deval(soln,t,6);
C_cvit=deval(soln,t,7);
C_hvit=deval(soln,t,8);
%Aqueous
C_vaq=deval(soln,t,9);
C_raq=deval(soln,t,10);
C_caq=deval(soln,t,11);
C_haq=deval(soln,t,12);

C_vret_Data(:,i) = C_vret;
C_vaq_Data(:,i) = C_vaq;
C_vvit_Data(:,i) = C_vvit;

C_rret_Data(:,i) = C_rret;
C_raq_Data(:,i) = C_raq;
C_rvit_Data(:,i) = C_rvit;


%Calculate index of lowest value for VEGF
[lowest_vret, Index_vret] = min(C_vret); 
[lowest_vvit, Index_vvit] = min(C_vvit);
[lowest_vaq, Index_vaq] = min(C_vaq);
Index_min = 10;

if lowest_vret <= 0.5 * v_ret_Initial || lowest_vvit <= 0.5 * v_vit_Initial || lowest_vaq <= 0.5 * v_aq_Initial
    Index_min = min([Index_vret, Index_vvit, Index_vaq]);
else
    beep
    'Drug loading is not enough to reduce VEGF to 50% of its original value'
end


%Calculates 10% Free VEGF Suppression Time
editedC_vret = C_vret(Index_min:end);
editedC_vvit = C_vvit(Index_min:end);
editedC_vaq = C_vaq(Index_min:end);
editedt = t(Index_min:end);

% Target concentrations for 10% suppression
target_concentration_ret_10 = 0.1 * v_ret_Initial;
target_concentration_vit_10 = 0.1 * v_vit_Initial;
target_concentration_aq_10 = 0.1 * v_aq_Initial;

% Find the times for 10% suppression
index_ret_10 = find(editedC_vret >= target_concentration_ret_10, 1);
time_at_target_ret_10 = editedt(index_ret_10);
index_vit_10 = find(editedC_vvit >= target_concentration_vit_10, 1);
time_at_target_vit_10 = editedt(index_vit_10);
index_aq_10 = find(editedC_vaq >= target_concentration_aq_10, 1);
time_at_target_aq_10 = editedt(index_aq_10);

% Target concentrations for 50% suppression
target_concentration_ret_50 = 0.5 * v_ret_Initial;
target_concentration_vit_50 = 0.5 * v_vit_Initial;
target_concentration_aq_50 = 0.5 * v_aq_Initial;

% Find the times for 50% suppression
index_ret_50 = find(editedC_vret >= target_concentration_ret_50, 1);
time_at_target_ret_50 = editedt(index_ret_50);
index_vit_50 = find(editedC_vvit >= target_concentration_vit_50, 1);
time_at_target_vit_50 = editedt(index_vit_50);
index_aq_50 = find(editedC_vaq >= target_concentration_aq_50, 1);
time_at_target_aq_50 = editedt(index_aq_50);

fprintf([' k_OffStore = %.4f  | 10%% VEGF suppression times (days): ', ...
         'Retina: %.4f, Vitreous: %.1f, Aqueous: %.1f\n'], ...
k_OffStore(i), time_at_target_ret_10, time_at_target_vit_10, time_at_target_aq_10);
fprintf([' k_OffStore = %.4f  | 50%% VEGF suppression times (days): ', ...
         'Retina: %.4f, Vitreous: %.1f, Aqueous: %.1f\n'], ...
        k_OffStore(i), time_at_target_ret_50, time_at_target_vit_50, time_at_target_aq_50);


Data_time_at_target_ret_10(i) = time_at_target_ret_10;
Data_time_at_target_ret_50(i) = time_at_target_ret_50;
Data_time_at_target_vit_10(i) = time_at_target_vit_10; 
Data_time_at_target_vit_50(i) = time_at_target_vit_50;
Data_time_at_target_aq_10(i) = time_at_target_aq_10; 
Data_time_at_target_aq_50(i) = time_at_target_aq_50;

end

save('without_DDS_koff.mat', ...
    't', ...
    'C_vret_Data', 'C_rret_Data', ...
    'C_vvit_Data', 'C_rvit_Data', ...
    'C_vaq_Data', 'C_raq_Data', ...
    'Data_time_at_target_ret_10', 'Data_time_at_target_ret_50', ...
    'Data_time_at_target_vit_10', 'Data_time_at_target_vit_50', ...
    'Data_time_at_target_aq_10', 'Data_time_at_target_aq_50', ...
    'k_OffStore', 'dose_in', 'DDS_geometry');


function derivVector=ODEs(t,y, k_Off_value) %must be (denominator, numerators)

%----------------------------------------------------------------
%Parameters
k_off = k_Off_value; % day^-1

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
k_D = 19000; %pM
V_in = 18.5; %pmol/day
tr_half = 7.5; %days


%----------------------------------------------------------------
% Initial Conditions
v_ret = y(1);   %VEGF
r_ret = y(2);   %ranibuzumab
c_ret = y(3);   %VEGF-ranibuzumab complex
h_ret = y(4);   %ranibizumab-VEGF-ranibuzumab complex
v_vit = y(5);   %VEGF
r_vit = y(6);   %ranibuzumab
c_vit = y(7);   %VEGF-ranibuzumab complex
h_vit = y(8);   %ranibizumab-VEGF-ranibuzumab complex
v_aq = y(9);    %VEGF
r_aq = y(10);   %ranibuzumab
c_aq = y(11);   %VEGF-ranibuzumab complex
h_aq = y(12);   %ranibizumab-VEGF-ranibuzumab complex

%----------------------------------------------------------------
% Algebraic relationships
%k_D = k_off/k_on
k_on = k_off/k_D;

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


%----------------------------------------------------------------
% ODEs
%Retina
dvret_dt = (k_off*c_ret - 2*k_on*v_ret*r_ret) - (S_ret/V_ret)*(pv_ILM + pv_RPE)*v_ret + (S_ret/V_ret)*pv_ILM*v_vit + (V_in/V_ret);
drret_dt = (k_off*c_ret - 2*k_on*v_ret*r_ret) + (2*k_off*h_ret - k_on*r_ret*c_ret) - (S_ret/V_ret)*(pr_ILM + pr_RPE)*r_ret + (S_ret/V_ret)*pr_ILM*r_vit;
dcret_dt = -(k_off*c_ret - 2*k_on*v_ret*r_ret) + (2*k_off*h_ret - k_on*r_ret*c_ret) - (S_ret/V_ret)*(pc_ILM + pc_RPE)*c_ret + (S_ret/V_ret)*pc_ILM*c_vit;
dhret_dt = -(2*k_off*h_ret - k_on*r_ret*c_ret) - (S_ret/V_ret)*(ph_ILM + ph_RPE)*h_ret + (S_ret/V_ret)*ph_ILM*h_vit;

%Vitreous
dvvit_dt = (k_off*c_vit - 2*k_on*v_vit*r_vit) + (S_ret/V_vit)*pv_ILM*v_ret - (S_ret/V_vit)*pv_ILM*v_vit - kv_el*v_vit;
drvit_dt = (k_off*c_vit - 2*k_on*v_vit*r_vit) + (2*k_off*h_vit - k_on*r_vit*c_vit) + (S_ret/V_vit)*pr_ILM*r_ret - (S_ret/V_vit)*pr_ILM*r_vit - kr_el*r_vit;
dcvit_dt = -(k_off*c_vit - 2*k_on*v_vit*r_vit) + (2*k_off*h_vit - k_on*r_vit*c_vit) + (S_ret/V_vit)*pc_ILM*c_ret - (S_ret/V_vit)*pc_ILM*c_vit - kc_el*c_vit;
dhvit_dt = -(2*k_off*h_vit - k_on*r_vit*c_vit) + (S_ret/V_vit)*ph_ILM*h_ret - (S_ret/V_vit)*ph_ILM*h_vit - kh_el*h_vit;

%Aqueous
dvaq_dt = (k_off*c_aq - 2*k_on*v_aq*r_aq) + (V_vit/V_aq)*kv_el*v_vit - (C_L/V_aq)*v_aq;
draq_dt = (k_off*c_aq - 2*k_on*v_aq*r_aq) + (2*k_off*h_aq - k_on*r_aq*c_aq) + (V_vit/V_aq)*kr_el*r_vit - (C_L/V_aq)*r_aq;
dcaq_dt = -(k_off*c_aq - 2*k_on*v_aq*r_aq) + (2*k_off*h_aq - k_on*r_aq*c_aq) + (V_vit/V_aq)*kc_el*c_vit - (C_L/V_aq)*c_aq;
dhaq_dt = -(2*k_off*h_aq - k_on*r_aq*c_aq) + (V_vit/V_aq)*kh_el*h_vit - (C_L/V_aq)*h_aq;


% Derivatives vector
derivVector = [dvret_dt, drret_dt, dcret_dt, dhret_dt, dvvit_dt, drvit_dt, dcvit_dt, dhvit_dt, dvaq_dt, draq_dt, dcaq_dt, dhaq_dt]';
end