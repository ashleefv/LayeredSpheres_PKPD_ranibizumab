clc; clear;

dose_common = [2 2 2 2 2 2 2];

cases = [
    % struct( ...
    %     "caseName","Without_DDS_dose", ...
    %     "DDS_geometry","NaN", ...
    %     "dose_in",(0.05: 0.01: 2), ...
    %     "radius_scale",0, ...
    %     "thickness_scale", 0)
    % 
    % struct( ...
    %     "caseName","With_DDS_dose", ...
    %     "DDS_geometry","Chitosan_PCL", ...
    %     "dose_in",(0.05: 0.01: 2), ...
    %     "radius_scale", 1, ...
    %     "thickness_scale", 1)

    struct( ...
        "caseName","Chitosan_single", ...
        "DDS_geometry","Chitosan", ...
        "dose_in",dose_common, ...
        "radius_scale",[0.0001 0.1 0.5 1 3 5 10], ...
        "thickness_scale",ones(1,7))

    struct( ...
       "caseName","PCL_single", ...
       "DDS_geometry","PCL", ...
       "dose_in",dose_common, ...
       "radius_scale",[0.0001 1 30 50 70 90 100], ...
       "thickness_scale",ones(1,7))

    struct( ...
        "caseName", "BiLayer_Changing_Chitosan", ...
        "DDS_geometry", "Chitosan_PCL", ...      % if your bilayer uses this flag
        "dose_in", dose_common, ...
        "radius_scale", [0.5 1 1.5 2 3 4 5], ...
        "thickness_scale", ones(1,7) )

    struct( ...
        "caseName", "BiLayer_Changing_PCL", ...
        "DDS_geometry", "Chitosan_PCL", ...      % or whatever your bilayer flag is
        "dose_in", dose_common, ...
        "radius_scale", ones(1,7), ...
        "thickness_scale", [0.01 0.1 1 10 20 25 30] )

    struct( ...
        "caseName", "BiLayer_changing_Both", ...
        "DDS_geometry", "Chitosan_PCL", ...      % or your bilayer geometry flag
        "dose_in", dose_common, ...
        "radius_scale", [0.1 0.5 1.5 2 3 4 5], ...
        "thickness_scale", [10 2 3 5 6 2 0.5] )

   struct( ...
        "caseName","Without_DDS_changing_kD", ...
        "DDS_geometry","NaN", ...
        "dose_in",[2 2 2 2], ...
        "radius_scale",0, ...
        "thickness_scale", 0)

    struct( ...
        "caseName","With_DDS_changing_kD", ...
        "DDS_geometry","Chitosan_PCL", ...
        "dose_in",[2 2 2 2], ...
        "radius_scale", ones(1,7), ...
        "thickness_scale", ones(1,7))

    struct( ...
        "caseName","Without_DDS_changing_koff", ...
        "DDS_geometry","NaN", ...
        "dose_in",[2 2 2 2], ...
        "radius_scale",0, ...
        "thickness_scale", 0)

    struct( ...
        "caseName","With_DDS_changing_Koff", ...
        "DDS_geometry","Chitosan_PCL", ...
        "dose_in",[2 2 2 2], ...
        "radius_scale", ones(1,7), ...
        "thickness_scale", ones(1,7))
];

% Time
tinitial = 0; tfinal = 2000; time_interval = 10;
t = linspace(tinitial, tfinal, tfinal*time_interval);
nt = numel(t);


for c = 1:numel(cases)

    DDS_geometry    = cases(c).DDS_geometry;
    dose_in         = cases(c).dose_in;
    radius_scale    = cases(c).radius_scale;
    thickness_scale = cases(c).thickness_scale;

 %===========================
 %For diffent K_D, K_off
 %===========================

if ismember(cases(c).caseName, ["With_DDS_changing_koff","Without_DDS_changing_koff"])
    k_off = [0.3456 0.864 1.728 3.456];
else
    k_off = 0.864;
end

if ismember(cases(c).caseName, ["With_DDS_changing_kD","Without_DDS_changing_kD"])   
    k_D = [19000 9500 4750 1900];
else
    k_D = 19000;
end

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

Data_time_at_target_ret_10 = [];
Data_time_at_target_ret_50 = [];
Data_time_at_target_vit_10 = []; 
Data_time_at_target_vit_50 = [];
Data_time_at_target_aq_10 = []; 
Data_time_at_target_aq_50 = [];


fprintf("\n=== Running case: %s | geometry: %s ===\n", cases(c).caseName, DDS_geometry);

nd = numel(dose_in);
radius_scale    = repmat(radius_scale(:).', 1, ceil(nd/numel(radius_scale)));
thickness_scale = repmat(thickness_scale(:).', 1, ceil(nd/numel(thickness_scale)));

radius_scale    = radius_scale(1:nd);
thickness_scale = thickness_scale(1:nd);

% geometry outputs
R1_all   = NaN(1,nd);
R2_all   = NaN(1,nd);
delR_all = NaN(1,nd);

% suppression outputs
Data_time_at_target_ret_10 = NaN(1,nd);
Data_time_at_target_ret_50 = NaN(1,nd);
Data_time_at_target_vit_10 = NaN(1,nd);
Data_time_at_target_vit_50 = NaN(1,nd);
Data_time_at_target_aq_10  = NaN(1,nd);
Data_time_at_target_aq_50  = NaN(1,nd);

% time series (optional but you save them, so they must exist)
C_vret_Data = NaN(nt, nd);
C_vvit_Data = NaN(nt, nd);
C_vaq_Data  = NaN(nt, nd);
C_rret_Data = NaN(nt, nd);
C_rvit_Data = NaN(nt, nd);
C_raq_Data  = NaN(nt, nd);

drug_dose_all   = cell(1, nd);   % each cell holds a vector
RealTime_all    = cell(1, nd);   % each cell holds a vector (days)


for i = 1:nd
    dose_specific = dose_in(i);
    par.k_off = k_off(min(i,end));
    par.k_D   = k_D(min(i,end));

    isNoDDS = strcmpi(DDS_geometry,"none") || strcmpi(DDS_geometry,"nan") || isempty(DDS_geometry);

    if isNoDDS
        % -------------------------
        % WITHOUT DDS
        % -------------------------
        Dose = (dose_specific*1e-3)/((48.35*1000)*(4.5E-3)) * 1e12;
        % no release profile
        rpar = [];   
        Ci = [v_ret_Initial,0,0,0,  v_vit_Initial,Dose,0,0,  v_aq_Initial,0,0,0];
        soln = ode45(@(tt,yy) ODEs(tt,yy,rpar, par), [tinitial tfinal], Ci);

        drug_dose_all{i}    = [];
        RealTime_all{i}     = [];

    else
        % -------------------------
        % WITH DDS
        % -------------------------
        [time_sec, drug_dose, initial_drug_dose, R1, R2, delR] = ...
        solve_FD_spheres_variable_diffusivity(dose_specific, tfinal, DDS_geometry, radius_scale(i), thickness_scale(i));

        RealTime_days = time_sec/(60*60*24);

        Dose_profile = (drug_dose*1e-3)/((48.35*1000)*(4.5E-3)) * 1e12;           % pM (or pM/day depending on your meaning)
        initDose_pM  = (initial_drug_dose*1e-3)/((48.35*1000)*(4.5E-3)) * 1e12;

        rpar = [Dose_profile(:), RealTime_days(:)];  % [amount, time]
        Ci = [v_ret_Initial,0,0,0,  v_vit_Initial,initDose_pM,0,0,  v_aq_Initial,0,0,0];
        soln = ode45(@(tt,yy) ODEs(tt,yy,rpar,par), [tinitial tfinal], Ci);

        R1_all(i)   = R1;
        R2_all(i)   = R2;
        delR_all(i) = delR;
        drug_dose_all{i}    = drug_dose;       
        RealTime_all{i}     = RealTime_days; 
end

    % ---- evaluate solution on your common t-grid ----
    C_vret = deval(soln,t,1);
    C_rret = deval(soln,t,2);
    C_vvit = deval(soln,t,5);
    C_rvit = deval(soln,t,6);
    C_vaq  = deval(soln,t,9);
    C_raq  = deval(soln,t,10);

    C_vret_Data(:,i) = C_vret;
    C_vvit_Data(:,i) = C_vvit;
    C_vaq_Data(:,i)  = C_vaq;
    C_rret_Data(:,i) = C_rret;
    C_rvit_Data(:,i) = C_rvit;
    C_raq_Data(:,i)  = C_raq;

    % ---- suppression time calculation (same for both) ----
    [lowest_vret, Iret] = min(C_vret);
    [lowest_vvit, Ivit] = min(C_vvit);
    [lowest_vaq,  Iaq ] = min(C_vaq);

    Index_min = 10;
    if lowest_vret <= 0.5*v_ret_Initial || lowest_vvit <= 0.5*v_vit_Initial || lowest_vaq <= 0.5*v_aq_Initial
        Index_min = min([Iret, Ivit, Iaq]);
end

    editedt   = t(Index_min:end);
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

    if ~isempty(ir10), Data_time_at_target_ret_10(i) = editedt(ir10); end
    if ~isempty(iv10), Data_time_at_target_vit_10(i) = editedt(iv10); end
    if ~isempty(ia10), Data_time_at_target_aq_10(i)  = editedt(ia10); end

    if ~isempty(ir50), Data_time_at_target_ret_50(i) = editedt(ir50); end
    if ~isempty(iv50), Data_time_at_target_vit_50(i) = editedt(iv50); end
    if ~isempty(ia50), Data_time_at_target_aq_50(i)  = editedt(ia50); end


    tret = NaN; tvit = NaN; taq = NaN;
    if ~isempty(ir10), tret = editedt(ir10); end
    if ~isempty(iv10), tvit = editedt(iv10); end
    if ~isempty(ia10), taq  = editedt(ia10); end

    fprintf("k_D=%g, k_off=%g | 10%% times (days): Ret=%.2f Vit=%.2f Aq=%.2f\n", ...
        par.k_D, par.k_off, tret, tvit, taq);

end

    fname = cases(c).caseName + ".mat";
    save(fname, ...
    'radius_scale','dose_in', ...
    'drug_dose_all','RealTime_all',...
    'Data_time_at_target_ret_10','Data_time_at_target_vit_10','Data_time_at_target_aq_10', ...
    'Data_time_at_target_ret_50','Data_time_at_target_vit_50','Data_time_at_target_aq_50', ...
    'C_vret_Data','C_vvit_Data','C_vaq_Data', ...
    'C_rret_Data','C_rvit_Data','C_raq_Data', ...
    'v_ret_Initial','v_vit_Initial','v_aq_Initial', ...
    't','thickness_scale','k_off','k_D','DDS_geometry', ...
    'R1_all','R2_all','delR_all');


     end
    
   
  

% ==========================================================================
                %ODEs function%
% ==========================================================================
function derivVector = ODEs(t,y,rpar,par)
         DrugAmountAtTime = 0;  % default: no DDS input

if nargin >= 3 && ~isempty(rpar)
    timeArray        = rpar(:,2);
    drugReleaseArray = rpar(:,1);
    DrugAmountAtTime = interp1(timeArray, drugReleaseArray, t, 'linear', 0);
end

%----------------------------------------------------------------
%Parameters

k_off = par.k_off;
k_D   = par.k_D;
k_on  = k_off./k_D;

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
%k_D = 19000; %pM
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
drvit_dt = (k_off*c_vit - 2*k_on*v_vit*r_vit) + (2*k_off*h_vit - k_on*r_vit*c_vit) + (S_ret/V_vit)*pr_ILM*r_ret - (S_ret/V_vit)*pr_ILM*r_vit - kr_el*r_vit + DrugAmountAtTime;
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



