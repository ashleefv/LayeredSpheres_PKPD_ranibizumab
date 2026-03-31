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


