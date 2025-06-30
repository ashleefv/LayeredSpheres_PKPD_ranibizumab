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

%Alter K_D values here
%k_D Values
k_D1 = 19000;    %pM
k_D2 = 9500;     %pM
k_D3 = 4750;     %pM
k_D4 = 1900;     %pM

k_DStore = [k_D1, k_D2, k_D3, k_D4];


% Time range
tinitial = 0; %days
tfinal=1000; %time (days)
time_interval = 10;
numberOfTimes = tfinal*time_interval;
t=linspace(tinitial, tfinal, numberOfTimes); %t is the time

%Parameters
k_off = 0.864; % day^-1

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


for i = 1:length(k_DStore) %for loop for each dose

k_D_value = k_DStore(i);
dose_specific = dose_in(i);
Dose = (dose_specific*10^-3)/((48.35*1000)*(4.5E-3)) * 10^12;

Ci=[v_ret_Initial,0,0,0,v_vit_Initial,Dose,0,0,v_aq_Initial,0,0,0]; %[VEGF, ranibizumab, VEGF-ranibizumab, ranibizumab-VEGF-ranibizumab] repeated 3 times for each chamber (Retina, Vitreous, Aqueous)
odeFunc = @(t,Ci) ODEs(t,Ci, k_D_value);
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

%disp(Index_min)
%disp(['The lowest value is:', num2str(lowest_vret), ' at index:', num2str(Index_vret)]);
%disp(['The lowest value is:', num2str(lowest_vvit), ' at index:', num2str(Index_vvit)]);
%disp(['The lowest value is:', num2str(lowest_vaq), ' at index:', num2str(Index_vaq)]);

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
%fprintf('With DDS 10 percent VEGF suppression for the retina chamber is: %.2f\n With DDS 10 percent VEGF suppression for the vitreous chamber is: %.2f\n With DDS 10 percent VEGF suppression for the aqueous chamber is: %.2f\n', time_at_target_ret_10, time_at_target_vit_10, time_at_target_aq_10);

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
%fprintf('With DDS 50 percent VEGF suppression for the retina chamber is: %.2f\n With DDS 50 percent VEGF suppression for the vitreous chamber is: %.2f\n With DDS 50 percent VEGF suppression for the aqueous chamber is: %.2f\n', time_at_target_ret_50, time_at_target_vit_50, time_at_target_aq_50);

Data_time_at_target_ret_10(i) = time_at_target_ret_10;
Data_time_at_target_ret_50(i) = time_at_target_ret_50;
Data_time_at_target_vit_10(i) = time_at_target_vit_10; 
Data_time_at_target_vit_50(i) = time_at_target_vit_50;
Data_time_at_target_aq_10(i) = time_at_target_aq_10; 
Data_time_at_target_aq_50(i) = time_at_target_aq_50;

end


% bar plot 10%
barWidth = 1;
figure(figure_count);
hold on
hBar1 = bar([Data_time_at_target_ret_10', Data_time_at_target_vit_10',Data_time_at_target_aq_10'], 'grouped', 'LineWidth', 1.5); 
xticks(1:length(dose_in)); 
xticklabels({'19000', '9500', '4750', '1900'});
ylabel('Pharmacodynamic Suppression Time (Days)');
xlabel('k_{D}');
set(hBar1, 'BarWidth', barWidth);
set(gca, 'FontSize', 12)
legend({'Retina 10%', 'Vitreous 10%', 'Aqueous 10%'}, 'FontSize',14, 'Location', 'northwest');
ylim([0,700])
pbaspect([1 1 1])
axis square
box on
exportgraphics(figure(figure_count),sprintf('bar_dose_response_non_DDS_10.png'), 'Resolution', 300)
figure_count = figure_count+1;
hold off

% bar plot 50%
barWidth = 1;
figure(figure_count);
hold on
hBar1 = bar([Data_time_at_target_ret_50', Data_time_at_target_vit_50',Data_time_at_target_aq_50'], 'grouped', 'LineWidth', 1.5); 
xticks(1:length(dose_in)); 
xticklabels({'19000', '9500', '4750', '1900'});
ylabel('Pharmacodynamic Suppression Time (Days)');
xlabel('k_{D}');
set(hBar1, 'BarWidth', barWidth);
set(gca, 'FontSize', 12)
legend({'Retina 50%', 'Vitreous 50%', 'Aqueous 50%'}, 'FontSize',14, 'Location', 'northwest');
ylim([0,700])
pbaspect([1 1 1])
axis square
box on
exportgraphics(figure(figure_count),sprintf('bar_dose_response_non_DDS_50.png'), 'Resolution', 300)
figure_count = figure_count+1;
hold off


% dose response VEGF
figure(figure_count);
hold on
for j = 1:length(dose_in)
plot(t,C_vret_Data(:,j), 'LineWidth', 2)
end
xlabel('Time(days)')
ylabel('Free VEGF (pM)')
title('Retina')
legend('k_{D} 19000', 'k_{D} 9500', 'k_{D} 4750', 'k_{D} 1900', 'FontSize',12, 'Location','southeast');
ylim([-0.6,60])
xlim([-3 300])
set(gca, 'FontSize', 20)
axis square
box on
exportgraphics(figure(figure_count),sprintf('Retina_dose_response_VEGF_nonDDS.png'), 'Resolution', 300)
figure_count = figure_count+1;
hold off

figure(figure_count);
hold on
for j = 1:length(dose_in)
plot(t,C_vvit_Data(:,j), 'LineWidth', 2)
end
xlabel('Time(days)')
ylabel('Free VEGF (pM)')
title('Vitreous')
ylim([-0.2,20])
xlim([-3 300])
set(gca, 'FontSize', 20)
axis square
box on
exportgraphics(figure(figure_count),sprintf('Vitreous_dose_response_VEGF_nonDDS.png'), 'Resolution', 300)
figure_count = figure_count+1;
hold off

figure(figure_count);
hold on
for j = 1:length(dose_in)
plot(t,C_vaq_Data(:,j), 'LineWidth', 2)
end
xlabel('Time(days)')
ylabel('Free VEGF (pM)')
title('Aqueous')
%legend('Dose 0.05 mg', 'Dose 0.1 mg', 'Dose 0.2 mg', 'Dose 0.5 mg', 'Dose 1 mg', 'Dose 2 mg', 'FontSize',12, 'Location','southeast');
ylim([-0.018,1.8])
xlim([-3 300])
set(gca, 'FontSize', 20)
axis square
box on
exportgraphics(figure(figure_count),sprintf('Aqueous_dose_response_VEGF_nonDDS.png'), 'Resolution', 300)
figure_count = figure_count+1;
hold off



save('without_DDS_kD', ...
    't', ...
    'C_vret_Data', 'C_rret_Data', ...
    'C_vvit_Data', 'C_rvit_Data', ...
    'C_vaq_Data', 'C_raq_Data', ...
    'Data_time_at_target_ret_10', 'Data_time_at_target_ret_50', ...
    'Data_time_at_target_vit_10', 'Data_time_at_target_vit_50', ...
    'Data_time_at_target_aq_10', 'Data_time_at_target_aq_50', ...
    'k_DStore', 'dose_in', 'DDS_geometry', ...
    'v_ret_Initial', 'v_vit_Initial', 'v_aq_Initial');


function derivVector=ODEs(t,y, k_D_value) %must be (denominator, numerators)

%----------------------------------------------------------------
%Parameters
k_off = 0.864; % day^-1

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
k_D = k_D_value; %pM
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