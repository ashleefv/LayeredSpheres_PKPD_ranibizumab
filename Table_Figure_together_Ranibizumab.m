
clear all

colNames = {'10% Retina','50% Retina','10% Vitreous','50% Vitreous', ...
            '10% Aqueous','50% Aqueous','Dose','K_D Value','K_off Value', ...
            'R_C=5.1 Multiplier','del_R=1.25 Multiplier', 'DDS_Geometry', ...
            'Total R1','Total R2','Total delR', 'Fraction PCL', 'Dataset'};

T = cell2table(cell(0, numel(colNames)), 'VariableNames', colNames);

files = { ...
    "Without_DDS_dose.mat", ...
    "With_DDS_dose.mat",...
    "BiLayer_Changing_Chitosan.mat", ...
    "BiLayer_Changing_PCL.mat", ...
    "BiLayer_Changing_Both.mat", ...
    "PCL_single.mat", ...
    "Chitosan_single.mat",...
    "Without_DDS_changing_kD", ...
    "With_DDS_changing_kD", ...
    "Without_DDS_changing_koff", ...
    "With_DDS_changing_koff"
    };

for i = 1:numel(files)
    S = load(files{i});

    
    dose_in = S.dose_in(:);      
    n = numel(dose_in);

    [~, baseName, ~] = fileparts(files{i});
    dataset = repmat(string(baseName), n, 1);


  %  ----- time outputs (default 0 if missing) -----
    ret10 = getvec(S,"Data_time_at_target_ret_10", n);
    ret50 = getvec(S,"Data_time_at_target_ret_50", n);
    vit10 = getvec(S,"Data_time_at_target_vit_10", n);
    vit50 = getvec(S,"Data_time_at_target_vit_50", n);
    aq10  = getvec(S,"Data_time_at_target_aq_10",  n);
    aq50  = getvec(S,"Data_time_at_target_aq_50",  n);

  %  ----- parameters (default 0 if missing) -----
    k_D            = getscalar(S,"k_D",0);
    k_off          = getscalar(S,"k_off",0);
    radius_scale   = getrow(S,"radius_scale",n);      % row 1xn
    thickness_scale= getrow(S,"thickness_scale",n);   % row 1xn

    %----- DDS info (defaults for without_DDS) -----
    DDS_geometry = "none";
    if isfield(S,"DDS_geometry"), DDS_geometry = string(S.DDS_geometry); end

    R1_all = getvec(S,"R1_all",n);
    R2_all = getvec(S,"R2_all",n);
  
    
    % --- Radius rules ---
geo = lower(strtrim(string(DDS_geometry)));
total_thickness = max(R2_all - R1_all, 0);

if contains(geo,"chitosan_pcl") || (contains(geo,"chitosan") && contains(geo,"pcl"))
    delR_all = total_thickness;   
elseif contains(geo,"pcl")
    delR_all = R2_all;   
    thickness_scale = zeros(n,1);
    radius_scale = zeros(n,1);
elseif contains(geo,"chitosan")
    delR_all = zeros(n,1);
    thickness_scale = zeros(n,1);
    
else
    delR_all = zeros(n,1);
end


fractionPCL = delR_all./R2_all;
    

    kD_row = getparam(S,"k_D",n);   % returns n×1
    koff_row   = getparam(S,"k_off", n);
    koff_row   = repmat(k_off, n, 1);
    rscale_row = radius_scale(:);        % n x 1
    tscale_row = thickness_scale(:);     % n x 1
    DDS_type   = repmat(DDS_geometry, n, 1);

    newRow = table(ret10, ret50, vit10, vit50, aq10, aq50, ...
        dose_in, kD_row, koff_row, rscale_row, tscale_row, DDS_type, ...
        R1_all, R2_all, delR_all, fractionPCL, dataset,...
        'VariableNames', colNames);

    T = [T; newRow];
end

T
% Export to CSV
 writetable(T, 'VEGF_suppression_all.csv');

 %-------- helper functions --------
function v = getvec(S, fname, n)
    if isfield(S,fname)
        v = S.(fname);
        v = v(:);
    else
        v = zeros(n,1);
    end
    if numel(v) ~= n
        v = zeros(n,1); % fallback if sizes don't match
    end
end

function x = getscalar(S, fname, defaultVal)
    if isfield(S,fname) && isscalar(S.(fname))
        x = S.(fname);
    else
        x = defaultVal;
    end
end

function r = getrow(S, fname, n)
    if isfield(S,fname)
        tmp = S.(fname);
        tmp = tmp(:).';
        if numel(tmp) == 1
            r = repmat(tmp, 1, n);
        elseif numel(tmp) == n
            r = tmp;
        else
            r = zeros(1,n);
        end
    else
        r = zeros(1,n);
    end
end

function v = getparam(S, fname, n)
    if isfield(S,fname)
        tmp = S.(fname);
        tmp = tmp(:);
        if isscalar(tmp)
            v = repmat(tmp, n, 1);
        elseif numel(tmp) == n
            v = tmp;
        else
            v = zeros(n,1);
        end
    else
        v = zeros(n,1);
    end
end


%%%%%%%%%%%%%%%% figure 3 %%%%%%%%%%%%%%%%%%%%%%%%%%%

load ("Without_DDS_dose.mat")

% --- Load patient data once ---
M = readmatrix('Patient40Data.xlsx');
t_pat    = M(:,1);   % time (days)
vegf_pat = M(:,2);   % free VEGF (pM). If your data is in col 5, use M(:,5)

figure(3);
figname = 'figure3';
figure_count = 1;
% Figure 3 for 0.5 mg
for j = 1:length(dose_in)
if (dose_in(j) == 0.5)

subplot(2,3,figure_count)
yyaxis left
semilogy(t,C_rret_Data(:,j), 'LineWidth', 2)
xlabel('Time(days)')
ylabel('Ranibizumab (pM)')
title('Retina')
ylim([1, 1E7])
xlim([-3 300])

yyaxis right
plot(t,C_vret_Data(:,j), 'LineWidth', 2)
ylabel('Free VEGF (pM)')
ylim([0,60])
set(gca, 'FontSize', 16)

% reference levels
y10 = 0.10 * v_ret_Initial;
y50 = 0.50 * v_ret_Initial;
% draw horizontal lines on the RIGHT axis
h10 = yline(y10, '--', '', ...
    'LabelHorizontalAlignment','right', 'LabelVerticalAlignment','bottom', ...
    'HandleVisibility','off', 'LineWidth',2.2);
h50 = yline(y50, ':',  '', ...
    'LabelHorizontalAlignment','right', 'LabelVerticalAlignment','bottom', ...
    'HandleVisibility','off', 'LineWidth',2.2);

axis square
box on
figure_count = figure_count+1;

subplot(2,3,figure_count)
yyaxis left
semilogy(t,C_rvit_Data(:,j), 'LineWidth', 2)
xlabel('Time(days)')
ylabel('Ranibizumab (pM)')
title('Vitreous')
ylim([1, 1E7]);
xlim([-3 300])

yyaxis right
plot(t,C_vvit_Data(:,j), 'LineWidth', 2)
ylabel('Free VEGF (pM)')
set(gca, 'FontSize', 16)

% reference levels
y10 = 0.10 * v_vit_Initial;
y50 = 0.50 * v_vit_Initial;
% draw horizontal lines on the RIGHT axis
h10 = yline(y10, '--', '10% of baseline', ...
    'LabelHorizontalAlignment','right', 'LabelVerticalAlignment','bottom', ...
    'HandleVisibility','off', 'LineWidth',2.2);
h50 = yline(y50, ':',  '50% of baseline', ...
    'LabelHorizontalAlignment','right', 'LabelVerticalAlignment','bottom', ...
    'HandleVisibility','off', 'LineWidth',2.2);

axis square
box on
figure_count = figure_count+1;

subplot(2,3,figure_count)
yyaxis left
semilogy(t,C_raq_Data(:,j), 'LineWidth', 2)
xlabel('Time(days)')
ylabel('Ranibizumab (pM)')
title('Aqueous')
ylim([1, 1E7]);
xlim([-3 300])

%--------------------------------------
yyaxis right
plot(t, C_vaq_Data(:,j), 'LineWidth', 2, 'DisplayName','Free VEGF')  % your curve
hold on
scatter(t_pat, vegf_pat, 36, 'o', ...
    'MarkerFaceColor','none', ...     % <-- no fill
    'MarkerEdgeColor','k', ...        % edge color (pick what you like)
    'LineWidth',1.2, ...
    'DisplayName','Patient data');        % overlay points
ylabel('Free VEGF (pM)')
set(gca, 'FontSize', 16)
% reference levels 
y10 = 0.10 * v_aq_Initial;
y50 = 0.50 * v_aq_Initial;

% draw horizontal lines on the RIGHT axis
h10 = yline(y10, '--', '', ...
    'LabelHorizontalAlignment','right', 'LabelVerticalAlignment','bottom', ...
    'HandleVisibility','off', 'LineWidth',2.2);
h50 = yline(y50, ':',  '', ...
    'LabelHorizontalAlignment','right', 'LabelVerticalAlignment','bottom', ...
    'HandleVisibility','off', 'LineWidth',2.2);
axis square
box on
legend({'Ranibizumab', 'Free VEGF', 'Patitent Data'}, 'FontSize',10, 'Location', 'southeast');
%exportgraphics(figure(figure_count),sprintf('Retina%d.png',dose_in(j)), 'Resolution', 300)
figure_count = figure_count+1;
hold off
end
end 


% Figure 3 for 0.5 mg
load("With_DDS_dose.mat")
for j = 1:length(dose_in)
if (dose_in(j) == 0.5)

subplot(2,3,figure_count)
yyaxis left
semilogy(t,C_rret_Data(:,j), 'LineWidth', 2)
xlabel('Time(days)')
ylabel('Ranibizumab (pM)')
title('Retina with DDS')
ylim([1, 1E7])
xlim([-3 300])

yyaxis right
plot(t,C_vret_Data(:,j), 'LineWidth', 2)
ylabel('Free VEGF (pM)')
ylim([0,60])
set(gca, 'FontSize', 16)

% reference levels
y10 = 0.10 * v_ret_Initial;
y50 = 0.50 * v_ret_Initial;
% draw horizontal lines on the RIGHT axis
h10 = yline(y10, '--', '', ...
    'LabelHorizontalAlignment','right', 'LabelVerticalAlignment','bottom', ...
    'HandleVisibility','off', 'LineWidth',2.2);
h50 = yline(y50, ':',  '', ...
    'LabelHorizontalAlignment','right', 'LabelVerticalAlignment','bottom', ...
    'HandleVisibility','off', 'LineWidth',2.2);

axis square
box on
figure_count = figure_count+1;

subplot(2,3,figure_count)
yyaxis left
semilogy(t,C_rvit_Data(:,j), 'LineWidth', 2)
xlabel('Time(days)')
ylabel('Ranibizumab (pM)')
title('Vitreous with DDS')
ylim([1, 1E7])
xlim([-3 300])

yyaxis right
plot(t,C_vvit_Data(:,j), 'LineWidth', 2)
ylabel('Free VEGF (pM)')
set(gca, 'FontSize', 16)

% reference levels
y10 = 0.10 * v_vit_Initial;
y50 = 0.50 * v_vit_Initial;
% draw horizontal lines on the RIGHT axis
h10 = yline(y10, '--', '', ...
    'LabelHorizontalAlignment','right', 'LabelVerticalAlignment','bottom', ...
    'HandleVisibility','off', 'LineWidth',2.2);
h50 = yline(y50, ':',  '', ...
    'LabelHorizontalAlignment','right', 'LabelVerticalAlignment','bottom', ...
    'HandleVisibility','off', 'LineWidth',2.2);
axis square
box on
figure_count = figure_count+1;

subplot(2,3,figure_count)
yyaxis left
semilogy(t,C_raq_Data(:,j), 'LineWidth', 2)
xlabel('Time(days)')
ylabel('Ranibizumab (pM)')
title('Aqueous with DDS')
ylim([1, 1E7])
xlim([-3 300])

yyaxis right
plot(t,C_vaq_Data(:,j), 'LineWidth', 2)
ylabel('Free VEGF (pM)')
set(gca, 'FontSize', 16)

% reference levels
y10 = 0.10 * v_aq_Initial;
y50 = 0.50 * v_aq_Initial;
% draw horizontal lines on the RIGHT axis
h10 = yline(y10, '--', '', ...
    'LabelHorizontalAlignment','right', 'LabelVerticalAlignment','bottom', ...
    'HandleVisibility','off', 'LineWidth',2.2);
h50 = yline(y50, ':',  '', ...
    'LabelHorizontalAlignment','right', 'LabelVerticalAlignment','bottom', ...
    'HandleVisibility','off', 'LineWidth',2.2);

axis square
box on
figure_count = figure_count+1;
end
end



labelstring = {'A', 'B', 'C', 'D','E','F'};
for v = 1:6
    subplot(2,3,v)
    hold on
    text(-0.2, 1.1, labelstring(v)', 'Units', 'normalized', 'FontWeight', 'bold','FontSize', 14)
end

% widthInches = 15;
% heightInches = 10;
% run('ScriptForExportingImages.m')



%%%%%%%%%%%%%%%% figure 4 %%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(4);
figname = 'figure4';
figure_count = 1;

load ("Without_DDS_dose.mat")


% bar plot 10%

subplot(2,2,figure_count);
hold on
scatter(dose_in, Data_time_at_target_ret_10, 60, 'red', 'LineWidth', 1.2);
scatter(dose_in, Data_time_at_target_vit_10, 60, 'cyan', 'LineWidth', 1.2);
scatter(dose_in, Data_time_at_target_aq_10,  60, 'black', 'LineWidth', 1.2);

ylabel('Pharmacodynamic Suppression Time (Days)');
xlabel('Drug Amount (mg)');
title('Without DDS');
set(gca, 'FontSize', 12)
legend({'Retina 10%','Vitreous 10%','Aqueous 10%'}, 'FontSize',14, 'Location','northwest');

xlim([min(dose_in) max(2.1)]);
ylim([0 450]);
pbaspect([1 1 1])
axis square
box on
figure_count = figure_count+1;
hold off

% bar plot 50%

subplot(2,2,figure_count);
hold on
dose_all = dose_in(:);
Data_time_at_target_ret_50 = Data_time_at_target_ret_50(:);
Data_time_at_target_vit_50 = Data_time_at_target_vit_50(:);
Data_time_at_target_aq_50  = Data_time_at_target_aq_50(:);

scatter(dose_in, Data_time_at_target_ret_50, 60, 'red', 'LineWidth', 1.2);
scatter(dose_in, Data_time_at_target_vit_50, 60, 'cyan','LineWidth', 1.2);
scatter(dose_in, Data_time_at_target_aq_50,  60, 'black', 'LineWidth', 1.2);

ylabel('Pharmacodynamic Suppression Time (Days)');
xlabel('Drug Amount (mg)');
title('Without DDS');
set(gca, 'FontSize', 12)
legend({'Retina 50%','Vitreous 50%','Aqueous 50%'}, 'FontSize',14, 'Location','northwest');

xlim([min(dose_in) 2.1]);   
ylim([0 450]);
axis square; box on
pbaspect([1 1 1])
axis square
box on
figure_count = figure_count+1;
hold off



load("With_DDS_dose.mat")


% bar plot 10%
subplot(2,2,figure_count);
hold on
scatter(dose_in, Data_time_at_target_ret_10, 60, 'red', 'LineWidth', 1.2);
scatter(dose_in, Data_time_at_target_vit_10, 60, 'cyan', 'LineWidth', 1.2);
scatter(dose_in, Data_time_at_target_aq_10,  60, 'black', 'LineWidth', 1.2);

ylabel('Pharmacodynamic Suppression Time (Days)');
xlabel('Drug Amount (mg)');
title('With DDS');
set(gca, 'FontSize', 12)
legend({'Retina 10%','Vitreous 10%','Aqueous 10%'}, 'FontSize',14, 'Location','northwest');

xlim([min(dose_in) max(2.1)]);
ylim([0 450]);
pbaspect([1 1 1])
axis square
box on
figure_count = figure_count+1;
hold off

% bar plot 50%
subplot(2,2,figure_count);
hold on
scatter(dose_in, Data_time_at_target_ret_50, 60, 'red', 'LineWidth', 1.2);
scatter(dose_in, Data_time_at_target_vit_50, 60, 'cyan','LineWidth', 1.2);
scatter(dose_in, Data_time_at_target_aq_50,  60, 'black', 'LineWidth', 1.2);

ylabel('Pharmacodynamic Suppression Time (Days)');
xlabel('Drug Amount (mg)');
title('With DDS');
set(gca, 'FontSize', 12)
legend({'Retina 50%','Vitreous 50%','Aqueous 50%'}, 'FontSize',14, 'Location','northwest');

xlim([min(dose_in) 2.1]);   
ylim([0 450]);
axis square; box on
pbaspect([1 1 1])
axis square
box on
figure_count = figure_count+1;
hold off



labelstring = {'A', 'B', 'C', 'D'};
for v = 1:4
    subplot(2,2,v)
    hold on
    text(-0.2, 1.1, labelstring(v)', 'Units', 'normalized', 'FontWeight', 'bold','FontSize', 14)
end

% widthInches = 10;
% heightInches = 10;
% run('ScriptForExportingImages.m')

%%%%%%%%%%%%%%%% figure 5 %%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(5);
figname = 'figure5';
figure_count = 1;

load ("Without_DDS_dose.mat")

x= [0.05,0.1,0.5,1,2];
[tfound, idx] = ismember(x, dose_in);
if ~all(tfound)
    error("Some requested doses are not in dose_in: %s", mat2str(x(~tfound)));
end

% dose response VEGF
subplot(2,3,figure_count);
hold on
for k = 1:numel(x)
    j = idx(k);
    plot(t, C_vret_Data(:,j), 'LineWidth', 2)
end

xlabel('Time(days)')
ylabel('Free VEGF (pM)')
title('Retina')
legend('0.05 mg', '0.1 mg', '0.5 mg', '1 mg', '2 mg', 'FontSize',12, 'Location','southeast');
ylim([-0.6,60])
xlim([-3 300])
set(gca, 'FontSize', 16)

% reference levels
y10 = 0.10 * v_ret_Initial;
y50 = 0.50 * v_ret_Initial;
% draw horizontal lines on the RIGHT axis
h10 = yline(y10, '--', '', ...
    'LabelHorizontalAlignment','right', 'LabelVerticalAlignment','bottom', ...
    'HandleVisibility','off', 'LineWidth',2.2);
h50 = yline(y50, ':',  '', ...
    'LabelHorizontalAlignment','right', 'LabelVerticalAlignment','bottom', ...
    'HandleVisibility','off', 'LineWidth',2.2);

axis square
box on
figure_count = figure_count+1;
hold off

subplot(2,3,figure_count);
hold on
for k = 1:numel(x)
    j = idx(k);
    plot(t, C_vvit_Data(:,j), 'LineWidth', 2)
end
xlabel('Time(days)')
ylabel('Free VEGF (pM)')
title('Vitreous')
ylim([-0.2,20])
xlim([-3 300])
set(gca, 'FontSize', 16)

% reference levels
y10 = 0.10 * v_vit_Initial;
y50 = 0.50 * v_vit_Initial;
% draw horizontal lines on the RIGHT axis
h10 = yline(y10, '--', '10% of baseline', ...
    'LabelHorizontalAlignment','right', 'LabelVerticalAlignment','bottom', ...
    'HandleVisibility','off', 'LineWidth',2.2);
h50 = yline(y50, ':',  '50% of baseline', ...
    'LabelHorizontalAlignment','right', 'LabelVerticalAlignment','bottom', ...
    'HandleVisibility','off', 'LineWidth',2.2);

axis square
box on
figure_count = figure_count+1;
hold off

subplot(2,3,figure_count);
hold on
for k = 1:numel(x)
    j = idx(k);
    plot(t, C_vaq_Data(:,j), 'LineWidth', 2)
end
xlabel('Time(days)')
ylabel('Free VEGF (pM)')
title('Aqueous')
ylim([-0.018,1.8])
xlim([-3 300])
set(gca, 'FontSize', 16)

% reference levels
y10 = 0.10 * v_aq_Initial;
y50 = 0.50 * v_aq_Initial;
% draw horizontal lines on the RIGHT axis
h10 = yline(y10, '--', '', ...
    'LabelHorizontalAlignment','right', 'LabelVerticalAlignment','bottom', ...
    'HandleVisibility','off', 'LineWidth',2.2);
h50 = yline(y50, ':',  '', ...
    'LabelHorizontalAlignment','right', 'LabelVerticalAlignment','bottom', ...
    'HandleVisibility','off', 'LineWidth',2.2);

axis square
box on
figure_count = figure_count+1;
hold off

% Figure 5 
load("With_DDS_dose.mat")

x= [0.05,0.1,0.5,1,2];
[tfound, idx] = ismember(x, dose_in);
if ~all(tfound)
    error("Some requested doses are not in dose_in: %s", mat2str(x(~tfound)));
end

% dose response VEGF
subplot(2,3,figure_count);
hold on
for k = 1:numel(x)
    j = idx(k);
    plot(t, C_vret_Data(:,j), 'LineWidth', 2)
end
xlabel('Time(days)')
ylabel('Free VEGF (pM)')
title('Retina with DDS')
ylim([-0.6,60])
xlim([-3 300])
set(gca, 'FontSize', 16)

% reference levels
y10 = 0.10 * v_ret_Initial;
y50 = 0.50 * v_ret_Initial;
% draw horizontal lines on the RIGHT axis
h10 = yline(y10, '--', '', ...
    'LabelHorizontalAlignment','right', 'LabelVerticalAlignment','bottom', ...
    'HandleVisibility','off', 'LineWidth',2.2);
h50 = yline(y50, ':',  '', ...
    'LabelHorizontalAlignment','right', 'LabelVerticalAlignment','bottom', ...
    'HandleVisibility','off', 'LineWidth',2.2);

axis square
box on
figure_count = figure_count+1;
hold off

subplot(2,3,figure_count);
hold on
for k = 1:numel(x)
    j = idx(k);
    plot(t, C_vvit_Data(:,j), 'LineWidth', 2)
end
xlabel('Time(days)')
ylabel('Free VEGF (pM)')
title('Vitreous with DDS')
ylim([-0.2,20])
xlim([-3 300])
set(gca, 'FontSize', 16)

% reference levels
y10 = 0.10 * v_vit_Initial;
y50 = 0.50 * v_vit_Initial;
% draw horizontal lines on the RIGHT axis
h10 = yline(y10, '--', '', ...
    'LabelHorizontalAlignment','right', 'LabelVerticalAlignment','bottom', ...
    'HandleVisibility','off', 'LineWidth',2.2);
h50 = yline(y50, ':',  '', ...
    'LabelHorizontalAlignment','right', 'LabelVerticalAlignment','bottom', ...
    'HandleVisibility','off', 'LineWidth',2.2);

axis square
box on
figure_count = figure_count+1;
hold off


subplot(2,3,figure_count);
hold on
for k = 1:numel(x)
    j = idx(k);
    plot(t, C_vaq_Data(:,j), 'LineWidth', 2)
end
xlabel('Time(days)')
ylabel('Free VEGF (pM)')
title('Aqueous with DDS')
ylim([-0.018,1.8])
xlim([-3 300])
set(gca, 'FontSize', 16)

% reference levels
y10 = 0.10 * v_aq_Initial;
y50 = 0.50 * v_aq_Initial;
% draw horizontal lines on the RIGHT axis
h10 = yline(y10, '--', '', ...
    'LabelHorizontalAlignment','right', 'LabelVerticalAlignment','bottom', ...
    'HandleVisibility','off', 'LineWidth',2.2);
h50 = yline(y50, ':',  '', ...
    'LabelHorizontalAlignment','right', 'LabelVerticalAlignment','bottom', ...
    'HandleVisibility','off', 'LineWidth',2.2);

axis square
box on
figure_count = figure_count+1;
hold off



labelstring = {'A', 'B', 'C', 'D','E','F'};
for v = 1:6
    subplot(2,3,v)
    hold on
    text(-0.2, 1.1, labelstring(v)', 'Units', 'normalized', 'FontWeight', 'bold','FontSize', 14)
end
% 
% widthInches = 15;
% heightInches = 10;
% run('ScriptForExportingImages.m')

% %%%%%%%%%%%%%%%% figure 6 %%%%%%%%%%%%%%%%%%%%%%%%%%%


figure(6);
figname = 'figure6';
figure_count = 1;

load("With_DDS_dose.mat")
x = [0.05, 0.1, 0.5, 1, 2];

% Find which runs correspond to x doses
[tfound, idx] = ismembertol(x, dose_in, 1e-12);
if ~all(tfound)
    error("Some requested doses are not in dose_in: %s", mat2str(x(~tfound)));
end

for k = 1:numel(idx)
    i = idx(k);   % <-- this is the run index in your saved arrays

    RealTime  = RealTime_all{i};
    drug_dose = drug_dose_all{i};

    subplot(1,5,figure_count)
    set(gca, 'FontSize', 10)
    hold on; box on

    plot(RealTime(2:end), drug_dose(2:end), 'LineWidth', 2)
    xlim([-2 200])
    xlabel('Time (days)')

    legend(sprintf('Dose %0.2f mg', dose_in(i)), 'Location', 'northeast')
    axis square
    pbaspect([1 1 1])
    figure_count = figure_count + 1;
    hold off
end

labelstring = {'A', 'B', 'C', 'D','E',};
for v = 1:5
    subplot(1,5,v)
    hold on
    text(-0.2, 1.1, labelstring(v)', 'Units', 'normalized', 'FontWeight', 'bold','FontSize', 10)
end

% widthInches = 25;
% heightInches = 5;
% run('ScriptForExportingImages.m')


% %%%%%%%%%%%%%%%% figure 7 & 9  %%%%%%%%%%%%%%%%%%%%%%%%%%%


figure(7);
figname = 'figure7';
figure_count = 1;

T=readtable('VEGF_suppression_all.csv');
FilteredT = T(strcmp(T.Dataset, 'Chitosan_single')|strcmp(T.Dataset, 'PCL_single')|strcmp(T.Dataset,'BiLayer_Changing_PCL')|strcmp(T.Dataset,'BiLayer_Changing_Chitosan') |strcmp(T.Dataset,'BiLayer_Changing_Both'),:)

T = FilteredT;

datasets = { ...
    'Chitosan_single', ...
    'PCL_single', ...
    'BiLayer_Changing_PCL', ...
    'BiLayer_Changing_Chitosan', ...
    'BiLayer_Changing_Both' ...
};

C = lines(numel(datasets));

plots = { ...
    "x10_Retina",   "10% VEGF Suppression in Retina"; ...
    "x10_Vitreous", "10% VEGF Suppression in Vitreous"; ...
    "x10_Aqueous",  "10% VEGF Suppression in Aqueous"; ...
    "x50_Retina",   "50% VEGF Suppression in Retina"; ...
    "x50_Vitreous", "50% VEGF Suppression in Vitreous"; ...
    "x50_Aqueous",  "50% VEGF Suppression in Aqueous" ...
};


legendLabels = { ...
    'Chitosan Base Single-layer DDS', ...
    'PCL Base Single-layer DDS', ...
    'Bilayer DDS-Changing PCL layer', ...
    'Bilayer DDS-Changing Chitosan layer', ...
    'Bilayer DDS Changing Both' ...
};

%figure;
tl = tiledlayout(2,3, "Padding","compact", "TileSpacing","compact");

labelFS_main  = 12;   % main subplot label font size
labelFS_inset = 9;    % inset label font size
markers = {'o','s','^','d','v'};

hLeg = gobjects(numel(datasets),1);   


for p = 1:size(plots,1)
    zName = plots{p,1};
    ttl   = plots{p,2};

    ax = nexttile(tl);
    axList (p) = ax;
    hold(ax,'on'); grid(ax,'on'); box(ax,'on');

    for k = 1:numel(datasets)
        idx = strcmp(T.Dataset, datasets(k));

        x = T{idx,"FractionPCL"};
        y = T{idx,"TotalR2"} * 1e4;   
        z = T{idx, zName};

        mask = isfinite(x) & isfinite(y) & isfinite(z);
        if ~any(mask), continue; end

        c = rescale(z(mask), 0, max(z(mask)));
        scatter(ax, x(mask), y(mask), 25, c, 'filled', 'Marker', markers{k});
    end

   
    colorbar(ax);
    colormap(ax, jet);

    xlabel(ax,"Fraction PCL", "FontWeight","bold", "FontSize",labelFS_main);
    ylabel(ax,"Total Radius (\mum)", "FontWeight","bold", "FontSize",labelFS_main);
    title(ax, ttl);
end
    

% ----  bottom legend ----
axLeg = axes('Position',[0 0 1 1], 'Visible','off', 'HitTest','off');
hold(axLeg,'on');
hDummy = gobjects(numel(datasets),1);

for k = 1:numel(datasets)
    hDummy(k) = plot(axLeg, nan, nan, markers{k}, 'MarkerSize',10, 'LineStyle','none','LineWidth',1.2);
end
lgd = legend(axLeg, hDummy, legendLabels, ...
    'Orientation','horizontal','NumColumns',numel(legendLabels), ...
    'Location','southoutside');
lgd.Box = 'off';

lgd.FontWeight = 'bold';
lgd.FontSize   = labelFS_main; 
% ----  Add insets ----
drawnow;  

tl.TileSpacing = "compact";
tl.Padding     = "compact";
tl.Units = 'normalized';
tlPos = tl.Position;
tlPos(2) = tlPos(2) + 0.04;    
tlPos(4) = tlPos(4) - 0.04;    
tl.Position = tlPos;
%-----------
%--------

insetW = 200; insetH = 130; padR = 8; padT = 8;

for p = 1:numel(axList)
    ax = axList(p);

label = char('A' + p - 1);
t = text(ax, -0.10, 1.10, label, 'Units','normalized', ...
    'FontWeight','bold','FontSize',14, ...
    'HorizontalAlignment','left','VerticalAlignment','top');
     t.Clipping = 'off';


    ax.Units = 'pixels';
    pos = ax.Position;
    ti  = ax.TightInset;

    inner = [pos(1)+ti(1), pos(2)+ti(2), ...
             pos(3)-ti(1)-ti(3), pos(4)-ti(2)-ti(4)];

    x0 = inner(1) + inner(3) - insetW - padR;
    y0 = inner(2) + inner(4) - insetH - padT;

    axInset = axes('Parent', gcf, 'Units','pixels', 'Position',[x0 y0 insetW insetH]);
    hold(axInset,'on'); box(axInset,'on'); grid(axInset,'on');

    
    hSc = findobj(ax, 'Type','Scatter');
    copyobj(hSc, axInset);


    xlim(axInset, xlim(ax));
    ylim(axInset, [0 0.003*1e4]);

    xlabel(axInset,"Fraction PCL", "FontWeight","bold", "FontSize",labelFS_inset);
    ylabel(axInset,"Total Radius (\mum)", "FontWeight","bold", "FontSize",labelFS_inset);

    colormap(axInset, colormap(ax));
    caxis(axInset, caxis(ax));
    axInset.FontSize = 7;
    uistack(axInset,'top');
end



% widthInches = 15;
% heightInches = 10;
% run('ScriptForExportingImages.m')



%%%%%%%%%%%%%%%% figure 8 %%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(8);
figname = 'figure8';

figure_count = 1;


load("Chitosan_single.mat")

for i = 1:length(dose_in)
    RealTime = RealTime_all{i};    
    drug_dose = drug_dose_all{i}; 

% plot_DDS_drug_release_dynamics
subplot(2,7,figure_count)


set(gca, 'FontSize', 12)
hold on
box on
plot(RealTime(2:end), drug_dose(2:end), 'LineWidth', 2)
xlim([-2 200])
xlabel('Time (days)')
%ylabel('DDS Drug Release Rate (mg/s)');
legend(sprintf('%0.4fR_C', radius_scale(i)), 'Location', 'northeast')
pbaspect([1 1 1])
axis square
figure_count = figure_count+1;
hold off
end



load ("PCL_single.mat")

% plot_DDS_drug_release_dynamics
for i = 1:length(dose_in)
    RealTime = RealTime_all{i};     % Extract ith drug release time profile
    drug_dose = drug_dose_all{i}; 
subplot(2,7,figure_count)
set(gca, 'FontSize', 12)
hold on
box on
plot(RealTime(2:end), drug_dose(2:end), 'LineWidth', 2)
xlim([-2 200])
xlabel('Time (days)')
%ylabel('DDS Drug Release Rate (mg/s)');
legend(sprintf('%0.4fR_P', radius_scale(i)), 'Location', 'northeast')
pbaspect([1 1 1])
axis square
figure_count = figure_count+1;
hold off
end

labelstring = {'A', 'B', 'C', 'D','E', 'F','G','H', 'I','J','K','L','M','N'};
for v = 1:14
    subplot(2,7,v)
    hold on
    text(-0.2, 1.1, labelstring(v)', 'Units', 'normalized', 'FontWeight', 'bold','FontSize', 10)
end

% widthInches = 35;
% heightInches = 10;
% run('ScriptForExportingImages.m')


% %%%%%%%%%%%%%%%% figure 9 %%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(9);
figname = 'figure9';

figure_count = 1;

load("BiLayer_Changing_Chitosan.mat")


for i = 1:length(dose_in)
    RealTime = RealTime_all{i};    
    drug_dose = drug_dose_all{i}; 

    % plot_DDS_drug_release_dynamics
subplot(3,7,figure_count);
set(gca, 'FontSize', 10)
hold on
box on
plot(RealTime(2:end), drug_dose(2:end), 'LineWidth', 2)
xlim([-2 200])
xlabel('Time (days)')
%ylabel('DDS Drug Release Rate (mg/s)');
legend(sprintf('%0.1fR_C', radius_scale(i)), 'Location', 'northeast')
pbaspect([1 1 1])
axis square
figure_count = figure_count+1;
hold off
end


load("BiLayer_Changing_PCL.mat")
for i = 1:length(dose_in)
    RealTime = RealTime_all{i};    
    drug_dose = drug_dose_all{i}; 

% plot_DDS_drug_release_dynamics
subplot(3,7,figure_count);
set(gca, 'FontSize', 10)
hold on
box on
plot(RealTime(2:end), drug_dose(2:end), 'LineWidth', 2)
xlim([-2 200])
xlabel('Time (days)')
%ylabel('DDS Drug Release Rate (mg/s)');
legend(sprintf('%0.2fΔR', thickness_scale(i)), 'Location', 'northeast')
pbaspect([1 1 1])
axis square
figure_count = figure_count+1;
hold off

end


load("BiLayer_changing_Both.mat")
for i = 1:length(dose_in)
    RealTime = RealTime_all{i};    
    drug_dose = drug_dose_all{i}; 

   % plot_DDS_drug_release_dynamics
subplot(3,7,figure_count);
set(gca, 'FontSize', 10)
hold on
box on
plot(RealTime(2:end), drug_dose(2:end), 'LineWidth', 2)
xlim([-2 200])
xlabel('Time (days)')
%ylabel('DDS Drug Release Rate (mg/s)');
legend(sprintf('%0.1fR_C, %0.1fΔR', radius_scale(i),thickness_scale(i)), 'Location', 'northeast')
pbaspect([1 1 1])
axis square
figure_count = figure_count+1;
hold off
end

labelstring = {'A', 'B', 'C', 'D','E', 'F','G','H', 'I','J','K','L','M','N','O','P','Q','R','S','T','U'};
for v = 1:21
    subplot(3,7,v)
    hold on
    text(-0.2, 1.1, labelstring(v)', 'Units', 'normalized', 'FontWeight', 'bold','FontSize', 10)
end

% widthInches = 35;
% heightInches = 15;
% run('ScriptForExportingImages.m')

%%%%%%%%%%%%%%%% figure 10 %%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(10);
figname = 'figure10';

figure_count = 1;
load("Without_DDS_changing_kD.mat")


% dose response VEGF
subplot(2,3,figure_count);
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

% reference levels
y10 = 0.10 * v_ret_Initial;
y50 = 0.50 * v_ret_Initial;
% draw horizontal lines on the RIGHT axis
h10 = yline(y10, '--', '', ...
    'LabelHorizontalAlignment','right', 'LabelVerticalAlignment','bottom', ...
    'HandleVisibility','off', 'LineWidth',2.2);
h50 = yline(y50, ':',  '', ...
    'LabelHorizontalAlignment','right', 'LabelVerticalAlignment','bottom', ...
    'HandleVisibility','off', 'LineWidth',2.2);

axis square
box on
figure_count = figure_count+1;
hold off

subplot(2,3,figure_count);
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

% reference levels
y10 = 0.10 * v_vit_Initial;
y50 = 0.50 * v_vit_Initial;
% draw horizontal lines on the RIGHT axis
h10 = yline(y10, '--', '10% of baseline', ...
    'LabelHorizontalAlignment','right', 'LabelVerticalAlignment','bottom', ...
    'HandleVisibility','off', 'LineWidth',2.2);
h50 = yline(y50, ':',  '50% of baseline', ...
    'LabelHorizontalAlignment','right', 'LabelVerticalAlignment','bottom', ...
    'HandleVisibility','off', 'LineWidth',2.2);

axis square
box on
figure_count = figure_count+1;
hold off

subplot(2,3,figure_count);
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

% reference levels
y10 = 0.10 * v_aq_Initial;
y50 = 0.50 * v_aq_Initial;
% draw horizontal lines on the RIGHT axis
h10 = yline(y10, '--', '', ...
    'LabelHorizontalAlignment','right', 'LabelVerticalAlignment','bottom', ...
    'HandleVisibility','off', 'LineWidth',2.2);
h50 = yline(y50, ':',  '', ...
    'LabelHorizontalAlignment','right', 'LabelVerticalAlignment','bottom', ...
    'HandleVisibility','off', 'LineWidth',2.2);

axis square
box on
figure_count = figure_count+1;
hold off

load("With_DDS_changing_kD.mat")

% dose response VEGF
subplot(2,3,figure_count);
hold on
for j = 1:length(dose_in)
plot(t,C_vret_Data(:,j), 'LineWidth', 2)
end
xlabel('Time(days)')
ylabel('Free VEGF (pM)')
title('Retina with DDS')
ylim([-0.6,60])
xlim([-3 300])
set(gca, 'FontSize', 20)

% reference levels
y10 = 0.10 * v_ret_Initial;
y50 = 0.50 * v_ret_Initial;
% draw horizontal lines on the RIGHT axis
h10 = yline(y10, '--', '', ...
    'LabelHorizontalAlignment','right', 'LabelVerticalAlignment','bottom', ...
    'HandleVisibility','off', 'LineWidth',2.2);
h50 = yline(y50, ':',  '', ...
    'LabelHorizontalAlignment','right', 'LabelVerticalAlignment','bottom', ...
    'HandleVisibility','off', 'LineWidth',2.2);

axis square
box on
figure_count = figure_count+1;
hold off

subplot(2,3,figure_count);
hold on
for j = 1:length(dose_in)
plot(t,C_vvit_Data(:,j), 'LineWidth', 2)
end
xlabel('Time(days)')
ylabel('Free VEGF (pM)')
title('Vitreous with DDS')
ylim([-0.2,20])
xlim([-3 300])
set(gca, 'FontSize', 20)

% reference levels
y10 = 0.10 * v_vit_Initial;
y50 = 0.50 * v_vit_Initial;
% draw horizontal lines on the RIGHT axis
h10 = yline(y10, '--', '', ...
    'LabelHorizontalAlignment','right', 'LabelVerticalAlignment','bottom', ...
    'HandleVisibility','off', 'LineWidth',2.2);
h50 = yline(y50, ':',  '', ...
    'LabelHorizontalAlignment','right', 'LabelVerticalAlignment','bottom', ...
    'HandleVisibility','off', 'LineWidth',2.2);

axis square
box on
figure_count = figure_count+1;
hold off

subplot(2,3,figure_count);
hold on
for j = 1:length(dose_in)
plot(t,C_vaq_Data(:,j), 'LineWidth', 2)
end
xlabel('Time(days)')
ylabel('Free VEGF (pM)')
title('Aqueous with DDS')
ylim([-0.018,1.8])
xlim([-3 300])
set(gca, 'FontSize', 20)

% reference levels
y10 = 0.10 * v_aq_Initial;
y50 = 0.50 * v_aq_Initial;
% draw horizontal lines on the RIGHT axis
h10 = yline(y10, '--', '', ...
    'LabelHorizontalAlignment','right', 'LabelVerticalAlignment','bottom', ...
    'HandleVisibility','off', 'LineWidth',2.2);
h50 = yline(y50, ':',  '', ...
    'LabelHorizontalAlignment','right', 'LabelVerticalAlignment','bottom', ...
    'HandleVisibility','off', 'LineWidth',2.2);

axis square
box on
figure_count = figure_count+1;
hold off

labelstring = {'A', 'B', 'C', 'D','E', 'F'};
for v = 1:6
    subplot(2,3,v)
    hold on
    text(-0.2, 1.1, labelstring(v)', 'Units', 'normalized', 'FontWeight', 'bold','FontSize', 14)
end

% widthInches = 15;
% heightInches = 10;
% run('ScriptForExportingImages.m')




%%%%%%%%%%%%%%%% figure 11 %%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(11);
figname = 'figure11';

figure_count = 1;



load("Without_DDS_changing_koff.mat")
% dose response VEGF
subplot(2,3,figure_count);
hold on
for j = 1:length(dose_in)
plot(t,C_vret_Data(:,j), 'LineWidth', 2)
end
xlabel('Time(days)')
ylabel('Free VEGF (pM)')
title('Retina')
legend('k_{off} 0.345', 'k_{off} 0.864', 'k_{off} 1.728', 'k_{off} 3.456', 'FontSize',12, 'Location','southeast');
ylim([-0.6,60])
xlim([-3 300])
set(gca, 'FontSize', 20)
axis square
box on
figure_count = figure_count+1;
hold off

subplot(2,3,figure_count);
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
figure_count = figure_count+1;
hold off

subplot(2,3,figure_count);
hold on
for j = 1:length(dose_in)
plot(t,C_vaq_Data(:,j), 'LineWidth', 2)
end
xlabel('Time(days)')
ylabel('Free VEGF (pM)')
title('Aqueous')
ylim([-0.018,1.8])
xlim([-3 300])
set(gca, 'FontSize', 20)
axis square
box on
figure_count = figure_count+1;
hold off




load("With_DDS_changing_Koff.mat")


% dose response VEGF
subplot(2,3,figure_count);
hold on
for j = 1:length(dose_in)
plot(t,C_vret_Data(:,j), 'LineWidth', 2)
end
xlabel('Time(days)')
ylabel('Free VEGF (pM)')
title('Retina with DDS')
ylim([-0.6,60])
xlim([-3 300])
set(gca, 'FontSize', 20)
axis square
box on
figure_count = figure_count+1;
hold off

subplot(2,3,figure_count);
hold on
for j = 1:length(dose_in)
plot(t,C_vvit_Data(:,j), 'LineWidth', 2)
end
xlabel('Time(days)')
ylabel('Free VEGF (pM)')
title('Vitreous with DDS')
ylim([-0.2,20])
xlim([-3 300])
set(gca, 'FontSize', 20)
axis square
box on
figure_count = figure_count+1;
hold off

subplot(2,3,figure_count);
hold on
for j = 1:length(dose_in)
plot(t,C_vaq_Data(:,j), 'LineWidth', 2)
end
xlabel('Time(days)')
ylabel('Free VEGF (pM)')
title('Aqueous with DDS')
%legend('Dose 0.05 mg', 'Dose 0.1 mg', 'Dose 0.2 mg', 'Dose 0.5 mg', 'Dose 1 mg', 'Dose 2 mg', 'FontSize',12, 'Location','southeast');
ylim([-0.018,1.8])
xlim([-3 300])
set(gca, 'FontSize', 20)
axis square
box on
figure_count = figure_count+1;
hold off



labelstring = {'A', 'B', 'C', 'D','E', 'F'};
for v = 1:6
    subplot(2,3,v)
    hold on
    text(-0.2, 1.1, labelstring(v)', 'Units', 'normalized', 'FontWeight', 'bold','FontSize', 14)
end

% widthInches = 15;
% heightInches = 10;
% 
% 
% run('ScriptForExportingImages.m')



%%%%%%%%%%%%%%%% figure 12 %%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(12);
figname = 'figure12';

figure_count = 1;

load("Without_DDS_changing_kD.mat")
% bar plot 10%
barWidth = 1;
subplot(2,4,figure_count);
hold on
hBar1 = bar([Data_time_at_target_ret_10', Data_time_at_target_vit_10',Data_time_at_target_aq_10'], 'grouped', 'LineWidth', 1.5); 
xticks(1:length(dose_in)); 
xticklabels({'19000', '9500', '4750', '1900'});
ylabel('Pharmacodynamic Suppression Time (Days)');
xlabel('k_{D}');
set(hBar1, 'BarWidth', barWidth);
set(gca, 'FontSize', 11)
legend({'Retina 10%', 'Vitreous 10%', 'Aqueous 10%'}, 'FontSize',14, 'Location', 'northwest');
ylim([0,700])
pbaspect([1 1 1])
axis square
box on
figure_count = figure_count+1;
hold off

% bar plot 50%
barWidth = 1;
subplot(2,4,figure_count);
hold on
hBar1 = bar([Data_time_at_target_ret_50', Data_time_at_target_vit_50',Data_time_at_target_aq_50'], 'grouped', 'LineWidth', 1.5); 
xticks(1:length(dose_in)); 
xticklabels({'19000', '9500', '4750', '1900'});
ylabel('Pharmacodynamic Suppression Time (Days)');
xlabel('k_{D}');
set(hBar1, 'BarWidth', barWidth);
set(gca, 'FontSize', 11)
legend({'Retina 50%', 'Vitreous 50%', 'Aqueous 50%'}, 'FontSize',14, 'Location', 'northwest');
ylim([0,700])
pbaspect([1 1 1])
axis square
box on
figure_count = figure_count+1;
hold off

load("Without_DDS_changing_koff.mat")
% bar plot 10%
barWidth = 1;
subplot(2,4,figure_count);
hold on
hBar1 = bar([Data_time_at_target_ret_10', Data_time_at_target_vit_10',Data_time_at_target_aq_10'], 'grouped', 'LineWidth', 1.5); 
xticks(1:length(dose_in)); 
xticklabels({'0.345', '0.864', '1.728', '3.456'});
ylabel('Pharmacodynamic Suppression Time (Days)');
xlabel('k_{off}');
set(hBar1, 'BarWidth', barWidth);
set(gca, 'FontSize', 11)
legend({'Retina 10%', 'Vitreous 10%', 'Aqueous 10%'}, 'FontSize',14, 'Location', 'northwest');
ylim([0,700])
pbaspect([1 1 1])
axis square
box on
figure_count = figure_count+1;
hold off

% bar plot 50%
barWidth = 1;
subplot(2,4,figure_count);
hold on
hBar1 = bar([Data_time_at_target_ret_50', Data_time_at_target_vit_50',Data_time_at_target_aq_50'], 'grouped', 'LineWidth', 1.5); 
xticks(1:length(dose_in)); 
xticklabels({'0.345', '0.864', '1.728', '3.456'});
ylabel('Pharmacodynamic Suppression Time (Days)');
xlabel('k_{off}');
set(hBar1, 'BarWidth', barWidth);
set(gca, 'FontSize', 11)
legend({'Retina 50%', 'Vitreous 50%', 'Aqueous 50%'}, 'FontSize',14, 'Location', 'northwest');
ylim([0,700])
pbaspect([1 1 1])
axis square
box on
figure_count = figure_count+1;
hold off



load("With_DDS_changing_kD.mat")
% bar plot 10%
barWidth = 1;
subplot(2,4,figure_count);
hold on
hBar1 = bar([Data_time_at_target_ret_10', Data_time_at_target_vit_10',Data_time_at_target_aq_10'], 'grouped', 'LineWidth', 1.5); 
xticks(1:length(dose_in)); 
xticklabels({'19000', '9500', '4750', '1900'});
ylabel('Pharmacodynamic Suppression Time (Days)');
xlabel('k_{D}');
set(hBar1, 'BarWidth', barWidth);
set(gca, 'FontSize', 11)
legend({'Retina 10%', 'Vitreous 10%', 'Aqueous 10%'}, 'FontSize',14, 'Location', 'northwest');
ylim([0,700])
pbaspect([1 1 1])
axis square
box on
figure_count = figure_count+1;
hold off

% bar plot 50%
barWidth = 1;
subplot(2,4,figure_count);
hold on
hBar1 = bar([Data_time_at_target_ret_50', Data_time_at_target_vit_50',Data_time_at_target_aq_50'], 'grouped', 'LineWidth', 1.5); 
xticks(1:length(dose_in)); 
xticklabels({'19000', '9500', '4750', '1900'});
ylabel('Pharmacodynamic Suppression Time (Days)');
xlabel('k_{D}');
set(hBar1, 'BarWidth', barWidth);
set(gca, 'FontSize', 11)
legend({'Retina 50%', 'Vitreous 50%', 'Aqueous 50%'}, 'FontSize',14, 'Location', 'northwest');
ylim([0,700])
pbaspect([1 1 1])
axis square
box on
figure_count = figure_count+1;
hold off




load("With_DDS_changing_Koff.mat")

% bar plot 10%
barWidth = 1;
subplot(2,4,figure_count);
hold on
hBar1 = bar([Data_time_at_target_ret_10', Data_time_at_target_vit_10',Data_time_at_target_aq_10'], 'grouped', 'LineWidth', 1.5); 
xticks(1:length(dose_in)); 
xticklabels({'0.345', '0.864', '1.728', '3.456'});
ylabel('Pharmacodynamic Suppression Time (Days)');
xlabel('k_{off}');
set(hBar1, 'BarWidth', barWidth);
set(gca, 'FontSize', 11)
legend({'Retina 10%', 'Vitreous 10%', 'Aqueous 10%'}, 'FontSize',14, 'Location', 'northwest');
ylim([0,700])
pbaspect([1 1 1])
axis square
box on
figure_count = figure_count+1;
hold off

% bar plot 50%
barWidth = 1;
subplot(2,4,figure_count);
hold on
hBar1 = bar([Data_time_at_target_ret_50', Data_time_at_target_vit_50',Data_time_at_target_aq_50'], 'grouped', 'LineWidth', 1.5); 
xticks(1:length(dose_in)); 
xticklabels({'0.345', '0.864', '1.728', '3.456'});
ylabel('Pharmacodynamic Suppression Time (Days)');
xlabel('k_{off}');
set(hBar1, 'BarWidth', barWidth);
set(gca, 'FontSize', 11)
legend({'Retina 50%', 'Vitreous 50%', 'Aqueous 50%'}, 'FontSize',14, 'Location', 'northwest');
ylim([0,700])
pbaspect([1 1 1])
axis square
box on
figure_count = figure_count+1;
hold off



labelstring = {'A', 'B', 'C', 'D','E', 'F','G','H'};
for v = 1:8
    subplot(2,4,v)
    hold on
    text(-0.2, 1.1, labelstring(v)', 'Units', 'normalized', 'FontWeight', 'bold','FontSize', 14)
end

% widthInches = 20;
% heightInches = 10;
% Add row titles
annotation('textbox', [0.05, 0.94, 0.9, 0.05], 'String', 'Without DDS', ...
    'FontSize', 16, 'FontWeight', 'bold', 'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center');

annotation('textbox', [0.05, 0.47, 0.9, 0.05], 'String', 'With DDS', ...
    'FontSize', 16, 'FontWeight', 'bold', 'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center');

%run('ScriptForExportingImages.m')
