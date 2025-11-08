clc;
clear;



%%%%%%%%%%%%%%%% figure 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%

load ("without_DDS.mat")

% --- Load patient data ---
M = readmatrix('Patient40Data.xlsx','Sheet','Sheet1');
t_pat    = M(:,1); 
vegf_pat = M(:,6);

figname = 'figure2';
figure_count = 1;
% Figure 2 for 0.5 mg
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
axis square
box on
%exportgraphics(figure(figure_count),sprintf('Retina%d.png',dose_in(j)), 'Resolution', 300)
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
axis square
box on
%exportgraphics(figure(figure_count),sprintf('Vitreous%d.png',dose_in(j)), 'Resolution', 300)
figure_count = figure_count+1;

subplot(2,3,figure_count)
yyaxis left
semilogy(t,C_raq_Data(:,j), 'LineWidth', 2)
xlabel('Time(days)')
ylabel('Ranibizumab (pM)')
title('Aqueous')
ylim([1, 1E7]);
xlim([-3 300])

yyaxis right
plot(t, C_vaq_Data(:,j), 'LineWidth', 2, 'DisplayName','Free VEGF')  
hold on
scatter(t_pat, vegf_pat, 36, 'o', ...
    'MarkerFaceColor','none', ...     
    'MarkerEdgeColor','k', ...        
    'LineWidth',1.2, ...
    'DisplayName','Patient data');        
ylabel('Free VEGF (pM)')
set(gca, 'FontSize', 16)

axis square
box on
legend({'Ranibizumab', 'Free VEGF', 'Patitent Data'}, 'FontSize',10, 'Location', 'southeast');
%exportgraphics(figure(figure_count),sprintf('Retina%d.png',dose_in(j)), 'Resolution', 300)
figure_count = figure_count+1;
hold off
end
end  

% Figure 2 for 0.5 mg
load("DDS_dose.mat")
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
axis square
box on
%exportgraphics(figure(figure_count),sprintf('Retina%d.png',dose_in(j)), 'Resolution', 300)
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
axis square
box on
%exportgraphics(figure(figure_count),sprintf('Vitreous%d.png',dose_in(j)), 'Resolution', 300)
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
axis square
box on
%exportgraphics(figure(figure_count),sprintf('Aqueous%d.png',dose_in(j)), 'Resolution', 300)
figure_count = figure_count+1;
end
end



labelstring = {'A', 'B', 'C', 'D','E','F'};
for v = 1:6
    subplot(2,3,v)
    hold on
    text(-0.2, 1.1, labelstring(v)', 'Units', 'normalized', 'FontWeight', 'bold','FontSize', 14)
end

widthInches = 15;
heightInches = 10;
run('ScriptForExportingImages.m')

%%%%%%%%%%%%%%%% figure 3 %%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(3);
figname = 'figure3';
figure_count = 1;

load ("without_DDS.mat")


% bar plot 10%
barWidth = 1;
subplot(2,2,figure_count);
hold on
hBar1 = bar([Data_time_at_target_ret_10', Data_time_at_target_vit_10',Data_time_at_target_aq_10'], 'grouped', 'LineWidth', 1.5); 
xticks(1:length(dose_in)); 
xticklabels({'0.05', '0.1', '0.5', '1', '2'});
ylabel('Pharmacodynamic Suppression Time (Days)');
xlabel('Drug Amount (mg)');
title('Without DDS')
set(hBar1, 'BarWidth', barWidth);
set(gca, 'FontSize', 12)
legend({'Retina 10%', 'Vitreous 10%', 'Aqueous 10%'}, 'FontSize',14, 'Location', 'northwest');
ylim([0,450])
pbaspect([1 1 1])
axis square
box on
%exportgraphics(figure(figure_count),sprintf('bar_dose_response_10.png'), 'Resolution', 300)
figure_count = figure_count+1;
hold off

% bar plot 50%
barWidth = 1;
subplot(2,2,figure_count);
hold on
hBar1 = bar([Data_time_at_target_ret_50', Data_time_at_target_vit_50',Data_time_at_target_aq_50'], 'grouped', 'LineWidth', 1.5); 
xticks(1:length(dose_in)); 
xticklabels({'0.05', '0.1', '0.5', '1', '2'});
ylabel('Pharmacodynamic Suppression Time (Days)');
xlabel('Drug Amount (mg)');
title('Without DDS')
set(hBar1, 'BarWidth', barWidth);
set(gca, 'FontSize', 12)
legend({'Retina 50%', 'Vitreous 50%', 'Aqueous 50%'}, 'FontSize',14, 'Location', 'northwest');
ylim([0,450])
pbaspect([1 1 1])
axis square
box on
%exportgraphics(figure(figure_count),sprintf('bar_dose_response_50.png'), 'Resolution', 300)
figure_count = figure_count+1;
hold off


% Figure 2 for 0.5 mg
load("DDS_dose.mat")
% bar plot 10%
barWidth = 1;
subplot(2,2,figure_count);
hold on
hBar1 = bar([Data_time_at_target_ret_10', Data_time_at_target_vit_10',Data_time_at_target_aq_10'], 'grouped', 'LineWidth', 1.5); 
xticks(1:length(dose_in)); 
xticklabels({'0.05', '0.1', '0.5', '1', '2'});
ylabel('Pharmacodynamic Suppression Time (Days)');
xlabel('Drug Amount (mg)');
title('With DDS')
set(hBar1, 'BarWidth', barWidth);
set(gca, 'FontSize', 12)
legend({'Retina 10%', 'Vitreous 10%', 'Aqueous 10%'}, 'FontSize',14, 'Location', 'northwest');
ylim([0 450]);
pbaspect([1 1 1])
axis square
box on
%exportgraphics(figure(figure_count),sprintf('bar_dose_response_10.png'), 'Resolution', 300)
figure_count = figure_count+1;
hold off

% bar plot 50%
barWidth = 1;
subplot(2,2,figure_count);
hold on
hBar1 = bar([Data_time_at_target_ret_50', Data_time_at_target_vit_50',Data_time_at_target_aq_50'], 'grouped', 'LineWidth', 1.5); 
xticks(1:length(dose_in)); 
xticklabels({'0.05', '0.1', '0.5', '1', '2'});
ylabel('Pharmacodynamic Suppression Time (Days)');
xlabel('Drug Amount (mg)');
title('With DDS')
set(hBar1, 'BarWidth', barWidth);
set(gca, 'FontSize', 12)
legend({'Retina 50%', 'Vitreous 50%', 'Aqueous 50%'}, 'FontSize',14, 'Location', 'northwest');
ylim([0 450]);
pbaspect([1 1 1])
axis square
box on
%exportgraphics(figure(figure_count),sprintf('bar_dose_response_50.png'), 'Resolution', 300)
figure_count = figure_count+1;
hold off



labelstring = {'A', 'B', 'C', 'D'};
for v = 1:4
    subplot(2,2,v)
    hold on
    text(-0.2, 1.1, labelstring(v)', 'Units', 'normalized', 'FontWeight', 'bold','FontSize', 14)
end

widthInches = 10;
heightInches = 10;
run('ScriptForExportingImages.m')



%%%%%%%%%%%%%%%% figure 4 %%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(4);
figname = 'figure4';
figure_count = 1;
load ("without_DDS.mat")


% dose response VEGF
subplot(2,3,figure_count);
hold on
for j = 1:length(dose_in)
plot(t,C_vret_Data(:,j), 'LineWidth', 2)
end
xlabel('Time(days)')
ylabel('Free VEGF (pM)')
title('Retina')
legend('0.05 mg', '0.1 mg', '0.5 mg', '1 mg', '2 mg', 'FontSize',12, 'Location','southeast');
ylim([-0.6,60])
xlim([-3 300])
set(gca, 'FontSize', 16)
axis square
box on
%exportgraphics(figure(figure_count),sprintf('Retina_dose_response_VEGF.png'), 'Resolution', 300)
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
set(gca, 'FontSize', 16)
axis square
box on
%exportgraphics(figure(figure_count),sprintf('Vitreous_dose_response_VEGF.png'), 'Resolution', 300)
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
set(gca, 'FontSize', 16)
axis square
box on
%exportgraphics(figure(figure_count),sprintf('Aqueous_dose_response_VEGF.png'), 'Resolution', 300)
figure_count = figure_count+1;
hold off


load("DDS_dose.mat")
for j = 1:length(dose_in)
if (dose_in(j) == 0.5)

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
set(gca, 'FontSize', 16)
axis square
box on
%exportgraphics(figure(figure_count),sprintf('Retina_dose_response_VEGF.png'), 'Resolution', 300)
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
set(gca, 'FontSize', 16)
axis square
box on
%exportgraphics(figure(figure_count),sprintf('Vitreous_dose_response_VEGF.png'), 'Resolution', 300)
figure_count = figure_count+1;
hold off

%line_style = ['-', '--', ':', '-.', '-', '--', ':', '-.','-', '--', ':', '-.'];
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
set(gca, 'FontSize', 16)
axis square
box on
%exportgraphics(figure(figure_count),sprintf('Aqueous_dose_response_VEGF.png'), 'Resolution', 300)
figure_count = figure_count+1;
hold off
end
end



labelstring = {'A', 'B', 'C', 'D','E','F'};
for v = 1:6
    subplot(2,3,v)
    hold on
    text(-0.2, 1.1, labelstring(v)', 'Units', 'normalized', 'FontWeight', 'bold','FontSize', 14)
end

widthInches = 15;
heightInches = 10;
run('ScriptForExportingImages.m')



%%%%%%%%%%%%%%%% figure 5 %%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(5);
figname = 'figure5';
figure_count = 1;
load("DDS_dose.mat")


for i = 1:length(dose_in)
    RealTime = RealTime_all{i};    
    drug_dose = drug_dose_all{i}; 
% plot_DDS_drug_release_dynamics
subplot(1,5,figure_count)
set(gca, 'FontSize', 10)
hold on
box on
plot(RealTime(2:end), drug_dose(2:end), 'LineWidth', 2)
xlim([-2 200])
xlabel('Time (days)')
%ylabel('DDS Drug Release Rate (mg/s)');
legend(sprintf('Dose %0.2f mg', dose_in(i)), 'Location', 'northeast')
pbaspect([1 1 1])
axis square
%exportgraphics(figure(figure_count),sprintf('drug_release_DDS%d.png', dose_in(i)), 'Resolution', 300)
figure_count = figure_count+1;
hold off
end




labelstring = {'A', 'B', 'C', 'D','E',};
for v = 1:5
    subplot(1,5,v)
    hold on
    text(-0.2, 1.1, labelstring(v)', 'Units', 'normalized', 'FontWeight', 'bold','FontSize', 10)
end

widthInches = 25;
heightInches = 5;
run('ScriptForExportingImages.m')


%%%%%%%%%%%%%%%% figure 6 %%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(6);
figname = 'figure6';
figure_count = 1;
load("Chitosan_single.mat")


% bar plot
barWidth = 1;
subplot(2,2,figure_count)
hold on
hBar1 = bar([Data_time_at_target_ret_10', Data_time_at_target_vit_10',Data_time_at_target_aq_10'], 'grouped', 'LineWidth', 1.5); 
xticks(1:length(radius_scale)); 
xticklabels({'0.0001R_{C}', '0.1R_{C}', '0.5R_{C}', 'R_{C}', '3R_{C}','5R_{C}', '10R_{C}'});
ylabel('Pharmacodynamic Suppression Time (Days)');
set(hBar1, 'BarWidth', barWidth);
title('Chitosan Single-Layered Varied R_{C}');
set(gca, 'FontSize', 12)
legend({'Retina 10%', 'Vitreous 10%', 'Aqueous 10%'}, 'FontSize',14, 'Location', 'northwest');
ylim([0 2000])
pbaspect([1 1 1])
axis square
box on
%exportgraphics(figure(figure_count),sprintf('barsingleC_10.png'), 'Resolution', 300)
figure_count = figure_count+1;
hold off

barWidth = 1;
subplot(2,2,figure_count)
hold on
hBar1 = bar([Data_time_at_target_ret_50', Data_time_at_target_vit_50',Data_time_at_target_aq_50'], 'grouped', 'LineWidth', 1.5); 
xticks(1:length(radius_scale)); 
xticklabels({'0.0001R_{C}', '0.1R_{C}', '0.5R_{C}', 'R_{C}', '3R_{C}','5R_{C}', '10R_{C}'});
ylabel('Pharmacodynamic Suppression Time (Days)');
set(hBar1, 'BarWidth', barWidth);
title('Chitosan Single-Layered Varied R_{C}');
set(gca, 'FontSize', 12)
%legend({'10% Suppression', '50% Suppression'},'FontSize',12, 'Location', 'northwest');
legend({'Retina 50%', 'Vitreous 50%', 'Aqueous 50%'}, 'FontSize',14, 'Location', 'northwest');
ylim([0 2000])
pbaspect([1 1 1])
axis square
box on
%exportgraphics(figure(figure_count),sprintf('barsingleC_50.png'), 'Resolution', 300)
figure_count = figure_count+1;
hold off








% bar plot
load ("PCL_single.mat")
barWidth = 1;
subplot(2,2,figure_count);
hold on
hBar1 = bar([Data_time_at_target_ret_10', Data_time_at_target_vit_10',Data_time_at_target_aq_10'], 'grouped', 'LineWidth', 1.5); 
xticks(1:length(radius_scale)); 
xticklabels({'0.0001R_{P}', 'R_{P}', '30R_{P}', '50R_{P}', '70R_{P}','90R_{P}', '100R_{P}'});
ylabel('Pharmacodynamic Suppression Time (Days)');
set(hBar1, 'BarWidth', barWidth);
title('PCL Single-Layered Varied R_{P}');
set(gca, 'FontSize', 12)
legend({'Retina 10%', 'Vitreous 10%', 'Aqueous 10%'}, 'FontSize',14, 'Location', 'northwest');
ylim([0 2000])
pbaspect([1 1 1])
axis square
box on
%exportgraphics(figure(figure_count),sprintf('barsinglePCL_10.png'), 'Resolution', 300)
figure_count = figure_count+1;
hold off

barWidth = 1;
subplot(2,2,figure_count);
hold on
hBar1 = bar([Data_time_at_target_ret_50', Data_time_at_target_vit_50',Data_time_at_target_aq_50'], 'grouped', 'LineWidth', 1.5); 
xticks(1:length(radius_scale)); 
xticklabels({'0.0001R_{P}', 'R_{P}', '30R_{P}', '50R_{P}', '70R_{P}','90R_{P}', '100R_{P}'});
ylabel('Pharmacodynamic Suppression Time (Days)');
set(hBar1, 'BarWidth', barWidth);
title('PCL Single-Layered Varied R_{P}');
set(gca, 'FontSize', 12)
%legend({'10% Suppression', '50% Suppression'},'FontSize',12, 'Location', 'northwest');
legend({'Retina 50%', 'Vitreous 50%', 'Aqueous 50%'}, 'FontSize',14, 'Location', 'northwest');
ylim([0 2000])
pbaspect([1 1 1])
axis square
box on
%exportgraphics(figure(figure_count),sprintf('barsinglePCL_50.png'), 'Resolution', 300)
figure_count = figure_count+1;
hold off


labelstring = {'A', 'B', 'C', 'D'};
for v = 1:4
    subplot(2,2,v)
    hold on
    text(-0.2, 1.1, labelstring(v)', 'Units', 'normalized', 'FontWeight', 'bold','FontSize', 14)
end

widthInches = 10;
heightInches = 10;
run('ScriptForExportingImages.m')



%%%%%%%%%%%%%%%% figure 7 %%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(7);
figname = 'figure7';

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
%legend(sprintf('%0.2fR_C, %0.2fΔR', radius_scale(i),thickness_scale(i)), 'Location', 'northeast')
%fontsize(figure(3), 17, "points")
pbaspect([1 1 1])
axis square
%exportgraphics(figure(figure_count),sprintf('drug_release_DDS%d.png', radius_scale(i)), 'Resolution', 300)
%exportgraphics(figure(figure_count),sprintf('drug_release_DDS%d.png', radius_scale(i)), 'Resolution', 300)
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
%legend(sprintf('%0.2fR_C, %0.2fΔR', radius_scale(i),thickness_scale(i)), 'Location', 'northeast')
%fontsize(figure(3), 17, "points")
pbaspect([1 1 1])
axis square
%exportgraphics(figure(figure_count),sprintf('drug_release_DDS%d.png', radius_scale(i)), 'Resolution', 300)
%exportgraphics(figure(figure_count),sprintf('drug_release_DDS%d.png', radius_scale(i)), 'Resolution', 300)
figure_count = figure_count+1;
hold off
end

labelstring = {'A', 'B', 'C', 'D','E', 'F','G','H', 'I','J','K','L','M','N'};
for v = 1:14
    subplot(2,7,v)
    hold on
    text(-0.2, 1.1, labelstring(v)', 'Units', 'normalized', 'FontWeight', 'bold','FontSize', 10)
end

widthInches = 35;
heightInches = 10;
run('ScriptForExportingImages.m')

%%%%%%%%%%%%%%%% figure 8 %%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(8);
figname = 'figure8';

figure_count = 1;
load("bi_layer_chitosan.mat")

% bar plot
barWidth = 1;
subplot(3,2,figure_count);
hold on
hBar1 = bar([Data_time_at_target_ret_10', Data_time_at_target_vit_10',Data_time_at_target_aq_10'], 'grouped', 'LineWidth', 1.5); 
xticks(1:length(radius_scale)); 
xticklabels({'0.5R_{C}', 'R_{C}', '1.5R_{C}', '2R_{C}', '3R_{C}','4R_{C}', '5R_{C}'});
ylabel('Pharmacodynamic Suppression Time (Days)');
set(hBar1, 'BarWidth', barWidth);
title('Bi-Layered Varied R_{C}, Constant \DeltaR');
set(gca, 'FontSize', 8)
legend({'Retina 10%', 'Vitreous 10%', 'Aqueous 10%'}, 'FontSize',08, 'Location', 'northwest');
ylim([0,2000])
pbaspect([1 1 1])
axis square
box on
%exportgraphics(figure(figure_count),sprintf('barbichitosan_10.png'), 'Resolution', 300)
figure_count = figure_count+1;
hold off

barWidth = 1;
subplot(3,2,figure_count);
hold on
hBar1 = bar([Data_time_at_target_ret_50', Data_time_at_target_vit_50',Data_time_at_target_aq_50'], 'grouped', 'LineWidth', 1.5); 
xticks(1:length(radius_scale)); 
xticklabels({'0.5R_{C}', 'R_{C}', '1.5R_{C}', '2R_{C}', '3R_{C}','4R_{C}', '5R_{C}'});
ylabel('Pharmacodynamic Suppression Time (Days)');
set(hBar1, 'BarWidth', barWidth);
title('Bi-Layered Varied R_{C}, Constant \DeltaR');
set(gca, 'FontSize', 8)
%legend({'10% Suppression', '50% Suppression'},'FontSize',12, 'Location', 'northwest');
legend({'Retina 50%', 'Vitreous 50%', 'Aqueous 50%'}, 'FontSize',8, 'Location', 'northwest');
ylim([0,2000])
pbaspect([1 1 1])
axis square
box on
%exportgraphics(figure(figure_count),sprintf('barbichitosan_50.png'), 'Resolution', 300)
figure_count = figure_count+1;
hold off


load("bi_layered_PCL.mat")
% bar plot
barWidth = 1;
subplot(3,2,figure_count);
hold on
hBar1 = bar([Data_time_at_target_ret_10', Data_time_at_target_vit_10',Data_time_at_target_aq_10'], 'grouped', 'LineWidth', 1.5); 
xticks(1:length(radius_scale)); 
xticklabels({'0.01\DeltaR', '0.1\DeltaR', '\DeltaR', '10\DeltaR', '20\DeltaR','25\DeltaR', '30\DeltaR'});
ylabel('Pharmacodynamic Suppression Time (Days)');
set(hBar1, 'BarWidth', barWidth);
title('Bi-Layered Varied \DeltaR, Constant R_{C}');
set(gca, 'FontSize', 08)
legend({'Retina 10%', 'Vitreous 10%', 'Aqueous 10%'}, 'FontSize',8, 'Location', 'northwest');
ylim([0,2000])
pbaspect([1 1 1])
axis square
box on
%exportgraphics(figure(figure_count),sprintf('barbiPCL_10.png'), 'Resolution', 300)
figure_count = figure_count+1;
hold off

barWidth = 1;
subplot(3,2,figure_count);
hold on
hBar1 = bar([Data_time_at_target_ret_50', Data_time_at_target_vit_50',Data_time_at_target_aq_50'], 'grouped', 'LineWidth', 1.5); 
xticks(1:length(radius_scale)); 
xticklabels({'0.01\DeltaR', '0.1\DeltaR', '\DeltaR', '10\DeltaR', '20\DeltaR','25\DeltaR', '30\DeltaR'});
ylabel('Pharmacodynamic Suppression Time (Days)');
set(hBar1, 'BarWidth', barWidth);
title('Bi-Layered Varied \DeltaR, Constant R_{C}');
set(gca, 'FontSize', 08)
%legend({'10% Suppression', '50% Suppression'},'FontSize',12, 'Location', 'northwest');
legend({'Retina 50%', 'Vitreous 50%', 'Aqueous 50%'}, 'FontSize',8, 'Location', 'northwest');
ylim([0,2000])
pbaspect([1 1 1])
axis square
box on
%exportgraphics(figure(figure_count),sprintf('barbiPCL_50.png'), 'Resolution', 300)
figure_count = figure_count+1;
hold off

load("bi_layer_changing_both.mat")
% bar plot
barWidth = 1;
subplot(3,2,figure_count);
hold on
hBar1 = bar([Data_time_at_target_ret_10', Data_time_at_target_vit_10',Data_time_at_target_aq_10'], 'grouped', 'LineWidth', 1.5); 
xticks(1:length(radius_scale)); 
xticklabels({'0.1R_{C},10\DeltaR', '0.5R_{C},2\DeltaR', '1.5R_{C},3\DeltaR', '2R_{C},5\DeltaR', '3R_{C},6\DeltaR', '4R_{C},2\DeltaR', '5R_{C},0.5\DeltaR'});
ylabel('Pharmacodynamic Suppression Time (Days)');
set(hBar1, 'BarWidth', barWidth);
title('Bi-Layered Varied R_{C} and \DeltaR');
set(gca, 'FontSize', 08)
legend({'Retina 10%', 'Vitreous 10%', 'Aqueous 10%'}, 'FontSize',8, 'Location', 'northwest');
ylim([0,2000])
pbaspect([1 1 1])
axis square
box on
%exportgraphics(figure(figure_count),sprintf('barbiPCL_10.png'), 'Resolution', 300)
figure_count = figure_count+1;
hold off

barWidth = 1;
subplot(3,2,figure_count);
hold on
hBar1 = bar([Data_time_at_target_ret_50', Data_time_at_target_vit_50',Data_time_at_target_aq_50'], 'grouped', 'LineWidth', 1.5); 
xticks(1:length(radius_scale)); 
xticklabels({'0.1R_{C},10\DeltaR', '0.5R_{C},2\DeltaR', '1.5R_{C},3\DeltaR', '2R_{C},5\DeltaR', '3R_{C},6\DeltaR', '4R_{C},2\DeltaR', '5R_{C},0.5\DeltaR'});
ylabel('Pharmacodynamic Suppression Time (Days)');
set(hBar1, 'BarWidth', barWidth);
title('Bi-Layered Varied R_{C} and \DeltaR');
set(gca, 'FontSize', 08)
%legend({'10% Suppression', '50% Suppression'},'FontSize',12, 'Location', 'northwest');
legend({'Retina 50%', 'Vitreous 50%', 'Aqueous 50%'}, 'FontSize',8, 'Location', 'northwest');
ylim([0,2000])
pbaspect([1 1 1])
axis square
box on
%exportgraphics(figure(figure_count),sprintf('barbiPCL_50.png'), 'Resolution', 300)
figure_count = figure_count+1;
hold off


labelstring = {'A', 'B', 'C', 'D','E', 'F'};
for v = 1:6
    subplot(3,2,v)
    hold on
    text(-0.2, 1.1, labelstring(v)', 'Units', 'normalized', 'FontWeight', 'bold','FontSize', 10)
end

widthInches = 10;
heightInches = 15;
run('ScriptForExportingImages.m')








%%%%%%%%%%%%%%%% figure 9 %%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(9);
figname = 'figure9';

figure_count = 1;

load("bi_layer_chitosan.mat")


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
%legend(sprintf('%0.2fR_C, %0.2fΔR', radius_scale(i),thickness_scale(i)), 'Location', 'northeast')
%fontsize(figure(3), 17, "points")
pbaspect([1 1 1])
axis square
%exportgraphics(figure(figure_count),sprintf('drug_release_DDS%d.png', radius_scale(i)), 'Resolution', 300)
%exportgraphics(figure(figure_count),sprintf('drug_release_DDS%d.png', radius_scale(i)), 'Resolution', 300)
figure_count = figure_count+1;
hold off
end


load("bi_layered_PCL.mat")
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
%legend(sprintf('%0.2fR_C, %0.2fΔR', radius_scale(i),thickness_scale(i)), 'Location', 'northeast')
%fontsize(figure(3), 17, "points")
pbaspect([1 1 1])
axis square
%exportgraphics(figure(figure_count),sprintf('drug_release_DDS%d.png', radius_scale(i)), 'Resolution', 300)
%exportgraphics(figure(figure_count),sprintf('drug_release_DDS%d.png', thickness_scale(i)), 'Resolution', 300)
figure_count = figure_count+1;
hold off

end


load("bi_layer_changing_both.mat")
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
%legend(sprintf('%0.1fΔR', thickness_scale(i)), 'Location', 'northeast')
legend(sprintf('%0.1fR_C, %0.1fΔR', radius_scale(i),thickness_scale(i)), 'Location', 'northeast')
%fontsize(figure(3), 17, "points")
pbaspect([1 1 1])
axis square
%exportgraphics(figure(figure_count),sprintf('drug_release_DDS%d.png', radius_scale(i)), 'Resolution', 300)
%exportgraphics(figure(figure_count),sprintf('drug_release_DDS%d.png', thickness_scale(i)), 'Resolution', 300)
figure_count = figure_count+1;
hold off
end

labelstring = {'A', 'B', 'C', 'D','E', 'F','G','H', 'I','J','K','L','M','N','O','P','Q','R','S','T','U'};
for v = 1:21
    subplot(3,7,v)
    hold on
    text(-0.2, 1.1, labelstring(v)', 'Units', 'normalized', 'FontWeight', 'bold','FontSize', 10)
end

widthInches = 35;
heightInches = 15;
run('ScriptForExportingImages.m')




%%%%%%%%%%%%%%%% figure 10 %%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(10);
figname = 'figure10';

figure_count = 1;
load("without_DDS_kD.mat")


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
axis square
box on
%exportgraphics(figure(figure_count),sprintf('Retina_dose_response_VEGF_nonDDS.png'), 'Resolution', 300)
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
%exportgraphics(figure(figure_count),sprintf('Vitreous_dose_response_VEGF_nonDDS.png'), 'Resolution', 300)
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
axis square
box on
%exportgraphics(figure(figure_count),sprintf('Aqueous_dose_response_VEGF_nonDDS.png'), 'Resolution', 300)
figure_count = figure_count+1;
hold off

load("DDS_doses_kD.mat")

% dose response VEGF
subplot(2,3,figure_count);
hold on
for j = 1:length(dose_in)
plot(t,C_vret_Data(:,j), 'LineWidth', 2)
end
xlabel('Time(days)')
ylabel('Free VEGF (pM)')
title('Retina with DDS')
%legend('19000', '9500', '4750', '1900', 'FontSize',12, 'Location','southeast');
ylim([-0.6,60])
xlim([-3 300])
set(gca, 'FontSize', 20)
axis square
box on
%exportgraphics(figure(figure_count),sprintf('Retina_dose_response_VEGF.png'), 'Resolution', 300)
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
%exportgraphics(figure(figure_count),sprintf('Vitreous_dose_response_VEGF.png'), 'Resolution', 300)
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
%exportgraphics(figure(figure_count),sprintf('Aqueous_dose_response_VEGF.png'), 'Resolution', 300)
figure_count = figure_count+1;
hold off

labelstring = {'A', 'B', 'C', 'D','E', 'F'};
for v = 1:6
    subplot(2,3,v)
    hold on
    text(-0.2, 1.1, labelstring(v)', 'Units', 'normalized', 'FontWeight', 'bold','FontSize', 14)
end

widthInches = 15;
heightInches = 10;
run('ScriptForExportingImages.m')





%%%%%%%%%%%%%%%% figure 11 %%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(11);
figname = 'figure11';

figure_count = 1;



load("without_DDS_koff.mat")
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
%exportgraphics(figure(figure_count),sprintf('Retina_dose_response_VEGF_nonDDS.png'), 'Resolution', 300)
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
%exportgraphics(figure(figure_count),sprintf('Vitreous_dose_response_VEGF_nonDDS.png'), 'Resolution', 300)
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
axis square
box on
%exportgraphics(figure(figure_count),sprintf('Aqueous_dose_response_VEGF_nonDDS.png'), 'Resolution', 300)
figure_count = figure_count+1;
hold off




load("DDS_Doses_koff.mat")


% dose response VEGF
subplot(2,3,figure_count);
hold on
for j = 1:length(dose_in)
plot(t,C_vret_Data(:,j), 'LineWidth', 2)
end
xlabel('Time(days)')
ylabel('Free VEGF (pM)')
title('Retina with DDS')
%legend('19000', '9500', '4750', '1900', 'FontSize',12, 'Location','southeast');
ylim([-0.6,60])
xlim([-3 300])
set(gca, 'FontSize', 20)
axis square
box on
%exportgraphics(figure(figure_count),sprintf('Retina_dose_response_VEGF.png'), 'Resolution', 300)
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
%exportgraphics(figure(figure_count),sprintf('Vitreous_dose_response_VEGF.png'), 'Resolution', 300)
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
%exportgraphics(figure(figure_count),sprintf('Aqueous_dose_response_VEGF.png'), 'Resolution', 300)
figure_count = figure_count+1;
hold off



labelstring = {'A', 'B', 'C', 'D','E', 'F'};
for v = 1:6
    subplot(2,3,v)
    hold on
    text(-0.2, 1.1, labelstring(v)', 'Units', 'normalized', 'FontWeight', 'bold','FontSize', 14)
end

widthInches = 15;
heightInches = 10;


run('ScriptForExportingImages.m')



%%%%%%%%%%%%%%%% figure 12 %%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(12);
figname = 'figure12';

figure_count = 1;

load("without_DDS_kD.mat")
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
%exportgraphics(figure(figure_count),sprintf('bar_dose_response_non_DDS_10.png'), 'Resolution', 300)
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
%exportgraphics(figure(figure_count),sprintf('bar_dose_response_non_DDS_50.png'), 'Resolution', 300)
figure_count = figure_count+1;
hold off

load("without_DDS_koff.mat")
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
%exportgraphics(figure(figure_count),sprintf('bar_dose_response_non_DDS_10.png'), 'Resolution', 300)
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
%exportgraphics(figure(figure_count),sprintf('bar_dose_response_non_DDS_50.png'), 'Resolution', 300)
figure_count = figure_count+1;
hold off



load("DDS_doses_kD.mat")
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
%exportgraphics(figure(figure_count),sprintf('bar_dose_response_10.png'), 'Resolution', 300)
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
%exportgraphics(figure(figure_count),sprintf('bar_dose_response_50.png'), 'Resolution', 300)
figure_count = figure_count+1;
hold off




load("DDS_Doses_koff.mat")

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
%exportgraphics(figure(figure_count),sprintf('bar_dose_response_10.png'), 'Resolution', 300)
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
%exportgraphics(figure(figure_count),sprintf('bar_dose_response_50.png'), 'Resolution', 300)
figure_count = figure_count+1;
hold off



labelstring = {'A', 'B', 'C', 'D','E', 'F','G','H'};
for v = 1:8
    subplot(2,4,v)
    hold on
    text(-0.2, 1.1, labelstring(v)', 'Units', 'normalized', 'FontWeight', 'bold','FontSize', 14)
end

widthInches = 20;
heightInches = 10;
% Add row titles
annotation('textbox', [0.05, 0.94, 0.9, 0.05], 'String', 'Without DDS', ...
    'FontSize', 16, 'FontWeight', 'bold', 'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center');

annotation('textbox', [0.05, 0.47, 0.9, 0.05], 'String', 'With DDS', ...
    'FontSize', 16, 'FontWeight', 'bold', 'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center');

run('ScriptForExportingImages.m')


