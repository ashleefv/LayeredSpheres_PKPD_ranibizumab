function [output1, output2, output3] = solve_FD_spheres_variable_diffusivity(dose_amount, Total_time, DDS_geometry, radius_scale, thickness_scale)

%% Solves the PDE for Fickian diffusion with variable diffusivity within a radially symmetric sphere 
% The PDE is solved numerically using method of lines, the spherical finite
% difference discretization scheme defined and analyzed in the reference, 
% and any built-in MATLAB ODE solver. 
%% --- Reference --- 
%  A. N. Ford Versypt & R. D. Braatz, Analysis of finite difference 
%  discretization schemes for diffusion in spheres with variable 
%  diffusivity, Computers & Chemical Engineering, 71 (2014) 241-252, 
%  https://doi.org/10.1016/j.compchemeng.2014.05.022
% close all
plotting = 'yes';
InitialDrugInOuterLayer = 0;
% User defined parameters for double-walled microspheres
% variable: definition, units
% R1: radius of the inner core, cm 
% R2: total radius including the inner and outer shells, cm
    % If R1 = 0, then sphere is all outer shell material without an inner core. dr1 = 0
    % If R2 = R1, then sphere is all the inner core material without an outer shell. dr2 = 0
% burst: initial burst, %
% DInner: drug diffusion coefficient for inner shell, cm2/s
% DOuter: drug diffusion coefficient for outer shell, cm2/s

%This values are for the default case, when there is no optimization
%going on. They can be modified manually.
if strcmp(DDS_geometry,"Chitosan_PCL")
    R1 = ((10.2e-4)/2)*radius_scale;
    delR = (((12.7e-4)/2) - ((10.2e-4)/2))*thickness_scale;
    R2 = (R1+delR); %If PCL/Chitosan
    DInner = 3e-15; % cm2/s;;
    DOuter = 3e-12; % cm2/s; for chiotosan 2.91e-15, for PLC 2.98e-11
    burst = 4; % %
    k = 1; %partition coefficient
    alpha0Inner = DInner/R2^2; % scale by outer radius in dimensional form
    alpha0Outer = DOuter/R2^2; % scale by outer radius in dimensional form
end
if strcmp(DDS_geometry,"Chitosan")
    R1 = 0;
    R2 = ((10.2e-4)/2)*radius_scale; %If PCL/Chitosan
    %DInner = 2.91e-15; % cm2/s;;
    DOuter = 6.5e-15; % cm2/s; for chiotosan 2.91e-15, for PLC 2.98e-11
    burst = 27; % %
    k = 1; %partition coefficient
    %alpha0Inner = DInner/R2^2; % scale by outer radius in dimensional form
    alpha0Outer = DOuter/R2^2; % scale by outer radius in dimensional form
end
if strcmp(DDS_geometry,"PCL")
    R1 = 0;
    R2 = ((12.7e-4)/2)*radius_scale; %If PCL/Chitosan
    %DInner = 2.91e-15; % cm2/s;;
    DOuter = 3e-12; % cm2/s; for chiotosan 2.91e-15, for PLC 2.98e-11
    burst = 4; % %
    k = 1; %partition coefficient
    %alpha0Inner = DInner/R2^2; % scale by outer radius in dimensional form
    alpha0Outer = DOuter/R2^2; % scale by outer radius in dimensional form
end


%% Diffusivity cases (defined computationally in FD_spheres_variable_diffusivity.m)
DCASE = 5; 
% bilayered microspheres

%% Parameters
%%% Chitosan-PCL core-shell microparticle
% Inner diameter: Chitosan monolayer particle 8.9 +/- 3.5 um
% Outer diameter: Chitosan-PCL bilayer particle 12.7 +/- 5.9 um
% Outer PCL layer 1-1.5um
Mdesired = 1000; % single layer
if R2>R1 
    if R1>0% two layers
        M = floor(R1/R2*Mdesired); %number of spatial intervals in Method of Lines for each layer
    else 
        M = Mdesired;
    end
else
    M = Mdesired;
end
NR = M+1; % number of spatial discretization points along r in each layer
initial_condition = dose_amount; % c(r,0) for 0 <= r < 1
boundary_condition = 0; % c(1,t) for t >= 0

% dr = 1/M; %\Delta r
% r = 0:dr:1; % dimensionless radius
%% Unit tests: with R2 = 1, R1 = 0 and 1 with M = 200 should give the same results as R1 = 0.5 with M = 100 if DInner = DOuter

if R2<R1
    beep
    'r2 must be greater than or equal to r1'
end
% initially
NR1 = NR;
NR2 = NR;
rInner = linspace(0,R1,NR1)./R2; % inner core from 0 to r1
dr1 = rInner(2)-rInner(1);

% to keep sundqvist FD approximately 2nd order
% On the Use of Nonuniform Grids in Finite-Difference Equations equation (4)

if 0<R1/R2 && R1/R2 <1 %Making sure R1 is not 0 and R1 is different to R2
    dr2approx = dr1^2+dr1;
    NR2 = ceil((R2-R1)./R2/dr2approx)+1;
end

rOuter = linspace(R1,R2,NR2)./R2;
dr2 = rOuter(end)-rOuter(end-1);
if abs(dr1-0) < eps
    r = rOuter;
elseif abs(dr2-0) < eps
    r = rInner;
else
    r = [rInner, rOuter]; % combined r for inner core and outer shell replicating r1, non-dimensionalized
end
Ntotal = length(r); % total numuber of points in two zones.
Mtotal = Ntotal - 1; % total number of intervals = 2*M
if dr1 >0 && dr2 > 0
    % two layers
    % inner BC at rr = 1, interior at rr=2:NR for inner core, 
    % interior at rr = NR+1:2*NR for outer shell
    % note for future, rr = NR will have dr not equal on both size of FD.
    c0(1:NR1-1) = initial_condition;
    c0(NR1+2:Ntotal-1) = InitialDrugInOuterLayer;
    % outer BC at 2*NR-1
    c0(Ntotal) = boundary_condition;
else
    c0(1:NR-1) = initial_condition;
    c0(NR) = boundary_condition;
end



if DCASE == 1
    timevector = [0 0.05 0.1 0.15 0.25 1]./alpha0;
    titlestring = 'Diffusion Case I';
elseif DCASE == 2
    timevector = [0 0.05 0.1175 0.14 0.15 0.2]./alpha0;
    titlestring = 'Diffusion Case II';
elseif DCASE == 3
    timevector = [0 0.05 0.1 0.25 0.5 1]./alpha0;
    titlestring = 'Diffusion Case III';
elseif DCASE == 4
    timevector = [0 0.02 0.05 0.07 0.1 0.28]./alpha0;
    titlestring = 'Diffusion Case IV';
elseif DCASE == 5
    timevector = [0:0.1:Total_time].*86400; % in seconds
    %timevector = [0:1:10].*3600; % scaled into hours.
    titlestring = 'Diffusion Case V';    
end

% CASE V: constant diffusivity with alpha = alpha0Inner 0<r<r1
            % and alpha = alpha0Outer r1<r<r2



% Make alpha a scalar for one layer and a vector of two values if two
% layers (inner and outer)
if dr1 >0 && dr2 > 0
%     r1 = (NR1-1)*dr1;
%     r2 = (NR2-1)*dr2;
%     alpha0Interface = r1*alpha0Inner+r2*alpha0Outer; % interface as the weighted average of the diffusivities in the two adjacent zones
%     alpha = [alpha0Inner, alpha0Interface, alpha0Outer];
    alpha = [alpha0Inner alpha0Outer];
    %Establish I.C at the interface
    gamma = alpha0Outer*dr1/(alpha0Inner*dr2);
    c0(NR1) = (initial_condition+InitialDrugInOuterLayer*gamma)/(1+gamma/k); %Interface condition from the core side
    c0(NR1+1) = c0(NR1)/k; %Interface condition from the shell side
elseif dr2 == 0
    alpha = alpha0Inner;
else
    alpha = alpha0Outer;   
end


%% ODE solver (ode45) call
[time,concentration] = ode23s(@(t,c)FD_spheres_variable_diffusivity_two_spheres(t,c,dr1,dr2,NR1,NR2,alpha,r,k),timevector,c0);
c=concentration';

%% drug concentration over time at different locations
% figure(1)
% plot(time/(60*60*24),c(NR1,:))
% xlabel('Time (days)')
% ylabel('Concentration')


if dr1>0
    denom1 = dr1*simps(c0(1:NR1).*rInner(1:NR1).*rInner(1:NR1));
    if dr2 >0
     	denom2 = dr2*simps(c0(NR1+1:Ntotal).*rOuter(1:NR2).*rOuter(1:NR2));         
    else
        denom2 = 0 ;  
    end
else
	denom1 = 0;
	denom2 = initial_condition/3*(rOuter(end)^3-rInner(end)^3);
end
cumulrel_num = ones(size(time));
    for tt = 1:length(time)
        if time(tt) == 0
            cumulrel_num(tt) = burst;
        else
            if dr1 >0
                % If X is a scalar spacing, then simps(X,Y) is equivalent to X*simps(Y).
                integrand1 = dr1*simps((concentration(tt,1:NR1)).*rInner(1:NR1).*rInner(1:NR1)); % from 0 to r1                
                if dr2 >0
                    integrand2 = dr2*simps((concentration(tt,NR1+1:Ntotal)).*rOuter(1:NR2).*rOuter(1:NR2)); % from r1 to r2
                else
                    integrand2 = 0;
                end
            else
                integrand1 = 0;
                integrand2 = dr2*simps((concentration(tt,1:NR)).*rOuter(1:NR).*rOuter(1:NR)); % from 0 to r2
            end  
            cumulrel_num(tt) = (100-burst)*(1-( integrand1 + integrand2 )/(denom1+denom2))+burst; % integrating from 0 to r1, then r1 to r2
        end
    end


% drug release percentages over time units of %/s
drug_release = zeros(size(cumulrel_num));
 for i = 1:length(cumulrel_num)   
    
    if i == 1
        drug_release(i) =  cumulrel_num(1); % do not use this as initial drug release RATE, instead this is the burst amount that determines initial drug dose AMOUNT
    else
        drug_release(i) = (cumulrel_num(i)-cumulrel_num(i-1)) / (time(i) - time(i-1));
    end
 end

% initial drug dose mg/s
initial_drug_dose = initial_condition*(drug_release(1)/100);

% drug release over time IC units (mg)/s
for i = 1:length(cumulrel_num) 
if i == 1
    drug_release(i) = initial_condition*(drug_release(2)/100)*60*60*24;
else
    drug_release(i) = initial_condition*(drug_release(i)/100)*60*60*24;
end
end
% figure(2)
% plot(time/(60*60*24),drug_release)
% xlabel('Time (days)')
% ylabel('Drug release')

output1 = time;
output2 = drug_release;
output3= initial_drug_dose;
end