function dcdt = FD_spheres_variable_diffusivity_two_spheres(t,c,dr1,dr2,NR1,NR2,alpha,r,k)
%   Finite difference discritization scheme for spheres
%   with variable diffusivity as described in the reference (scheme 1).
%
%% --- Reference --- 
%  A. N. Ford Versypt & R. D. Braatz, Analysis of finite difference 
%  discretization schemes for diffusion in spheres with variable 
%  diffusivity, Computers & Chemical Engineering, 71 (2014) 241-252, 
%  https://doi.org/10.1016/j.compchemeng.2014.05.022
% 
%% --- Input ---
% 	alpha(NR)   D(r,t)/R^2 for the species
% 	dr1         Spatial discretization size for inner core
%   dr2         Spatial discretization size for outer shell
% 	c(NR)       ODE solution vector
%	NR          Number of spatial discretizations
%% --- Output ---
%   dcdt(NR)      ODE derivative vector updated with diffusion contribution

%% Compute the diffusion term at r=0
    if dr1 == 0
        dr = dr2;
    else
        dr = dr1;
    end
    for rr = 1
     	dcdt(rr) = 6*alpha(1)*(c(rr+1)-c(rr))/dr^2;
    end
    
NR = NR1;    
%% Compute the diffusion terms at 0<r<r2
% --- Loop over all interior spatial discretizations 0<r<r2 ---
if dr1 >0 && dr2 > 0
    %    for rr = NR1 interface, in algebraic form, assumes continuity
%     rr = NR1;
%     gamma = alpha(end)*dr1/(alpha(1)*dr2);
%     c(rr) = (c(rr-1)+c(rr+1)*gamma)/(1+gamma);
    % inner layer
   for rr = 2:NR1-1
       dr = dr1;
        i = rr-1; % i indexing starts at 0 and rr indexing in MATLAB starts at 1
%       scheme 0 from the reference
        dcdt(rr) = alpha(1)/i/dr^2*...
        ((i+1)*c(rr+1)-2*i*c(rr)+(i-1)*c(rr-1));
   end 
  
	% outer layer
   for rr = NR1+1:NR1+NR2-2
       dr = dr2;
        i = rr-1; % i indexing starts at 0 and rr indexing in MATLAB starts at 1
%       scheme 0 from the reference
        dcdt(rr) = alpha(end)/i/dr^2*...
        ((i+1)*c(rr+1)-2*i*c(rr)+(i-1)*c(rr-1));
   end 
    
%       % interface for rr = NR goes here in ODE form
   for rr = NR1
       gamma = alpha(end)*dr1/(alpha(1)*dr2);
       dcdt(rr) = (dcdt(rr-1)+dcdt(rr+1)*gamma)/(k+gamma);
   end
else    
	for rr = 2:NR-1
        i = rr-1; % i indexing starts at 0 and rr indexing in MATLAB starts at 1
%       scheme 0 from the reference
        dcdt(rr) = alpha(1)/i/dr^2*...
        ((i+1)*c(rr+1)-2*i*c(rr)+(i-1)*c(rr-1));
	end 
end    

%% Constant boundary condition at surface r = R2 (outer surface sphere)  
if dr1 >0 && dr2 > 0
    for rr = NR1+NR2-1
        dcdt(rr) = 0;
    end
else
    for rr = NR
        dcdt(rr) = 0;
    end
end    
    dcdt = dcdt';
    

end