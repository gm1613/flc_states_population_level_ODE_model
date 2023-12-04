%% Simulations of full model for 10d + 2w cold

clc
clear
close all

tc0 = 10; % Pre-cold growth duration (days)
tcf = tc0 + 14; % Last timepoint of cold conditions (days)

fad0 = 0.75; %Initial fraction of ON dividingcopies (0.8 for vrn2, 0.8 for clf, 0.75 for ColFRI)
fnd0 = 0;    %Initial fraction of nucleated dividing copies
fsd0 = 1-fad0-fnd0; %Initial fraction of spread dividing copies
x0 = [fad0;fnd0;fsd0;0;0;0;0;0;0;0;0]; %Define initial conditions for ODE solver

tspan =0:0.01:tcf; % Array of timepoints for ODE solver (days)

options = odeset('AbsTol',1e-8,'RelTol',1e-8); %Setting tolerance parameters for solver

tex_flag = 0; % Flag set to 0 for wild-type and 1 for COOLAIR defective mutant
fluct_flag = 0; % Flag set to 0 for constant 5C conditions, 1 for FM, 2 for FS
vrn2_flag = 0; % Flag set to 0 for wild-type and 1 for a nucleation mutant

% Simulating full model using ode15s
[t,fvals] = ode15s(@(t,v) flc_state_model_fluct(t,v,tc0,tcf,tex_flag,fluct_flag,vrn2_flag), tspan,x0,options);



%% Process ODE solution to compute levels of histone modifications over time

% Extract individual model variables
VIN3 = fvals(:,7); %spliced
vin3 = fvals(:,8); %unspliced
L = fvals(:,9);

fad = fvals(:,1);
fnd = fvals(:,2);
fsd = fvals(:,3);
fa = fvals(:,4);
fn = fvals(:,5);
fs = fvals(:,6);
k36me3 = fvals(:,10);
flc_mrna = fvals(:,11);


% Compute total fraction of copies in each state
atotal = (fa + fad) ./ (1 + fa + fs + fn );
ntotal = (fn + fnd) ./ (1 + fa + fs + fn);
stotal = (fs + fsd) ./ (1 + fa + fs + fn);

% Compute levels of H3K27me3 (n: nucleation region and b: body)
k27n = 0*atotal + 1*ntotal + 1*stotal;

k27b = 0*atotal + 0*ntotal + 1*stotal;


% Plot changes over time
figure
plot(t,atotal,t,ntotal,t,stotal,'Linewidth',2)
box on
ax = gca;
ax.LineWidth = 2;
ax.FontSize = 20;
legend('atotal','ntotal','stotal')
ylim([0 1.2])

figure
plot(t,k36me3,t,k27n,'Linewidth',2)
box on
ax = gca;
ax.LineWidth = 2;
ax.FontSize = 20;
legend('k36n','k27n')

figure
plot(t,vin3,t,flc_mrna,'Linewidth',2)
box on
ax = gca;
ax.LineWidth = 2;
ax.FontSize = 20;
legend('VIN3')

%% Simulations of full model for 10d + 6w constant cold + 30d

clc
clear
close all

tc0 = 10; % Pre-cold growth duration (days)
tcf = tc0 + 42; % Last timepoint of cold conditions (days)
tend = tcf + 30;

fad0 = 0.75; %Initial fraction of ON dividingcopies (0.8 for vrn2, 0.8 for clf, 0.75 for ColFRI)
fnd0 = 0;    %Initial fraction of nucleated dividing copies
fsd0 = 1-fad0-fnd0; %Initial fraction of spread dividing copies
x0 = [fad0;fnd0;fsd0;0;0;0;0;0;0;0;0]; %Define initial conditions for ODE solver

tspan =0:0.01:tend; % Array of timepoints for ODE solver (days)

options = odeset('AbsTol',1e-8,'RelTol',1e-8); %Setting tolerance parameters for solver

tex_flag = 0; % Flag set to 0 for wild-type and 1 for COOLAIR defective mutant
fluct_flag = 0; % Flag set to 0 for constant 5C conditions, 1 for FM, 2 for FS
vrn2_flag = 0; % Flag set to 0 for wild-type and 1 for a nucleation mutant

% Simulating full model for wild-type, COOLAIR defective, and nucleation mutants
%[t,fvals] = ode15s(@(t,v) flc_state_model_fluct(t,v,tc0,tcf,tex_flag,fluct_flag,vrn2_flag), tspan,x0,options);

% Simulating full model for spreading mutant
[t,fvals] = ode15s(@(t,v) flc_state_model_fluct_spreading_mutant(t,v,tc0,tcf,tex_flag,fluct_flag,vrn2_flag), tspan,x0,options);


%% Process ODE solution to compute levels of histone modifications over time

% Extract individual model variables
VIN3 = fvals(:,7); %spliced
vin3 = fvals(:,8); %unspliced
L = fvals(:,9);

fad = fvals(:,1);
fnd = fvals(:,2);
fsd = fvals(:,3);
fa = fvals(:,4);
fn = fvals(:,5);
fs = fvals(:,6);
k36me3 = fvals(:,10);
flc_mrna = fvals(:,11);


% Compute total fraction of copies in each state
atotal = (fa + fad) ./ (1 + fa + fs + fn );
ntotal = (fn + fnd) ./ (1 + fa + fs + fn);
stotal = (fs + fsd) ./ (1 + fa + fs + fn);

% Compute levels of H3K27me3 (n: nucleation region and b: body)
k27n = 0*atotal + 1*ntotal + 1*stotal;

k27b = 0*atotal + 0*ntotal + 1*stotal;


% Plot changes over time
figure
plot(t,atotal,t,ntotal,t,stotal,'Linewidth',2)
box on
ax = gca;
ax.LineWidth = 2;
ax.FontSize = 20;
legend('atotal','ntotal','stotal')
ylim([0 1.2])

figure
plot(t,k36me3,t,k27n,'Linewidth',2)
box on
ax = gca;
ax.LineWidth = 2;
ax.FontSize = 20;
legend('k36n','k27n')

figure
plot(t,vin3,t,flc_mrna,'Linewidth',2)
box on
ax = gca;
ax.LineWidth = 2;
ax.FontSize = 20;
legend('VIN3')









