%% FLC state model for FS

function f = flc_state_model_fluct_spreading_mutant(t,v,tc0,tcf,tex_flag,fluct_flag,vrn2_flag)
%% In addition to spreading mutant, a double mutant can be simulated by setting the following flags
% tex_flag set to 0 for wild-type and 1 for COOLAIR defective mutant
% fluct_flag set to 0 for constant 5C conditions, 1 for FM, 2 for FS
% vrn2_flag set to 0 for wild-type and 1 for a nucleation mutant
%% Setting current temperature
if (fluct_flag==2)
T = 22*(t<tc0) + fluctuating_strong(t)*(tc0<=t && t<tcf) + 22*(t>=tcf);
elseif (fluct_flag==1)
T = 22*(t<tc0) + fluctuating_mild(t)*(tc0<=t && t<tcf) + 22*(t>=tcf);    
else
T = 22*(t<tc0) + 5*(tc0<=t && t<tcf) + 22*(t>=tcf);
end

VIN3 = v(7); %spliced
vin3 = v(8); %unspliced
L = v(9);


% See Antoniou-Kourounioti et al., Cell Systems, 2018 for the LSCD model
% of VIN3 dynamics
degv = 18;   %Degradation rate constant of spliced VIN3
sv = 4.4*degv; %VIN3 'splicing rate'
degL = 0.009;  %Degradation rate of long term input to VIN3
pc1 = 0.0315;  % LSCD model parameter for Current input (C)
pc2 = 0.0300;  % LSCD model parameter for Current input (C)
pD = 2.05;     % LSCD model parameter for Diurnal input (D)
tm = 6;        % LSCD model parameter for Diurnal input (D)
ps2 = 0.007*(vrn2_flag==0)+0*(vrn2_flag==1);   % VIN3 mediated nucleation rate parameter


%% Temperature and VIN3 dependent nucleation rate constant s2
T1 = -1;
T2 = 18;
if (T1<T)&&(T<T2)
s2 = ps2*VIN3*(T-(T1))*(T2-T);
else
s2 = 0;
end
%% Long-term input to VIN3: Production rate LP
if (T<17)
    LP = 1;
else
    LP = 0;
end

%% Current input to VIN3: C(T)
if (T<=8)
    C = pc1;
elseif (T<15.4)
    C = pc1 - pc2*(T-8)/(15.4-8);
else
    C = pc1 - pc2;
end
    
%% Short-term input to VIN3: S(T)
% Tmax and resetting can be ignored for both CC and FS conditions considered here
if (T<15)
    S = 1;
else
    S = 0.749;
end

%% Diurnal input to VIN3: D(t)
D = (pD + sin(2*pi*(t-(tm-1)/24)))^2;



%% Growth parameters (divisions per day)
g0 = 0.4; % Pre-cold growth
gc = 0.01; % Growth during cold
gp = 0.4;  % Post-cold growth

if (T>=15)
    g = g0;
else
    g = gc;
end

dn = 32; %See model description for details

%% Reactivation parameters (see model description for details)
clf = 0.0025; %(called 'delta' in model equations)
rct = 0.02; %(called 'gamma' in model equations)
%% The following parameter values are only valid for rct = 0.02 and clf = 0.0025
nfrac = 0.8937;
sfrac = 0.0124;
afrac = 0.0939;

%% Setting current values of remaining model variables (see model description for details)

fad = v(1);
fnd = v(2);
fsd = v(3);
fa = v(4);
fn = v(5);
fs = v(6);
k36me3 = v(10);
flc_mrna = v(11);

%% Compute fractions of active, nucleated and stpread copies
atotal = (fa + fad) ./ (1 + fa + fs + fn);
ntotal = (fn + fnd) ./ (1 + fa + fs + fn);
stotal = (fs + fsd) ./ (1 + fa + fs + fn);

% Parameters specifying contribution of each state to H3K36me3 levels
ra0 = 1;
rn = 0.3;
rs = 0;

% Parameters specifying dynamics of antisense pathway (see model
% description)
beta = 0.08;
beta2 = 0.12;
dfun1 = @(beta,t) 0.5 + 0.5*exp(-beta*(t-tc0));

% Compute temperature and time dependent silencing effect of antisense pathway
if (T<4 && tex_flag==0)
    qt = 0.05;
elseif (T>15 || tex_flag==1)
    qt = 1;
else
    qt = 1*(t<tc0) + (1 - (1-dfun1(beta,tcf))*exp(-beta2*(t-tcf))).*(t>=tcf) + dfun1(beta,t).*((tc0<=t).*( t<tcf));
end

% Compute temperature dependent Pol II elongation rate
if (T==5 && fluct_flag==0)
    vt = 0.6;
elseif(T<15 && fluct_flag>0)
    vt = 0.6;
else
    vt = 1;
end

kd_k36me3 = 1.21;%Turnover rate constant of H3K36me3
kd_flc = 2.77; %Turnover rate constant of FLC mRNA
%Compute normalised rate of co-transcriptional addition ofH3K36me3
k36me3_addition = kd_k36me3*(((ra0*qt).* atotal) + ((rn*qt).* ntotal) + ((rs).* stotal))./vt;
%Compute normalised initiation rate of FLC transcription
flc_initiation = kd_flc*(((ra0*qt).* atotal) + ((rn*qt).* ntotal) + ((rs).* stotal));    
 
%% Compute time derivatives using model equations
dfad = -s2*fad + rct*g*fnd;
dfnd = s2*fad - (clf+rct)*g*fnd;
dfsd = clf*g*fnd;
dfa = dn*g*fad - s2*fa + dn*rct*g*fnd + afrac*dn*(1-clf-rct)*g*fnd;
dfn = s2*fa + nfrac*dn*(1-clf-rct)*g*fnd;
dfs = dn*g*fsd + dn*clf*g*fnd + sfrac*dn*(1-clf-rct)*g*fnd;
dVIN3 = sv*vin3 - degv*VIN3;
dvin3 = L*S*C*D - sv*vin3;
dL = LP - degL*L;
dk36me3 = k36me3_addition - kd_k36me3*k36me3;
dflc_mrna =  flc_initiation - kd_flc*flc_mrna;


%% Return array of time derivatives to ODE solver
f = [dfad;dfnd;dfsd;dfa;dfn;dfs;dVIN3;dvin3;dL;dk36me3;dflc_mrna];

end
