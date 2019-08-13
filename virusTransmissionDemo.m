%% README
% By J. Freedland, 2019.

% This demo simulates the spread of HIV and dengue viruses within a host.
% It relies on two main functions:
%
%   modelSystem:        solves a system of ODEs to simulate the spread of the
%                       respective bloodborne virus.
%   virusAdaptation:    uses fminsearch to "adapt" our virus until it is
%                       capable of spreading within a certain period of time.
%
% A series of initial conditions are presented below based on
% referenced literature.

%% HIV (without adaptation)
% Modeling of the spread of HIV in a host.
% Source: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5072357/

% Initial concentrations of cells.
params.target      = 5e3;       % Healthy cells per uL
params.infected    = 0;         % Infected cells per uL
params.virus       = 0.4e-3;    % Viruses per uL

% Virus-specific spread parameters.
params.lambda  = 100;                   % Production of new healthy cells (per uL)
params.dT      = 0.1;                   % Healthy cell death rate (without virus)
params.dI      = 0.5;                   % Infected cell death rate
params.beta    = 10^-5;                 % Virus infectivity (per uL^-1)
params.p       = 1.5*10^3;              % Virus production rate
params.c       = 10;                    % Virus clearence rate (by immune system)
params.omega   = 10^-3;                 % Cell-to-cell infectivity (can set to 0 to consider virus-cell interactions only).

params.time_phase   = 30;      % Number of hours/days to consider (depends on units)
f0 = [1 1 1];                  % Adaptational parameters ([1 1 1] == no adaptations).

[t, T, I, V] = modelSystem(params,f0);  % Solve set of ODEs.

figure(1)
semilogy(t,T,t,I,t,V) 
xlabel('time')
ylabel('population')
legend({'healthy cells','infected cells','viruses'})

%% HIV (with adaptation)
% Modeling of the spread of an adapted virus in a host. 
% Source: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5072357/

% Interpretting our output:
% f0 is the output of our function. It is a set of three values ([x y z])
% that adjust the following parameters:
    % Virus infectivity:            x * params.beta;
    % Virus replication:            y * params.p;
    % Cell-to-cell infectivity:     z * params.omega;
% For instance, an output of [2 1 1] suggests that a 2x increase in virus
% infectivity can create a "successful infection" in the allotted time.

% Each simulation will produce different solutions, all of which are valid.
% Our code deliberately biases a specific pathway to ensure diversity in our output.

% Initial concentrations of cells.
params.target      = 5e3;       % Healthy cells per uL
params.infected    = 0;         % Infected cells per uL
params.virus       = 0.4e-3;    % Viruses per uL

% Virus-specific spread parameters.
params.lambda  = 100;                   % Production of new healthy cells (per uL)
params.dT      = 0.1;                   % Healthy cell death rate (without virus)
params.dI      = 0.5;                   % Infected cell death rate
params.beta    = 10^-5;                 % Virus infectivity (per uL^-1)
params.p       = 1.5*10^3;              % Virus production rate
params.c       = 10;                    % Virus clearence rate (by immune system)
params.omega   = 10^-3;                 % Cell-to-cell infectivity (can set to 0 to consider virus-cell interactions only).

params.time_phase   = 30;      % Number of hours/days to consider in the base simulation (depends on units).
nSimulations        = 1;       % Number of adaptations to run
time_optimized      = 1;       % Number of hours/days for our virus to fully infects in (depends on units)
percentInfected     = 1;       % Percent of healthy cells infected for us to conclude a "successful infection" occurred. (0 to 1)

fSolution = zeros(nSimulations,3);
for a = 1:nSimulations
    [f0, ~] = virusAdaptation(params,time_optimized,percentInfected); % Perform adaptation search.
    fSolution(a,:) = f0; % Save result.
    disp(f0)
end

% Simulate result
[t1, T1, I1, V1] = modelSystem(params,[1 1 1]); % Original virus.
[t2, T2, I2, V2] = modelSystem(params,f0);      % New virus

figure(1)
subplot(1,2,1)
semilogy(t1,T1,t1,I1,t1,V1) 
xlabel('time')
ylabel('population')
title('original virus')

subplot(1,2,2)
semilogy(t2,T2,t2,I2,t2,V2) 
legend({'healthy cells','infected cells','viruses'})
title('adapted virus')
%% Dengue Virus
% source: https://www.sciencedirect.com/science/article/pii/S0895717708002732

% Amount of solution considered (in ml)
mL = 27;

% Initial concentrations of cells.
params.target      = 5e3;     % Healthy cells per uL
params.infected    = 0;      % Infected cells per uL
params.virus       = 0.4e-3;  % Viruses per uL

% Values of model parameters passed in y
params.lambda  = 80;         % per uL
params.dT      = 0;
params.dI      = 0.5; 
params.beta    = 10^-3;  % per uL^-1
params.p       = 20;
params.c       = 0.8;
params.omega   = 10^-3;

params.time_phase   = 30;