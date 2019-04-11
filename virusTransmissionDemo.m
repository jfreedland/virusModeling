%% Virus transmission demo
% By J. Freedland, 2019.
%
% Run this script to model HIV, Hepatitis B, and Dengue virus.
% These models solve a set of three ordinary differential equations using
% published parameters derived from clinical studies.
%
% Our first model (HIV) simulates viral spread (virusTransmission.m)
% Then, it calculates how much more infective HIV would need to be to
% successfully infect ~100% of cells within 8 hours (findInfectivity.m)
%
% Our second model (Hepatitis B) simulates viral spread (virusTransmission.m)
% Then, it calculates how much more infective Hepatitis B would need to be
% successfully infect several different percentages of cells over 8 hours
% (infectivityConstants.m)
%
% Our third model (Dengue virus) only simulates viral spread. (virusTransmission.m)
%%%%%%%%%%%%%%%%%%%

%% HIV model
% source: https://www.mdpi.com/1999-4915/4/10/1984/htm

% Amount of solution considered (in mL)
mL = 27;

% Initial concentrations of cells.
params.target      = 10 * (mL*1e3);     % Healthy cells per uL
params.infected    = 0 * (mL*1e3);      % Infected cells per uL
params.virus       = 10^-6 * (mL*1e3);  % Viruses per uL

% Virus-specific spread parameters.
params.lambda  = 0.17 * (mL*1e3);       % Production of new healthy cells (per uL)
params.dT      = 0;                     % Healthy cell death rate (without virus)
params.dI      = 0.39;                  % Infected cell death rate
params.beta    = 6.5*10^-4 / (mL*1e3);  % Virus infectivity (per uL^-1)
params.p       = 850;                   % Virus production rate
params.c       = 3;                     % Virus clearence rate (by immune system)

params.time_phase   = 500;      % Number of hours/days to consider (depends on units)

% Model results for HIV virus.
[t, T, I, V]   = virusTransmission(params);

% Variant virus.
percentComplete     = 100;  % Percent of infection completed in ...
hours               = 8;    % ... ___ hours.
params.time_phase   = 1;    % Number of hours/days to model through (depends on units).
params.dI           = 0;    % Assume infected cells are immortal.
stableBeta          = params.beta;

% How much more infective does the new virus need to be?
% ex: if infectivityConstant = 2, the new virus must be 2x as infective.
infectivityConstant = findInfectivity(params,percentComplete,hours);

% Model using the new infectivity.
params.beta         = stableBeta * infectivityConstant;
[t2, T2, I2, V2]    = virusTransmission(params);
t2 = t2 * 24;

subplot(1,2,1)
semilogy(t,T,t,I,t,V,'LineWidth',3)
legend('healthy','infected','virus')
xlabel('time (days)')
ylabel('number of cells')
title('Standard model')

subplot(1,2,2)
semilogy(t2,T2,t2,I2,t2,V2,'LineWidth',3)
legend('healthy','infected','virus')
xlabel('time (hours)')
ylabel('number of cells')
title('Adjusted Model')

% For exporting:
A = [t', I'];
B = [t2', I2'];

%% Hepatitis B model
% source: https://www.sciencedirect.com/science/article/pii/S0022519307000938

% Amount of solution considered (in ml)
mL = 27;

% Initial concentrations of cells.
params.target      = (13.6*10^6) * mL;  % Healthy cells per mL
params.infected    = 0 * mL;            % Infected cells per mL
params.virus       = 0.33 * mL;         % Viruses per mL

% Values of model parameters passed in y
params.lambda  = 10 * mL;               % per mL
params.dT      = 0;
params.dI      = 0.5;
params.beta    = 7.57*10^-10 / mL;      % per mL^-1
params.p       = 102;
params.c       = 0.67;

params.time_phase   = 200;
[t, T, I, V] = virusTransmission(params);

% Variant virus.
percentRange        = 50:25:100;     % Percent of infection completed in ...
hours               = 8;            % ... ___ hours.
params.time_phase   = 1;            % Number of hours/days to model through (depends on units).
params.dI           = 0;            % Assume infected cells are immortal.

% How much more infective does the new virus need to be?
% ex: if infectivityConstant = 2, the new virus must be 2x as infective.
[p, infectivityConstant] = infectivityConstants(params,percentRange,hours);

subplot(1,2,1)
semilogy(t,T,t,I,t,V,'LineWidth',3)
legend('healthy','infected','virus')
xlabel('time (days)')
ylabel('number of cells')
title('Standard model')

subplot(1,2,2)
plot(p,infectivityConstant,'-o')
xlabel(['percent completed after ' num2str(hours) ' hours'])
ylabel('fold infectivity')
title('Adjusted Model')

% for exporting:
A = [t', I'];
B = [p', infectivityConstant'];

%% Dengue Virus
% source: https://www.sciencedirect.com/science/article/pii/S0895717708002732

% Amount of solution considered (in ml)
mL = 27;

% Initial concentrations of cells.
params.target      = 10 * (mL*1e3);     % Healthy cells per uL
params.infected    = 0 * (mL*1e3);      % Infected cells per uL
params.virus       = 10^-6 * (mL*1e3);  % Viruses per uL

% Values of model parameters passed in y
params.lambda  = 80 * (mL*1e3);         % per uL
params.dT      = 0;
params.dI      = 0.5; 
params.beta    = (1*10^-3) / (mL*1e3);  % per uL^-1
params.p       = 20;
params.c       = 0.8;

params.time_phase   = 15;
[t, T, I, V] = virusTransmission(params);

semilogy(t,T,t,I,t,V,'LineWidth',3)
legend('healthy','infected','virus')
xlabel('time (days)')
ylabel('number of cells')
title('Standard model')

% for exporting:
A = [t', I'];