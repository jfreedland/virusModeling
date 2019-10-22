%% README
% By J. Freedland, 2019.
%
% This demo simulates the spread of dengue within a host according to
% https://www.sciencedirect.com/science/article/pii/S0895717708002732.
%
% Functions:
%   modelSystem:        solves a system of ODEs to simulate the spread of a
%                       bloodborne virus given specific parameters.
%
% Inputs: 
%   We use a global variable (type) to indicate the type of simulation to run,
%       type = 1: simulation using standard parameters.
%       type = 2: simulation with zero immune cell activity.
%       type = 3: simulation with slowly degrading immune cell activity.
%           degredationTime: time until immunity is fully degraded (hours)
%       type = 4: simulation with slowly degrading immune cell activity
%                 UNTIL we reach a time when a "cure" kicks in.
%           degredationTime: time until immunity is fully degraded (hours)
%           healthResponse: time at which recovery begins (hours).
%
% This code is intended to supplement relevant literature on Buffy the
% Vampire Slayer by simulating "vampirism".


%% Replicating results for https://www.sciencedirect.com/science/article/pii/S0895717708002732
% see: Fig. 6 and Fig. 7.
clear
global type
type = 1; %#ok<NASGU>

params.S  = 250;        % Target cells
params.I  = 10;         % Infected cells
params.V  = 150;        % Viruses
params.Z  = 2000;       % Immune cells

% Values of model parameters passed in y
params.mu      = 80;    % Target cell production: cell/(day ul)
params.alpha   = 1/3;   % Life span (1/day)
params.a       = 0.001; % Infectivity
params.beta    = 0.5;   % Infection period
params.v       = 0.001; % Immune cell efficiency
params.k       = 20;    % Virus multiplication
params.gamma   = 0.8;   % Virus clearence rate
params.c       = 0.01;  % Rate of immune cell production (due to density of infected cells)
params.d       = 0.03;  % Rate of immune cell production (due to interactions with infected cells)
params.delta   = 1/365; % 1/day
params.n       = 0.265; % Base immune cell production: cell/(day ul)

params.time_phase = 15; % Time to run simulation.

[t, T, I, V, Z] = modelSystem(params);  % Solve set of ODEs.

figure(1)
subplot(2,2,1)
axis([0 250 0 15])
plot(t,T) 
subplot(2,2,2)
axis([0 15 0 15])
plot(t,I) 
subplot(2,2,3)
axis([0 250 0 15])
plot(t,V) 
subplot(2,2,4)
axis([0 10000 0 15])
plot(t,Z) 

%% Figure 1
% We repeat the above simulation using parameters observed in Darla's
% siring (Angel, S2E9).

clear
global type
type = 1; %#ok<NASGU>

% As measured in Angel, S2E9.
params.S  = 1e6;      % Target cells
params.I  = 5;        % Infected cells
params.V  = 5;        % Viruses
params.Z  = 2000;     % Immune cells

% These parameters are all the same as Dengue Virus.
params.mu      = 80;    % Target cell production: cell/(day ul)
params.alpha   = 1/3;   % Life span (1/day)
params.a       = 0.001; % Infectivity
params.beta    = 0.5;   % Infection period
params.v       = 0.001; % Immune cell efficiency
params.k       = 20;    % Virus multiplication
params.gamma   = 0.8;   % Virus clearence rate
params.c       = 0.01;  % Rate of immune cell production (due to density of infected cells)
params.d       = 0.03;  % Rate of immune cell production (due to interactions with infected cells)
params.delta   = 1/365; % 1/day
params.n       = 0.265; % Base immune cell production: cell/(day ul)

params.time_phase = 1;  % Simulate for 1 day.

[t, ~, I, ~, Z] = modelSystem(params);  % Solve set of ODEs.
t = t .* 24;            % Convert to hours

% First simulation
figure(1)
subplot(2,2,1)
semilogy(t,I)
axis([0 24 1 1e6])
xlabel('time (hours)')
ylabel('infected cells per uL blood')
subplot(2,2,3)
plot(t,Z)
axis([0 24 0 8e4])
xlabel('time (hours)')
ylabel('immune cells per uL blood')

type = 2; %#ok<NASGU>                   % Negate all immune cell activity
[t, ~, I, ~, Z] = modelSystem(params);  % Solve set of ODEs.
t = t .* 24;                            % Convert to hours

% Second simulation
figure(1)
subplot(2,2,2)
semilogy(t,I)
axis([0 24 1 1e6])
subplot(2,2,4)
plot(t,Z)
axis([0 24 0 8e4])

%% Figure 2
clear
close
clc
global type
global degredationTime
type = 3; %#ok<NASGU>

% As measured in Angel, S2E9.
params.S  = 1e6;      % Target cells
params.I  = 5;        % Infected cells
params.V  = 5;        % Viruses
params.Z  = 2000;     % Immune cells

% These parameters are all the same as Dengue Virus.
params.mu      = 80;    % Target cell production: cell/(day ul)
params.alpha   = 1/3;   % Life span (1/day)
params.a       = 0.001; % Infectivity
params.beta    = 0.5;   % Infection period
params.v       = 0.001; % Immune cell efficiency
params.k       = 20;    % Virus multiplication
params.gamma   = 0.8;   % Virus clearence rate
params.c       = 0.01;  % Rate of immune cell production (due to density of infected cells)
params.d       = 0.03;  % Rate of immune cell production (due to interactions with infected cells)
params.delta   = 1/365; % 1/day
params.n       = 0.265; % Base immune cell production: cell/(day ul)

params.time_phase = 1;  % Simulate for 1 day.
degradingTimes = [14.5 11 8.5];  % Hours until immunity is lost (representative of short, medium, long bite)

% Repeat simulation for multiple exposure rates 
volume = 1:10:500; % in mL

turningTime = zeros(length(volume),length(degradingTimes));
counter = 1;
for a = degradingTimes
    degredationTime = a; % Set immune system's degredation time
    counter2 = 1;
    
    for b = (volume  ./ 20) % 1 virus/uL blood == 20mL ingested volume
        
        params.I = b;   % Set relevant starting condition
        params.V = b;
        
        [t, ~, I, ~, Z] = modelSystem(params);  % Solve set of ODEs.
        t = t .* 24;                            % Convert to hours
    
        turningTime(counter2,counter) = max([0 t(find(I > 4000, 1 ))]); % Find when the 0.4% threshold is crossed
        counter2 = counter2 + 1;
    end
    counter = counter + 1;
end

turningTime(turningTime == 0) = NaN;

% Plot
plot(volume,turningTime(:,1),volume,turningTime(:,2),volume,turningTime(:,3))
legend({'short','medium','long'})
xlabel('volume of vampire blood ingested (mL)')
ylabel('time until fully turned (hrs)')

%% Supplementary Figure 1
clear
close
clc
global type
global degredationTime
type = 3; %#ok<NASGU>

% As measured in Angel, S2E9.
params.S  = 1e6;      % Target cells
params.I  = 5;        % Infected cells
params.V  = 5;        % Viruses
params.Z  = 2000;     % Immune cells

% These parameters are all the same as Dengue Virus.
params.mu      = 80;    % Target cell production: cell/(day ul)
params.alpha   = 1/3;   % Life span (1/day)
params.a       = 0.001; % Infectivity
params.beta    = 0.5;   % Infection period
params.v       = 0.001; % Immune cell efficiency
params.k       = 20;    % Virus multiplication
params.gamma   = 0.8;   % Virus clearence rate
params.c       = 0.01;  % Rate of immune cell production (due to density of infected cells)
params.d       = 0.03;  % Rate of immune cell production (due to interactions with infected cells)
params.delta   = 1/365; % 1/day
params.n       = 0.265; % Base immune cell production: cell/(day ul)

params.time_phase = 1;  % Simulate for 1 day.
degradingTimes = [14.5 13 11 8.5];  % Hours until immunity is lost (representative of short, medium, long bite)

% Repeat simulation for multiple exposure rates 
figure(1)
hold on
for a = degradingTimes
    
    degredationTime = a; % Set immune system's degredation time
        
    [t, ~, I, ~, Z] = modelSystem(params);  % Solve set of ODEs.
    t = t .* 24;                            % Convert to hours
    
    plot(t,I)
    line([0 24],[4e3 4e3],'Color','red','LineStyle','--')
end
hold off
set(gca, 'YScale', 'log')

%% Naive healing model for victim (calculation).
clear
close
clc
global type
global degredationTime
global healthResponse
type = 4; %#ok<NASGU>

% As measured in Angel, S2E9.
params.S  = 1e6;      % Target cells
params.I  = 5;        % Infected cells
params.V  = 5;        % Viruses
params.Z  = 2000;     % Immune cells

% These parameters are all the same as Dengue Virus.
params.mu      = 80;    % Target cell production: cell/(day ul)
params.alpha   = 1/3;   % Life span (1/day)
params.a       = 0.001; % Infectivity
params.beta    = 0.5;   % Infection period
params.v       = 0.001; % Immune cell efficiency
params.k       = 20;    % Virus multiplication
params.gamma   = 0.8;   % Virus clearence rate
params.c       = 0.01;  % Rate of immune cell production (due to density of infected cells)
params.d       = 0.03;  % Rate of immune cell production (due to interactions with infected cells)
params.delta   = 1/365; % 1/day
params.n       = 0.265; % Base immune cell production: cell/(day ul)

params.time_phase = 1;  % Simulate for 1 day.
degradingTimes = [14.5 13 11 8.5];  % Hours until immunity is lost (representative of short, medium, long bite)
healthResponse = 6;     % Time at which healing begins

% Repeat simulation for multiple exposure rates 
figure(1)
for a = degradingTimes
    
    degredationTime = a; % Set immune system's degredation time
        
    [t, ~, I, ~, Z] = modelSystem(params);  % Solve set of ODEs.
    t = t .* 24;                            % Convert to hours
    
    subplot(2,1,1)
    hold on
    plot(t,I)
    line([0 24],[4e3 4e3],'Color','red','LineStyle','--')
    
    subplot(2,1,2)
    hold on
    plot(t,Z)
end
subplot(2,1,1)
set(gca, 'YScale', 'log')
xlabel('time (hours)')
ylabel('infected cells per uL blood')
hold off

subplot(2,1,2)
xlabel('time (hours)')
ylabel('immune cells per uL blood')
hold off
