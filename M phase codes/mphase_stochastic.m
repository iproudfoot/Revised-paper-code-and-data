clear all
close all
%% Initial conditions
ICN = 20;
T = 24; % max time
delta = 0.01; % time step
%% Defining initial populations in cell cycle compartments
r = rem(0.05*ICN,3);
G = 0.95*ICN;
if r == 0
    M1 = 0.05*ICN/3
    M2 = 0.05*ICN/3
    M3 = 0.05*ICN/3
else
M1 = floor(0.05*ICN/3) + 1
if r == 1
M2 = floor(0.05*ICN/3)
end
if r == 2
M2 = floor(0.05*ICN/3) + 1
end
M3 = floor(0.05*ICN/3)
end
%% Defining the transition rates
k1 = 0.04; % transition rate into m phase (per hr)
lambda1 = 4.5; % transition rate between each mphase compartment.
%%
t = 0; % set time to 0

while t < T % loop over times less than end time
    
%% Reset number of transitions
no_transitions_G_M1 = 0;
no_transitions_M1_M2 = 0;
no_transitions_M2_M3 = 0;
no_transitions_M3_G = 0;
tot_trans_G_M1 = 0;
%% G2 to mphase transition
if G > 0 
event_log = zeros(1,G);
R = rand([1, G]);
event_log(R<exp(-delta/20))=1;
no_transitions_G_M1 = sum(event_log);
end

%% M1 to M2 transtition
if M1 > 0
event_log = zeros(1,M1);
R = rand([1,M1]);
event_log(R<exp(-lambda1*delta))=1;
no_transitions_M1_M2 = sum(event_log);
end

%% M2 to M3 transition
if M2 > 0 
event_log = zeros(1,M2);
R = rand([1,M2]);
event_log(R<exp(-lambda1*delta))=1;
no_transitions_M2_M3 = sum(event_log);
end

%% M3 to G1 transition
if M3 > 0 
event_log = zeros(1,M3);
R = rand([1,M3]);
event_log(R<exp(-lambda1*delta))=1;
no_transitions_M3_G = sum(event_log);
end

%% Update populations
G = G - no_transitions_G_M1 + 2*no_transitions_M3_G;
M1 = M1 - no_transitions_M1_M2 + no_transitions_G_M1;
M2 = M2 - no_transitions_M2_M3 + no_transitions_M1_M2;
M3 = M3 - no_transitions_M3_G + no_transitions_M2_M3;
tot_trans_G_M1 = tot_trans_G_M1 + no_transitions_G_M1;
%% Move to next time step
t = t + delta;
end
