% Workspace Preparation
% This script is meant to be run at the beginning of each script in this
% project to prepare MATLAB with paths and other code that is redundant in
% each script
%
% Matt Kmiecik
%
% Started 28 Jan 2020
%

% Sets working directory ----
'C:\Users\pains\Desktop\matt-eeglab'; % modify this path for your machine

% Starts EEGLAB ----
[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;

% Loads in participant information ----
% modify this path for your machine
[NUM,TXT,RAW] = xlsread('C:\Users\pains\Desktop\matt-eeglab\data\0-global-vars\vis-subj-info.xlsx');
