%% Parameter
clc;
clear;
close all;
num_training = 300;
%% Data
load airline_data
num_class = num_class;
num_flight = num_flight;
capacity = capacity;
price = price;
%% Training data
arrival_rate = arrival_rate_2011;
v = v_2011;
w = w_2011;
v0 = v0_2011;
%% SBLP with spike
SBLP_with_spike;
save data_post_sblp;
%% Testing
training;
%% Training data
% arrival_rate = arrival_rate_2012;
% v = v_2012;
% w = w_2012;
% v0 = v0_2012;
%% send to condor
