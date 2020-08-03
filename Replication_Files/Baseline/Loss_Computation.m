%% This finds the 20 best points from the inputted data file as in steps 4-5 of the estimation procedure.

%% Housekeeping
clear all;
close all;
clc;

addpath('Data');
addpath('Results');

%% Load moments from both data and estimation

datamoments = csvread('Data_Moments.csv');

load('Moments_seed1_grid6','StandardMom');

Momarray = StandardMom';

%% Set up weighting matrix W

ref_weight = datamoments;
ref_weight(14:15) = datamoments(12:13);
std_indices = [1,3,5,7,9];
ref_weight(std_indices) = datamoments(std_indices)/3;
W = diag((1./(ref_weight.^2)));

%% Compute Overall Loss

Loss = (Momarray-datamoments)'*W*(Momarray-datamoments);

%% Compute individual losses

losses = zeros(15,1);

for i = 1:15
   losses(i) = (Momarray(i)'-datamoments(i))' * W(i,i) * (Momarray(i)'-datamoments(i));
   losses(i) = losses(i) / Loss *100;
end

Loss_Save = sqrt(Loss) / 3 * 100;


save('Results/Losses','losses','Loss_Save')
