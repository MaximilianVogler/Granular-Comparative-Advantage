%% This finds the 20 best points from the inputted data file as in steps 4-5 of the estimation procedure.

%% Housekeeping
clear all;
close all;
clc;

addpath('Data');

%% Load moments from both data and estimation

datamoments = csvread('Data_Moments.csv');

load('estimation_seed1_grid5');

%% Set up weighting matrix W

ref_weight = datamoments;
ref_weight(14:15) = datamoments(12:13);
std_indices = [1,3,5,7,9];
ref_weight(std_indices) = datamoments(std_indices)/3;
W = diag((1./(ref_weight.^2)));

%% Compute loss function for each parameter combination

NumGrid = size(Momarray,1);
Loss = zeros(NumGrid,1);

for i = 1:NumGrid
    Loss(i) = (Momarray(i,:)'-datamoments)'*W*(Momarray(i,:)'-datamoments);
end

%% Find 20 best points (i.e. with lowest loss)

cutoff = 0.001;
indicesP = find(Loss<quantile(Loss,cutoff));

num = size(indicesP,1);

if num ~= 20
    error('We are not looking at the 20 best points')
end

paramsP = zeros(num,6);

for i = 1:num
    paramsP(i,1:5) = Paramarray(indicesP(i),1:5);
    paramsP(i,6) = Loss(indicesP(i));
end

disp(paramsP)

save('GridOptimization_seed1_grid5','paramsP')



