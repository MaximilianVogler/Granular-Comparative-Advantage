clear all;
close all;
clc;

rng(1);

tic
ALPHA = MCgenerate_vectorized(4);
toc