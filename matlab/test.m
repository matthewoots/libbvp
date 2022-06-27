clc
clear all
close all

s0(1:3,1) = [5; -5; 10];
s0(4:6,1) = [1; -1; -2];
s0(7:9,1) = [0; 0; 0];
 
d0(1:3,1) = [0; -0; 0.5];
d0(4:6,1) = [0; -0; 0];
d0(7:9,1) = [0; -0; 0];

timeint = 0.1;
runtime = 5;

%% Original BVP without finding landing optimal curve
[intervals_original, sOpt_original] = ...
    boundary_value_problem(timeint, runtime, d0, s0);