function [M, v, zer] = load_pija()
M = dlmread('m.csv');
v = dlmread('v.csv');
zer = dlmread('zeros.csv');