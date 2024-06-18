%% refresh
clear;
close all;
clc;
rng(0);

%%%%%%%%%%%%%%%%%%%%%%
% 第一组
% 32 0.0164    0.1315 阈值-1
% 第二组
% 29 0.0136    0.5341 阈值-1
% 第四组
% 30 0.0162    0.1499 阈值-1
% 第六组
% 30 0.0069    0.1775 阈值-1
% 第七组
% 31 0.0150    0.1797 阈值-1
%%%%%%%%%%%%%%%%%%%%
%仿真预设

% 30 0.3192    0.9433
% 50 0.1784    0.4000
% 70 0.1670    0.3453
% all 0.0375    0.1504
%26 0.0249    0.2604
%% add path for including some tool functions
addpath('func');

%% params

% fig4b.graph_file = './data/final_data/2_evenness_rhomboid_fig4b_test.mat';

fig4b.graph_file = './data/final_data/final_real_experiment6.mat';
%实物实验12367可用
disp('------------------------------------------------------------------');
disp('Plot Fig.4(b) in Fig. 2.');
fig4b.eps = 1e-2;
fig4b.fig.title = 'Fig.4(b)';
fig4b.fig.view_a = 30; fig4b.fig.view_e = 15;
calib_func1(fig4b);



