% =========================================================================
% Main Benchmark Script for RDCO
%
% This script runs the benchmark experiments for the paper:
% "Receptor-Driven Cycloidal Optimization (RDCO): A Novel Bio-Inspired 
%  Algorithm with Dynamic Geometric Search"
% =========================================================================

clc;
clear;
close all;

%% 1. Experiment Configuration
functions_to_test = 1:12; % CEC 2022 functions to test
num_runs = 30;            % Number of independent runs for statistics
nVar = 20;                % Problem dimensionality

%% 2. Initialization
num_functions = length(functions_to_test);
FinalBestCost = zeros(num_functions, num_runs);
ExecutionTimes = zeros(num_functions, num_runs);

%% 3. Main Experiment Loop
for func_idx = 1:num_functions
    func_num = functions_to_test(func_idx);
    
    fprintf('\n--- Evaluating RDCO on F%d (Dim=%d) ---\n', func_num, nVar);

    % Problem Definition
    fhd = str2func('cec22_test_func');
    VarMin = -100;
    VarMax = 100;

    % Algorithm Parameters
    MaxFEs = 45000;
    nPop = 30;

    % Perform Runs
    for r = 1:num_runs
        tic;
        BestSol = RDCO(fhd, nVar, VarMin, VarMax, nPop, MaxFEs, func_num);
        ExecutionTimes(func_idx, r) = toc;
        FinalBestCost(func_idx, r) = BestSol.Cost;
        fprintf('Run %2d/%d | Best Fitness: %.4e\n', r, num_runs, BestSol.Cost);
    end
end

%% 4. Display Final Results
fprintf('\n================== OVERALL RESULTS ==================\n');
results_table = table();
results_table.Function = functions_to_test';
results_table.Mean_Fitness = mean(FinalBestCost, 2);
results_table.Std_Dev = std(FinalBestCost, 0, 2);
results_table.Avg_Time_s = mean(ExecutionTimes, 2);

disp(results_table);