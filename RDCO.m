function [BestSol] = RDCO(fhd, nVar, VarMin, VarMax, nPop, MaxFEs, func_num)
% =========================================================================
% Core function for the Receptor-Driven Cycloidal Optimization (RDCO) algorithm
%
% As presented in the paper:
% "Receptor-Driven Cycloidal Optimization (RDCO): A Novel Bio-Inspired
%  Algorithm with Dynamic Geometric Search"
%
% Author: Lin Yang
% Affiliation: Mudanjiang Normal University
% Email: 15765057075@163.com
%
% Please cite the main paper if you use this code in your research.
% (Full citation details will be provided upon the paper's publication.)
% =========================================================================
%
% INPUTS
%   fhd         : Function handle for the objective function
%   nVar        : Problem dimensionality
%   VarMin/VarMax: Search space boundaries
%   nPop        : Population size
%   MaxFEs      : Maximum function evaluations (stopping criterion)
%   func_num    : CEC function number (if applicable)
%
% OUTPUTS
%   BestSol     : A struct containing the best solution (Position and Cost)
% =========================================================================

%% 0. Initialization
% --- Algorithm-specific constants ---
ALPHA_COEFF = 0.3; % Weighting factor for receptor-driven vs. random behavior
PHI_CONST = (sqrt(5) - 1) / 2; % Golden ratio, for epicycloid operator
EPSILON = 1e-8;     % To prevent division by zero

% --- Population and Global Best Initialization ---
empty_individual.Position = [];
empty_individual.Cost = [];
pop = repmat(empty_individual, nPop, 1);
BestSol.Cost = inf;

current_FEs = 0;
for i = 1:nPop
    pop(i).Position = unifrnd(VarMin, VarMax, [1, nVar]);
    pop(i).Cost = feval(fhd, pop(i).Position', func_num);
    current_FEs = current_FEs + 1;
    if pop(i).Cost < BestSol.Cost
        BestSol = pop(i);
    end
    if current_FEs >= MaxFEs; return; end
end

%% 1. Main Optimization Loop
iter = 1;
while current_FEs < MaxFEs

    % --- 1.1. Fitness Sorting and Global Best Update ---
    Costs = [pop.Cost];
    [~, SortIdx] = sort(Costs);
    if pop(SortIdx(1)).Cost < BestSol.Cost
        BestSol = pop(SortIdx(1));
    end
    gbest_pos = BestSol.Position;

    % --- 1.2. Population Tiering for Role Assignment ---
    N1 = floor(nPop/3);
    N2 = floor(nPop/3);
    S = min([N1, N2, (nPop - N1 - N2)]); % Number of collaborative subgroups

    if S > 0
        top_indices = SortIdx(1:S);
        mid_indices = SortIdx(N1+1 : N1+S);
        bot_indices = SortIdx(N1+N2+1 : N1+N2+S);
    end

    % --- 1.3. Receptor-Driven Control Mechanism (Hill Equation) ---
    % This section implements the adaptive control based on agent performance (rank).
    current_ranks = zeros(1, nPop);
    current_ranks(SortIdx) = 1:nPop;

    if S > 0
        % Normalized rank (L), serves as the "stimulus" for the Hill equation.
        L_mid = (nPop - current_ranks(mid_indices)) / (nPop - 1 + EPSILON);
        L_bot = (nPop - current_ranks(bot_indices)) / (nPop - 1 + EPSILON);
    end

    % Dynamic Hill model parameters adapt based on search progress (tau).
    tau = current_FEs / MaxFEs;
    n_mu = 2 - 0.8*tau + 0.1*randn();      % Hill coefficient for exploration (mu)
    n_kappa = 0.3 + 0.6*tau + 0.05*randn(); % Hill coefficient for exploitation (kappa)
    K_mu = 0.8 - 0.4*tau;                 % Dissociation constant for exploration
    K_kappa = 0.2 + 0.4*tau;              % Dissociation constant for exploitation

    % Calculate activation levels (E) based on the Hill equation.
    if S > 0
        E_mu_mid = (L_mid.^n_mu) ./ (K_mu^n_mu + L_mid.^n_mu + EPSILON);
        E_kappa_mid = (L_mid.^n_kappa) ./ (K_kappa^n_kappa + L_mid.^n_kappa + EPSILON);
        E_overall_mid = max(E_mu_mid, E_kappa_mid);

        E_mu_bot = (L_bot.^n_mu) ./ (K_mu^n_mu + L_bot.^n_mu + EPSILON);
        E_kappa_bot = (L_bot.^n_kappa) ./ (K_kappa^n_kappa + L_bot.^n_kappa + EPSILON);
        E_overall_bot = max(E_mu_bot, E_kappa_bot);
    end

    if S > 0
        x_subbest_indices = zeros(1,S);
        x_subbest_E_mu = zeros(1,S);
        x_subbest_E_kappa = zeros(1,S);
    end

    %% 1.4. Geometric Search Operators (Cycloids)
    for k = 1:S
        if current_FEs >= MaxFEs; break; end

        x_top = pop(top_indices(k)).Position;
        x_mid = pop(mid_indices(k)).Position;
        x_bot = pop(bot_indices(k)).Position;

        % --- Hypocycloid Operator for Local Exploitation (Middle-tier agent) ---
        R1 = x_top - x_mid;
        norm_R1 = norm(R1);
        d1 = R1 / (norm_R1 + EPSILON);
        u1 = randn(1, nVar) - dot(randn(1, nVar), d1) * d1; u1 = u1 / (norm(u1) + EPSILON);
        theta1 = 2*pi * (ALPHA_COEFF * E_overall_mid(k) + (1-ALPHA_COEFF)*rand());
        v1 = cos(theta1)*d1 + sin(theta1)*u1;
        k1 = 0.1 + 0.4*tau; % Hypocycloid parameter
        theta1_prime = (1/k1 - 1) * theta1;
        v1_prime = cos(theta1_prime)*d1 + sin(theta1_prime)*u1;
        new_mid_pos = x_top + (1-k1)*norm_R1*v1 + k1*norm_R1*v1_prime;

        % --- Epicycloid Operator for Global Exploration (Bottom-tier agent) ---
        R2 = x_top - x_bot;
        norm_R2 = norm(R2);
        d2 = R2 / (norm_R2 + EPSILON);
        u2 = randn(1, nVar) - dot(randn(1, nVar), d2) * d2; u2 = u2 / (norm(u2) + EPSILON);
        theta2 = 2*pi * (ALPHA_COEFF * E_overall_bot(k) + (1-ALPHA_COEFF)*rand());
        v2 = cos(theta2)*d2 + sin(theta2)*u2;
        k2 = PHI_CONST; % Epicycloid parameter
        theta2_prime = (1 + 1/k2) * theta2;
        v2_prime = cos(theta2_prime)*d2 + sin(theta2_prime)*u2;
        new_bot_pos = x_top + (1+k2)*norm_R2*v2 - k2*norm_R2*v2_prime;

        % --- Fusion Update for Top-tier agent ---
        new_top_pos = x_top + (new_mid_pos - x_top)*E_overall_mid(k) + (new_bot_pos - x_top)*E_overall_bot(k);

        % --- Evaluation and Greedy Selection ---
        new_mid_pos = max(min(new_mid_pos, VarMax), VarMin);
        new_mid_cost = feval(fhd, new_mid_pos', func_num); current_FEs=current_FEs+1; if current_FEs>=MaxFEs; break; end
        if new_mid_cost < pop(mid_indices(k)).Cost
            pop(mid_indices(k)).Position = new_mid_pos; pop(mid_indices(k)).Cost = new_mid_cost;
            if new_mid_cost < BestSol.Cost; BestSol = pop(mid_indices(k)); end
        end

        new_bot_pos = max(min(new_bot_pos, VarMax), VarMin);
        new_bot_cost = feval(fhd, new_bot_pos', func_num); current_FEs=current_FEs+1; if current_FEs>=MaxFEs; break; end
        if new_bot_cost < pop(bot_indices(k)).Cost
            pop(bot_indices(k)).Position = new_bot_pos; pop(bot_indices(k)).Cost = new_bot_cost;
            if new_bot_cost < BestSol.Cost; BestSol = pop(bot_indices(k)); end
        end

        new_top_pos = max(min(new_top_pos, VarMax), VarMin);
        new_top_cost = feval(fhd, new_top_pos', func_num); current_FEs=current_FEs+1; if current_FEs>=MaxFEs; break; end
        if new_top_cost < pop(top_indices(k)).Cost
            pop(top_indices(k)).Position = new_top_pos; pop(top_indices(k)).Cost = new_top_cost;
            if new_top_cost < BestSol.Cost; BestSol = pop(top_indices(k)); end
        end

        sub_indices = [top_indices(k), mid_indices(k), bot_indices(k)];
        sub_costs = [pop(top_indices(k)).Cost, pop(mid_indices(k)).Cost, pop(bot_indices(k)).Cost];
        [~, min_idx] = min(sub_costs);
        x_subbest_indices(k) = sub_indices(min_idx);

        rank_sub = current_ranks(x_subbest_indices(k)); L_sub = (nPop-rank_sub)/(nPop-1+EPSILON);
        x_subbest_E_mu(k) = (L_sub^n_mu) / (K_mu^n_mu + L_sub^n_mu + EPSILON);
        x_subbest_E_kappa(k) = (L_sub^n_kappa) / (K_kappa^n_kappa + L_sub^n_kappa + EPSILON);
    end
    if current_FEs >= MaxFEs; break; end

    %% 1.5. Update for Remaining (Non-subgroup) Individuals
    is_in_subgroup = false(1,nPop);
    if S>0; is_in_subgroup(top_indices)=true; is_in_subgroup(mid_indices)=true; is_in_subgroup(bot_indices)=true; end
    remaining_indices = find(~is_in_subgroup);
    for i=1:length(remaining_indices)
        if current_FEs >= MaxFEs; break; end
        idx = remaining_indices(i);
        pop(idx).Position = unifrnd(VarMin, VarMax, [1,nVar]);
        pop(idx).Cost = feval(fhd, pop(idx).Position', func_num); current_FEs=current_FEs+1;
        if pop(idx).Cost < BestSol.Cost; BestSol = pop(idx); end
    end
    if current_FEs >= MaxFEs; break; end

    %% 1.6. Global Guidance for Subgroup Leaders
    for k = 1:S
        if current_FEs >= MaxFEs; break; end
        idx_subbest = x_subbest_indices(k);
        pos_subbest = pop(idx_subbest).Position;
        E_mu = x_subbest_E_mu(k); E_kappa = x_subbest_E_kappa(k);
        R3 = gbest_pos - pos_subbest; norm_R3 = norm(R3); if norm_R3<EPSILON; continue; end
        d3 = R3 / norm_R3;
        u3 = randn(1,nVar) - dot(randn(1,nVar),d3)*d3; u3=u3/(norm(u3)+EPSILON);
        theta3 = 2*pi*(ALPHA_COEFF*max(E_mu,E_kappa) + (1-ALPHA_COEFF)*rand());
        v3 = cos(theta3)*d3 + sin(theta3)*u3;
        if E_mu > E_kappa % Hypocycloid towards global best
            theta3_prime = (1/k1-1)*theta3;
            v3_prime = cos(theta3_prime)*d3 + sin(theta3_prime)*u3;
            new_pos = gbest_pos + (1-k1)*norm_R3*v3 + k1*norm_R3*v3_prime;
        else % Epicycloid towards global best
            theta3_prime = (1+1/k2)*theta3;
            v3_prime = cos(theta3_prime)*d3 + sin(theta3_prime)*u3;
            new_pos = gbest_pos + (1+k2)*norm_R3*v3 - k2*norm_R3*v3_prime;
        end
        new_pos = max(min(new_pos, VarMax), VarMin);
        new_cost = feval(fhd, new_pos', func_num); current_FEs=current_FEs+1; if current_FEs>=MaxFEs; break; end
        if new_cost < pop(idx_subbest).Cost
            pop(idx_subbest).Position = new_pos; pop(idx_subbest).Cost = new_cost;
            if new_cost < BestSol.Cost; BestSol = pop(idx_subbest); end
        end
    end
    if current_FEs >= MaxFEs; break; end

    iter = iter + 1;
end
end