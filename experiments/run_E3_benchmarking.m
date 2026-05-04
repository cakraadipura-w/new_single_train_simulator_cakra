%% run_E3_benchmarking.m — Experiment E3: Strong Benchmarking (G3)
%
% Segment : IS04 (1642 m)
% Methods :
%   M1 — DP            (deterministic, loose slack only)
%   M2 — NSGA-II vanilla + original CC_CR
%   M3 — MOPSO         + original CC_CR
%   M4 — Improved DE   + original CC_CR
%   M5 — DQN           + original CC_CR
%   M6 — Proposed      (improved CC_CR + BRL-SDE-NSGA-II)
%
% Slack conditions:
%   Tight : 8.4%  → T_target = T_min × 1.084
%   Loose : 20%   → T_target = T_min × 1.20
%
% Outputs : E3_results.csv, E3_summary.txt, Pareto front plot

clc;
addpath(genpath(fullfile(fileparts(mfilename('fullpath')), '..')));

%% ===== CONFIG =====
RS_FILE     = 'rollingstock_Guangzhou_L7.m';
SLACK_TIGHT = 8.4;
SLACK_LOOSE = 20.0;
POP_SIZE    = 200;
ITERATIONS  = 300;

% Jika dipanggil dari main_experiments.m, gunakan ACTIVE_SEG dari workspace
if exist('ACTIVE_SEG','var') && ~isempty(ACTIVE_SEG)
    ROUTE_IS04 = ACTIVE_SEG.file;
    N_RUNS     = ACTIVE_NRUNS;
    OUT_DIR    = fullfile(fileparts(mfilename('fullpath')), '..', 'experiment_results', ACTIVE_SEG.name);
else
    ROUTE_IS04 = 'Guangzhou_Line7_IS04_5.200-6.842km.mat';
    N_RUNS     = 30;
    OUT_DIR    = fullfile(fileparts(mfilename('fullpath')), '..', 'experiment_results');
end
if ~exist(OUT_DIR,'dir'), mkdir(OUT_DIR); end

METHODS = struct( ...
    'id',          {'M1',    'M2',        'M3',    'M4',       'M5',  'M6'}, ...
    'name',        {'DP','NSGA-II vanilla','MOPSO','Improved DE','DQN','Proposed'}, ...
    'solver',      {'dp',   'nsga2',     'mopso', 'de_moea',  'dqn', 'nsga2'}, ...
    'use_improved',{false,   false,       false,   false,      false,  true}, ...
    'nsga2_var',   {'orig', 'original',  'none',  'none',     'none','rl_sde'}, ...
    'determ',      {true,    false,       false,   false,      false,  false}, ...
    'skip_tight',  {true,    false,       false,   false,      false,  false} );

SLACK_CONDS = [SLACK_TIGHT, SLACK_LOOSE];

%% ===== GLOBALS + SETUP =====
global use_improved nsga2_variant pop_size iterations dimension
global vel_profile time_obj_max parallel_use show_progress driving_strategy

show_progress    = false;
parallel_use     = true;
driving_strategy = "CC_CR";

cfg = struct('route_file', ROUTE_IS04, 'rollingstock_file', RS_FILE, ...
    'driving_strategy',"CC_CR", 'pop_size',POP_SIZE, 'iterations',ITERATIONS, ...
    'use_improved',true, 'nsga2_variant','rl_sde', 'parallel_use',true, 'sim',struct());
info = setup_project(cfg);
setup_parallel_pool(cfg, info);

[T_min, E_base] = flat_out_baseline_rk4(info.routePath, info.rollingstockPath);
fprintf('T_min=%.2f s | E_baseline=%.4f kWh\n', T_min, E_base);

T_TARGETS = struct('tight', T_min*(1+SLACK_TIGHT/100), 'loose', T_min*(1+SLACK_LOOSE/100));
fprintf('T_tight=%.2f s (%.1f%%) | T_loose=%.2f s (%.1f%%)\n', ...
    T_TARGETS.tight, SLACK_TIGHT, T_TARGETS.loose, SLACK_LOOSE);

%% ===== MAIN LOOP =====
all_rows = {};
raw_store = struct('method',{},'slack',{},'run',{},'F_nd',{},'F_all',{});

for sci = 1:2
    slack  = SLACK_CONDS(sci);
    T_tgt  = T_min * (1 + slack/100);
    time_obj_max = T_tgt * 1.50;   % lebih longgar — mode santai tetap bisa ditemukan
    slack_label  = ternary(slack==SLACK_TIGHT,'tight','loose');

    fprintf('\n===== Slack=%.4g%% (%s) | T_target=%.2f s =====\n', slack, slack_label, T_tgt);
    all_F_cond = zeros(0,2);

    for mi = 1:numel(METHODS)
        M = METHODS(mi);
        if slack==SLACK_TIGHT && M.skip_tight
            fprintf('  [%s] %s — skipped (tight slack)\n', M.id, M.name);
            continue;
        end

        n_runs_m = ternary(M.determ, 1, N_RUNS);
        fprintf('\n  [%s] %s | %d run(s)\n', M.id, M.name, n_runs_m);

        use_improved  = M.use_improved;
        nsga2_variant = M.nsga2_var;
        init_worker_globals(info.routePath, info.rollingstockPath, cfg.sim, use_improved);
        pool = gcp('nocreate');
        if ~isempty(pool)
            wait(parfevalOnAll(@init_worker_globals, 0, info.routePath, ...
                info.rollingstockPath, cfg.sim, use_improved));
        end
        pop_size   = POP_SIZE;
        iterations = ITERATIONS;

        % Push pop_size/iterations/time_obj_max/nsga2_variant to workers
        if ~isempty(pool)
            wait(parfevalOnAll(@push_worker_globals, 0, POP_SIZE, ITERATIONS, T_tgt*1.50, M.nsga2_var));
        end

        pops  = cell(n_runs_m,1);
        tims  = zeros(n_runs_m,1);

        parfor (run_id = 1:n_runs_m, 4)
            if M.determ, rng(42,'twister'); else, rng(run_id,'twister'); end
            t0_ = tic;
            pops{run_id}  = run_solver(M.solver, vel_profile);
            tims(run_id)  = toc(t0_);
        end

        for run_id = 1:n_runs_m
            pop_r = pops{run_id};
            if isempty(pop_r) || size(pop_r,2)<dimension+2
                all_rows{end+1}={M.id,M.name,slack,run_id,NaN,NaN,NaN,tims(run_id),NaN}; %#ok<AGROW>
                continue;
            end
            F_all = double(pop_r(:,dimension+1:dimension+2));
            valid = F_all(:,1)<1e5;
            F_nd  = nd_filter(F_all(valid,:));
            if ~isempty(F_nd), all_F_cond=[all_F_cond;F_nd]; end %#ok<AGROW>
            raw_store(end+1) = struct('method',M.id,'slack',slack,'run',run_id, ...
                'F_nd',F_nd,'F_all',F_all); %#ok<AGROW>
            all_rows{end+1}={M.id,M.name,slack,run_id,NaN,NaN,NaN,tims(run_id),NaN}; %#ok<AGROW>
            fprintf('    run %02d | F1=%d | %.1f s\n', run_id, size(F_nd,1), tims(run_id));
        end
    end

    % --- compute metrics vs global ref ---
    F_ref = nd_filter(all_F_cond);
    fprintf('  Reference front: %d pts\n', size(F_ref,1));
    idx_cond = find([raw_store.slack]==slack);
    ri_in_cond = 0;
    for ri_all = 1:numel(all_rows)
        r = all_rows{ri_all};
        if r{3}~=slack, continue; end
        ri_in_cond = ri_in_cond+1;
        if ri_in_cond > numel(idx_cond), break; end
        rs = raw_store(idx_cond(ri_in_cond));
        [hv_v,igd_v,sp_v,fr_v] = compute_metrics(rs.F_all, F_ref, T_tgt);
        all_rows{ri_all}{5}=hv_v; all_rows{ri_all}{6}=igd_v;
        all_rows{ri_all}{7}=sp_v; all_rows{ri_all}{9}=fr_v;
    end
end

%% ===== SAVE CSV =====
T_csv = cell2table(vertcat(all_rows{:}), ...
    'VariableNames',{'method_id','method_name','slack_pct','run_id', ...
                     'hv','igd','spread','time_s','feasible'});
writetable(T_csv, fullfile(OUT_DIR,'E3_results.csv'));

%% ===== SUMMARY + WILCOXON =====
fid = fopen(fullfile(OUT_DIR,'E3_summary.txt'),'w');
fprintf(fid,'Experiment E3 — Strong Benchmarking\n');
fprintf(fid,'Segment: IS04 | T_min=%.2f s\n\n', T_min);

for sci=1:2
    slack=SLACK_CONDS(sci);
    T_tgt=T_min*(1+slack/100);
    label=ternary(slack==SLACK_TIGHT,'TIGHT','LOOSE');
    fprintf(fid,'\n=== Slack=%.4g%% (%s) T_target=%.2f s ===\n',slack,label,T_tgt);
    fprintf(fid,'%-18s  %8s  %8s  %10s  %8s  %8s\n', ...
        'Method','HV_mean','HV_std','Time_mean','Feas','EnSav%%');

    hv_M6 = [];
    for mi=1:numel(METHODS)
        M=METHODS(mi);
        if slack==SLACK_TIGHT&&M.skip_tight, continue; end
        mask=strcmp(T_csv.method_id,M.id)&T_csv.slack_pct==slack;
        hv_v=T_csv.hv(mask); t_v=T_csv.time_s(mask); fr_v=T_csv.feasible(mask);

        E_best_arr=zeros(sum(mask),1);
        idx_rs=find(strcmp({raw_store.method},M.id)&[raw_store.slack]==slack);
        for ii=1:min(numel(idx_rs),numel(E_best_arr))
            F_nd=raw_store(idx_rs(ii)).F_nd;
            if ~isempty(F_nd)
                fnd=F_nd(F_nd(:,1)<=T_tgt,:);
                if ~isempty(fnd), E_best_arr(ii)=min(fnd(:,2)); end
            end
        end
        esav=mean(energy_saving_percent(E_best_arr(E_best_arr>0), E_base),'omitnan');

        line=sprintf('%-18s  %8.4f  %8.4f  %10.1f  %8.3f  %8.2f', ...
            [M.id ': ' M.name], mean(hv_v,'omitnan'),std(hv_v,'omitnan'), ...
            mean(t_v,'omitnan'),mean(fr_v,'omitnan'),esav);
        disp(line); fprintf(fid,'%s\n',line);
        if strcmp(M.id,'M6'), hv_M6=hv_v; end
    end

    if ~isempty(hv_M6)
        fprintf(fid,'\nWilcoxon signrank HV: M6 vs others — slack=%.4g%%\n',slack);
        for mi=1:numel(METHODS)
            M=METHODS(mi);
            if any(strcmp(M.id,{'M1','M6'})) || (slack==SLACK_TIGHT && M.skip_tight), continue; end
            mask=strcmp(T_csv.method_id,M.id)&T_csv.slack_pct==slack;
            hv_m=T_csv.hv(mask);
            [x,y]=match_paired(hv_M6,hv_m);
            if numel(x)>=5 && exist('signrank','file')==2
                p_val=signrank(x,y);
            else
                p_val=wilcoxon_srtest(x,y);
            end
            fprintf(fid,'  M6 vs %s (%s): p=%.4f%s\n',M.id,M.name,p_val, ...
                ternary(p_val<0.05,'  *',''));
        end
    end
end
fclose(fid);

%% ===== PLOTS: PARETO FRONTS FOR TIGHT SLACK =====
colors_m = lines(numel(METHODS));
figure(3001); clf; hold on;
for mi=1:numel(METHODS)
    M=METHODS(mi);
    if M.skip_tight, continue; end
    idx=strcmp({raw_store.method},M.id)&[raw_store.slack]==SLACK_CONDS(1);
    rs=raw_store(idx);
    F_u=zeros(0,2);
    for ri=1:numel(rs)
        if ~isempty(rs(ri).F_nd), F_u=[F_u;rs(ri).F_nd]; end %#ok<AGROW>
    end
    F_u=nd_filter(F_u);
    if ~isempty(F_u)
        [~,ord]=sort(F_u(:,1));
        lw=ternary(strcmp(M.id,'M6'),2.5,1.2);
        plot(F_u(ord,1),F_u(ord,2),'-o','Color',colors_m(mi,:),'LineWidth',lw, ...
            'DisplayName',[M.id ': ' M.name]);
    end
end
T_tight=T_min*(1+SLACK_TIGHT/100);
xline(T_tight,'--k','T_{target}','LabelVerticalAlignment','bottom');
xlabel('Running time (s)'); ylabel('Energy (kWh)');
title(sprintf('E3: Pareto Front — IS04 Tight Slack (%.4g%%)',SLACK_TIGHT));
legend('Location','northeast'); grid on;
saveas(gcf, fullfile(OUT_DIR,'E3_pareto_tight.png'));

fprintf('\nE3 complete → %s\n', OUT_DIR);

%% ===== HELPERS =====
function pop = run_solver(solver, vp)
    global vel_profile %#ok<NASGU>
    try
        switch lower(solver)
            case 'nsga2',   [pop,~,~]=nsga2_main(vp);
            case 'mopso',   pop=mopso_main(vp);
            case 'de_moea', [pop,~,~]=de_moea_main(vp);
            case 'dp',      [pop,~,~]=dp_main(vp);
            case 'dqn',     [pop,~,~]=dqn_main(vp);
            otherwise,      error('Unknown solver: %s',solver);
        end
    catch ME
        warning('[run_solver] %s failed: %s', solver, ME.message);
        pop = zeros(0,4);
    end
end

function F_nd = nd_filter(F)
    if isempty(F), F_nd=zeros(0,2); return; end
    n=size(F,1); k=true(n,1);
    for i=1:n
        if ~k(i),continue;end
        for j=1:n
            if i==j||~k(j),continue;end
            if all(F(j,:)<=F(i,:))&&any(F(j,:)<F(i,:)),k(i)=false;break;end
        end
    end
    F_nd=F(k,:);[~,o]=sort(F_nd(:,1));F_nd=F_nd(o,:);
end

function [x,y]=match_paired(a,b)
    n=min(numel(a),numel(b));
    x=a(1:n); y=b(1:n);
    ok=~isnan(x)&~isnan(y); x=x(ok); y=y(ok);
end

function v=ternary(c,a,b)
    if c,v=a;else,v=b;end
end
