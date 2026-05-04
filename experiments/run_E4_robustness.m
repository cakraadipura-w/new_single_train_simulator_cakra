%% run_E4_robustness.m — Experiment E4: Robustness Evaluation (G4)
%
% Segment  : IS04 (1642 m)
% Nominal  : AW2 load (mass_factor=1.27), η=0.85, normal speed limit, T_target=180 s
% Step 1   : Run Config D 30 times → 30 Pareto fronts
% Step 2   : Select best solution per run (min E feasible or min T)
% Step 3   : Evaluate each candidate under 18 uncertainty scenarios
%            (3 loads × 3 efficiencies × 2 speed limits = 18)
% Step 4   : Compute F_rate and energy degradation; robust = F_rate ≥ 0.9
% Outputs  : E4_results.csv, E4_summary.txt, boxplot, bar chart

clc;
addpath(genpath(fullfile(fileparts(mfilename('fullpath')), '..')));

%% ===== CONFIG =====
RS_FILE    = 'rollingstock_Guangzhou_L7.m';
FEAS_TOL   = 1.01;
ROBUST_THR = 0.9;
POP_SIZE   = 200;
ITERATIONS = 300;

% Jika dipanggil dari main_experiments.m, gunakan ACTIVE_SEG dari workspace
if exist('ACTIVE_SEG','var') && ~isempty(ACTIVE_SEG)
    ROUTE_IS04 = ACTIVE_SEG.file;
    T_TARGET   = ACTIVE_SEG.T_sched;
    N_RUNS     = ACTIVE_NRUNS;
    OUT_DIR    = fullfile(fileparts(mfilename('fullpath')), '..', 'experiment_results', ACTIVE_SEG.name);
else
    ROUTE_IS04 = 'Guangzhou_Line7_IS04_5.200-6.842km.mat';
    T_TARGET   = 180;
    N_RUNS     = 30;
    OUT_DIR    = fullfile(fileparts(mfilename('fullpath')), '..', 'experiment_results');
end

MASS_AW0   = 1.00;
MASS_AW2   = 1.27;
MASS_AW3   = 1.45;
MASS_NOMINAL = MASS_AW2;
ETA_NOMINAL  = 0.85;

LOADS   = struct('name',{'AW0','AW2','AW3'},'mf',{MASS_AW0,MASS_AW2,MASS_AW3});
ETAS    = struct('name',{'eta_lo','eta_nom','eta_hi'},'val',{0.75,0.85,0.90});
SPEEDS  = struct('name',{'spd_nom','spd_-20'},'red',{0,0.20});

if ~exist(OUT_DIR,'dir'), mkdir(OUT_DIR); end

%% ===== GLOBALS + SETUP =====
global use_improved nsga2_variant pop_size iterations dimension
global vel_profile time_obj_max parallel_use show_progress driving_strategy
global Mass lambda inertial_mass Davis

use_improved     = true;
nsga2_variant    = 'rl_sde';
show_progress    = false;
parallel_use     = true;
driving_strategy = "CC_CR";
time_obj_max     = T_TARGET * 1.50;   % lebih longgar untuk phase optimasi

cfg = struct('route_file',ROUTE_IS04,'rollingstock_file',RS_FILE, ...
    'driving_strategy',"CC_CR",'pop_size',POP_SIZE,'iterations',ITERATIONS, ...
    'use_improved',true,'nsga2_variant','rl_sde','parallel_use',true,'sim',struct());

info = setup_project(cfg);
Mass_AW0_nom = Mass;

% Scale to AW2 nominal
Mass          = Mass_AW0_nom * MASS_AW2;
inertial_mass = Mass * (1 + lambda);
fprintf('Nominal: AW2 mass=%.0f t | η=%.2f | T_target=%d s\n', Mass, ETA_NOMINAL, T_TARGET);

setup_parallel_pool(cfg, info);
pool = gcp('nocreate');
if ~isempty(pool)
    wait(parfevalOnAll(@init_worker_globals, 0, info.routePath, info.rollingstockPath, ...
        cfg.sim, true));
    wait(parfevalOnAll(@deal_mass, 0, Mass_AW0_nom * MASS_AW2));
end

[T_min, E_base_AW2] = flat_out_baseline_rk4(info.routePath, info.rollingstockPath);
E_NOMINAL = E_base_AW2 / ETA_NOMINAL;
fprintf('T_min=%.2f s | E_baseline(AW2,η=0.85)=%.4f kWh_elec\n', T_min, E_NOMINAL);

%% ===== STEP 1: GENERATE 30 PARETO FRONTS =====
fprintf('\n--- Step 1: Generating %d Pareto fronts (nominal AW2) ---\n', N_RUNS);
pop_size   = POP_SIZE;
iterations = ITERATIONS;

% Push to workers
if ~isempty(pool)
    wait(parfevalOnAll(@push_worker_globals, 0, POP_SIZE, ITERATIONS, T_TARGET*1.50));
end

all_pops  = cell(N_RUNS,1);
all_tims  = zeros(N_RUNS,1);
all_F_ref = zeros(0,2);

parfor (run_id = 1:N_RUNS, 4)
    rng(run_id,'twister');
    t0_=tic;
    [p,~,~]=nsga2_main(vel_profile);
    all_tims(run_id)=toc(t0_);
    all_pops{run_id}=p;
end

for run_id=1:N_RUNS
    F=double(all_pops{run_id}(:,dimension+1:dimension+2));
    valid=F(:,1)<1e5;
    if any(valid), all_F_ref=[all_F_ref;nd_filter(F(valid,:))]; end %#ok<AGROW>
    fprintf('  run %02d | F1=%d pts | %.1fs\n',run_id, ...
        sum(F(valid,1)<1e5&F(valid,2)<1e5),all_tims(run_id));
end

%% ===== STEP 2: SELECT BEST SOLUTION PER RUN =====
fprintf('\n--- Step 2: Selecting best solution per run ---\n');
F_ref_global = nd_filter(all_F_ref);
mins_r=min(F_ref_global);maxs_r=max(F_ref_global);
rng_r=maxs_r-mins_r+eps;
ref_pt=ones(1,2)+0.1;

cand_X     = zeros(N_RUNS,dimension);
cand_F_nom = zeros(N_RUNS,2);
cand_HV    = zeros(N_RUNS,1);

for run_id=1:N_RUNS
    pop_r=all_pops{run_id};
    F_all=double(pop_r(:,dimension+1:dimension+2));
    X_all=double(pop_r(:,1:dimension));
    valid=F_all(:,1)<1e5;
    F_v=F_all(valid,:); X_v=X_all(valid,:);
    F_nd=nd_filter(F_v);
    nd_idx=nondom_idx(F_v);
    X_nd=X_v(nd_idx,:);

    Fn=(F_nd-mins_r)./rng_r;
    cand_HV(run_id)=hv2d_local(Fn,ref_pt);

    feas_nd=F_nd(:,1)<=T_TARGET;
    if any(feas_nd)
        [~,bi]=min(F_nd(feas_nd,2));
        F_feas=F_nd(feas_nd,:); X_feas=X_nd(feas_nd,:);
    else
        [~,bi]=min(F_nd(:,1));
        F_feas=F_nd; X_feas=X_nd;
    end
    cand_X(run_id,:)    = X_feas(min(bi,end),:);
    cand_F_nom(run_id,:)= F_feas(min(bi,end),:);
end

E_nom_cand = cand_F_nom(:,2) / ETA_NOMINAL;
fprintf('Candidates: mean T=%.1fs, mean E=%.3f kWh_elec\n', ...
    mean(cand_F_nom(:,1),'omitnan'), mean(E_nom_cand,'omitnan'));

%% ===== STEP 3: BUILD SCENARIO LIST + EVALUATE =====
fprintf('\n--- Step 3: Evaluating 18 scenarios ---\n');
scenarios = struct('load_name',{},'eta_name',{},'spd_name',{}, ...
    'mass_factor',{},'eta_traction',{},'speed_reduction',{});

for li=1:numel(LOADS)
    for ei=1:numel(ETAS)
        for vi=1:numel(SPEEDS)
            s=struct('load_name',LOADS(li).name,'eta_name',ETAS(ei).name, ...
                'spd_name',SPEEDS(vi).name, ...
                'mass_factor',LOADS(li).mf,'eta_traction',ETAS(ei).val, ...
                'speed_reduction',SPEEDS(vi).red);
            scenarios(end+1)=s; %#ok<AGROW>
        end
    end
end
n_sc=numel(scenarios);

T_mat = zeros(N_RUNS,n_sc);
E_mat = zeros(N_RUNS,n_sc);
T_feas_limit = T_TARGET * FEAS_TOL;

for sol_id=1:N_RUNS
    X=cand_X(sol_id,:);
    fprintf('  Sol %02d/%d |', sol_id, N_RUNS);
    for sc_id=1:n_sc
        [T_sc,E_sc]=simulate_with_uncertainty(X,scenarios(sc_id));
        T_mat(sol_id,sc_id)=T_sc;
        E_mat(sol_id,sc_id)=E_sc;
        fprintf('.');
    end
    fprintf('|\n');
end

%% ===== STEP 4: ROBUSTNESS METRICS =====
feas_rate    = zeros(N_RUNS,1);
E_degrad_pct = nan(N_RUNS,1);

for sol_id=1:N_RUNS
    feas_sc=T_mat(sol_id,:)<=T_feas_limit;
    feas_rate(sol_id)=sum(feas_sc)/n_sc;
    E_nom=E_nom_cand(sol_id);
    if any(feas_sc) && E_nom>0
        E_degrad_pct(sol_id)=mean( ...
            (E_mat(sol_id,feas_sc)-E_nom)/E_nom*100);
    end
end

robust = feas_rate >= ROBUST_THR;
fprintf('\nRobust solutions (F_rate>=%.0f%%): %d/%d\n', ROBUST_THR*100,sum(robust),N_RUNS);
fprintf('Mean energy degradation (robust): %.2f%%\n', mean(E_degrad_pct(robust),'omitnan'));

%% ===== SAVE CSV =====
rows={};
for sol_id=1:N_RUNS
    for sc_id=1:n_sc
        sc=scenarios(sc_id);
        feas=T_mat(sol_id,sc_id)<=T_feas_limit;
        E_n=E_nom_cand(sol_id);
        deg=ternary(feas&&E_n>0,(E_mat(sol_id,sc_id)-E_n)/E_n*100,NaN);
        rows{end+1}={sol_id,sc.load_name,sc.eta_name,sc.spd_name, ...
            T_mat(sol_id,sc_id),E_mat(sol_id,sc_id),feas,deg}; %#ok<AGROW>
    end
end
T_csv=cell2table(vertcat(rows{:}),'VariableNames', ...
    {'solution_id','load','efficiency','speed','travel_time','energy','feasible','degradation_pct'});
writetable(T_csv, fullfile(OUT_DIR,'E4_results.csv'));
fprintf('Saved: %s\n', fullfile(OUT_DIR,'E4_results.csv'));

%% ===== SUMMARY =====
fid=fopen(fullfile(OUT_DIR,'E4_summary.txt'),'w');
fprintf(fid,'Experiment E4 — Robustness Evaluation\n');
fprintf(fid,'Segment IS04 | Nominal: AW2 mass_factor=%.2f | η=%.2f | T_target=%ds\n', ...
    MASS_AW2,ETA_NOMINAL,T_TARGET);
fprintf(fid,'Scenarios: %d | Feas threshold: T<=T_target×%.2f=%.1fs | Robust: F_rate>=%.0f%%\n\n', ...
    n_sc,FEAS_TOL,T_feas_limit,ROBUST_THR*100);
fprintf(fid,'Robust solutions: %d / %d\n',sum(robust),N_RUNS);
fprintf(fid,'Mean E_degrad (all)    : %.2f%%\n',mean(E_degrad_pct,'omitnan'));
fprintf(fid,'Mean E_degrad (robust) : %.2f%%\n',mean(E_degrad_pct(robust),'omitnan'));
fprintf(fid,'\n%-8s  %8s  %10s  %8s\n','Sol_ID','F_rate','E_degrad%%','Robust');
for sol_id=1:N_RUNS
    fprintf(fid,'%-8d  %8.3f  %10.2f  %8s\n',sol_id,feas_rate(sol_id), ...
        E_degrad_pct(sol_id),ternary(robust(sol_id),'YES','no'));
end
fclose(fid);

%% ===== PLOTS =====
figure(4001);clf;
deg_mat=zeros(N_RUNS,n_sc);
for sol_id=1:N_RUNS
    E_n=E_nom_cand(sol_id);
    if E_n>0, deg_mat(sol_id,:)=(E_mat(sol_id,:)-E_n)/E_n*100; end
end
sc_labels=arrayfun(@(s)sprintf('%s\n%s\n%s',s.load_name,s.eta_name,s.spd_name), ...
    scenarios,'UniformOutput',false);
boxplot(deg_mat,'Labels',sc_labels);
ylabel('Energy Degradation (%)');
title('E4: Energy Degradation Under 18 Uncertainty Scenarios');
xtickangle(45); grid on;
saveas(gcf, fullfile(OUT_DIR,'E4_energy_degradation.png'));

figure(4002);clf;
b=bar(1:N_RUNS, feas_rate, 0.7);
hold on;
b.FaceColor='flat';
for si=1:N_RUNS
    b.CData(si,:)=ternary(robust(si),[0.3 0.7 0.3],[0.9 0.4 0.4]);
end
yline(ROBUST_THR,'--k','Robust threshold','LineWidth',1.5);
xlabel('Candidate Solution'); ylabel('Feasibility Rate (F_{rate})');
title(sprintf('E4: F_{rate} per Candidate — %d scenarios | IS04 AW2 nominal',n_sc));
ylim([0 1.1]); grid on;
saveas(gcf, fullfile(OUT_DIR,'E4_feasibility_rate.png'));

fprintf('\nE4 complete → %s\n',OUT_DIR);

%% ===== HELPERS =====
function idx=nondom_idx(F)
    n=size(F,1); keep=true(n,1);
    for i=1:n
        if ~keep(i),continue;end
        for j=1:n
            if i==j||~keep(j),continue;end
            if all(F(j,:)<=F(i,:))&&any(F(j,:)<F(i,:)),keep(i)=false;break;end
        end
    end
    idx=find(keep);
end

function F_nd=nd_filter(F)
    if isempty(F),F_nd=zeros(0,2);return;end
    n=size(F,1);k=true(n,1);
    for i=1:n
        if ~k(i),continue;end
        for j=1:n
            if i==j||~k(j),continue;end
            if all(F(j,:)<=F(i,:))&&any(F(j,:)<F(i,:)),k(i)=false;break;end
        end
    end
    F_nd=F(k,:);[~,o]=sort(F_nd(:,1));F_nd=F_nd(o,:);
end

function hv=hv2d_local(F,ref)
    if isempty(F),hv=0;return;end
    n=size(F,1);k=true(n,1);
    for i=1:n,if~k(i),continue;end,for j=1:n,if i==j||~k(j),continue;end
    if all(F(j,:)<=F(i,:))&&any(F(j,:)<F(i,:)),k(i)=false;break;end,end,end
    F=F(k,:);[~,o]=sort(F(:,1));F=F(o,:);F(:,2)=cummin(F(:,2));
    hv=0;pE=ref(2);
    for i=1:size(F,1),if F(i,2)>pE,continue;end,hv=hv+(ref(1)-F(i,1))*(pE-F(i,2));pE=F(i,2);end
end

function v=ternary(c,a,b),if c,v=a;else,v=b;end,end

function deal_mass(m),global Mass lambda inertial_mass;Mass=m;inertial_mass=Mass*(1+lambda);end
