clear all
close all
clc


%% ── 0. Operating point ───────────────────────────────────────────────────
load TCL_MechModel_Parameters.mat   % provides TCL struct (Ys, Us, Samp_T, Xs)
Ys = TCL.Ys(:);
Us = TCL.Us(:);
Ts = TCL.Samp_T;   % 4 s

k_step1   = 101;
k_step2   = 301;
step_locs = [k_step1-1, k_step2-1];

%%  SECTION 1 ── CONTROLLER / DATASET REGISTRY

% ── Controller colours (shared between sim & exp) ────────────────────────
CTRL_COLORS = { [0.15 0.35 0.75];   % 1  PI      – blue
                [0.80 0.15 0.15];   % 2  MPC-BB  – red
                [0.10 0.60 0.10];   % 3  MPC-Mech– green
                [0.80 0.40 0.00] }; % 4  MPC-OT  – orange

% ────────────────────────────────────────────────────────────────────────
% SIM registry  (PATH A only: 2xN matrices already in .mat files)
% ────────────────────────────────────────────────────────────────────────
SIM(1).label          = 'PI';
SIM(1).file           = 'TCL_PI_Servo_SimResults.mat';
SIM(1).Yk_var         = 'Yk';
SIM(1).Uk_var         = 'Uk';
SIM(1).Rk_f_var       = 'Rk_f';
SIM(1).ek_f_var       = 'ek_f';
SIM(1).has_ss_targets = false;
SIM(1).Xs_var         = '';
SIM(1).Us_var         = '';

SIM(2).label          = 'MPC-BB';
SIM(2).file           = 'TCL_MPC_BB_SimResults_NoiseFree.mat';
SIM(2).Yk_var         = 'Yk';
SIM(2).Uk_var         = 'Uk';
SIM(2).Rk_f_var       = 'Rk_f';
SIM(2).ek_f_var       = 'ek_f';
SIM(2).has_ss_targets = true;
SIM(2).Xs_var         = 'Xs_k';
SIM(2).Us_var         = 'Us_k';

SIM(3).label          = 'MPC-Mech';
SIM(3).file           = 'TCL_MPC_Mech_SimResults_NoiseFree.mat';
SIM(3).Yk_var         = 'Yk';
SIM(3).Uk_var         = 'Uk';
SIM(3).Rk_f_var       = 'Rk_f';
SIM(3).ek_f_var       = 'ek_f';
SIM(3).has_ss_targets = true;
SIM(3).Xs_var         = 'Xs_k';
SIM(3).Us_var         = 'Us_k';

SIM(4).label          = 'MPC-OT';
SIM(4).file           = 'TCL_MPC_OT_SimResults_NoiseFree.mat';
SIM(4).Yk_var         = 'Yk_ot';
SIM(4).Uk_var         = 'Uk_ot';
SIM(4).Rk_f_var       = 'Rk_f_ot';
SIM(4).ek_f_var       = 'ek_f_ot';
SIM(4).has_ss_targets = true;
SIM(4).Xs_var         = 'Xs_k_ot';
SIM(4).Us_var         = 'Us_k_ot';

% ────────────────────────────────────────────────────────────────────────
% EXP registry
%   PATH A = variable is already a 2xN matrix  (Yk_var non-empty)
%   PATH B = assemble from individual scalar vectors (Yk_var empty)
% ────────────────────────────────────────────────────────────────────────
EXP(1).label          = 'PI';
EXP(1).file           = 'TCL_PI_Servo_ExpResults.mat';
EXP(1).Yk_var         = '';          % PATH B
EXP(1).Uk_var         = '';
EXP(1).Rk_f_var       = '';
EXP(1).ek_f_var       = '';
EXP(1).t1_var         = 't1s';
EXP(1).t2_var         = 't2s';
EXP(1).h1_var         = 'h1s';
EXP(1).h2_var         = 'h2s';
EXP(1).R1_var         = 'R1s';
EXP(1).R2_var         = 'R2s';
EXP(1).e1_var         = 'e1s';
EXP(1).e2_var         = 'e2s';
EXP(1).kT_var         = '';
EXP(1).has_ss_targets = false;
EXP(1).xs_vars        = {};
EXP(1).us_vars        = {};
EXP(1).xs_absolute    = false;
EXP(1).us_absolute    = false;

EXP(2).label          = 'MPC-BB';
EXP(2).file           = 'TCL_MPC_BB_ExpResults.mat';
EXP(2).Yk_var         = '';          % PATH B
EXP(2).Uk_var         = '';
EXP(2).Rk_f_var       = '';
EXP(2).ek_f_var       = '';
EXP(2).t1_var         = 't1s';
EXP(2).t2_var         = 't2s';
EXP(2).h1_var         = 'h1s';
EXP(2).h2_var         = 'h2s';
EXP(2).R1_var         = 'R1s';
EXP(2).R2_var         = 'R2s';
EXP(2).e1_var         = 'e1s';
EXP(2).e2_var         = 'e2s';
EXP(2).kT_var         = 'kT';
EXP(2).has_ss_targets = true;
EXP(2).xs_vars        = {'xs1s', 'xs2s'};
EXP(2).us_vars        = {'us1s', 'us2s'};
EXP(2).xs_absolute    = false;
EXP(2).us_absolute    = false;

EXP(3).label          = 'MPC-Mech';
EXP(3).file           = 'TCL_MPC_Mech_ExpResults.mat';
EXP(3).Yk_var         = 'Yk';       % PATH A
EXP(3).Uk_var         = 'Uk';
EXP(3).Rk_f_var       = 'Rk_f';
EXP(3).ek_f_var       = 'ek_f';
EXP(3).t1_var         = '';
EXP(3).t2_var         = '';
EXP(3).h1_var         = '';
EXP(3).h2_var         = '';
EXP(3).R1_var         = '';
EXP(3).R2_var         = '';
EXP(3).e1_var         = '';
EXP(3).e2_var         = '';
EXP(3).kT_var         = 'kT';
EXP(3).has_ss_targets = true;
EXP(3).xs_vars        = {'Xs_k'};
EXP(3).us_vars        = {'Us_k'};
EXP(3).xs_absolute    = false;
EXP(3).us_absolute    = false;

EXP(4).label          = 'MPC-OT';
EXP(4).file           = 'TCL_MPC_OT_ExpResults.mat';
EXP(4).Yk_var         = 'Yk_ot';    % PATH A
EXP(4).Uk_var         = 'Uk_ot';
EXP(4).Rk_f_var       = 'Rk_f_ot';
EXP(4).ek_f_var       = 'ek_f_ot';
EXP(4).t1_var         = '';
EXP(4).t2_var         = '';
EXP(4).h1_var         = '';
EXP(4).h2_var         = '';
EXP(4).R1_var         = '';
EXP(4).R2_var         = '';
EXP(4).e1_var         = '';
EXP(4).e2_var         = '';
EXP(4).kT_var         = 'kT';
EXP(4).has_ss_targets = true;
EXP(4).xs_vars        = {'Xs_k_ot'};
EXP(4).us_vars        = {'Us_k_ot'};
EXP(4).xs_absolute    = true;
EXP(4).us_absolute    = true;

nC = length(SIM);   % same for EXP

%%  SECTION 2 ── LOAD DATA  (helper inline functions below)

fprintf('\n── Loading SIMULATION data ──────────────────────────────────\n');
for c = 1:nC
    fprintf('  [SIM %d] %s  (%s)\n', c, SIM(c).label, SIM(c).file);
    raw = load(SIM(c).file);

    SIM(c).Yk   = raw.(SIM(c).Yk_var);
    SIM(c).Uk   = raw.(SIM(c).Uk_var);
    SIM(c).Rk_f = raw.(SIM(c).Rk_f_var);
    SIM(c).ek_f = raw.(SIM(c).ek_f_var);
    SIM(c).kT   = (0 : size(SIM(c).Yk, 2) - 1)';

    if SIM(c).has_ss_targets && ~isempty(SIM(c).Xs_var) && isfield(raw, SIM(c).Xs_var)
        SIM(c).Xs_k = raw.(SIM(c).Xs_var);
        SIM(c).Us_k = raw.(SIM(c).Us_var);
    else
        SIM(c).Xs_k = [];
        SIM(c).Us_k = [];
    end

    SIM(c).N = length(SIM(c).kT);
    fprintf('     %d samples  (%.1f min)\n', SIM(c).N, SIM(c).N*Ts/60);
end

fprintf('\n── Loading EXPERIMENT data ───────────────────────────────────\n');
for c = 1:nC
    fprintf('  [EXP %d] %s  (%s)\n', c, EXP(c).label, EXP(c).file);
    raw = load(EXP(c).file);

    % --- Yk ---------------------------------------------------------------
    if ~isempty(EXP(c).Yk_var) && isfield(raw, EXP(c).Yk_var)
        EXP(c).Yk = raw.(EXP(c).Yk_var);                   % PATH A
    else
        EXP(c).Yk = [raw.(EXP(c).t1_var)(:)'; ...
                     raw.(EXP(c).t2_var)(:)'];              % PATH B
    end

    % --- Uk ---------------------------------------------------------------
    if ~isempty(EXP(c).Uk_var) && isfield(raw, EXP(c).Uk_var)
        EXP(c).Uk = raw.(EXP(c).Uk_var);
    else
        EXP(c).Uk = [raw.(EXP(c).h1_var)(:)'; ...
                     raw.(EXP(c).h2_var)(:)'];
    end

    % --- Rk_f -------------------------------------------------------------
    if ~isempty(EXP(c).Rk_f_var) && isfield(raw, EXP(c).Rk_f_var)
        EXP(c).Rk_f = raw.(EXP(c).Rk_f_var);
    else
        EXP(c).Rk_f = [raw.(EXP(c).R1_var)(:)'; ...
                       raw.(EXP(c).R2_var)(:)'];
    end

    % --- ek_f -------------------------------------------------------------
    if ~isempty(EXP(c).ek_f_var) && isfield(raw, EXP(c).ek_f_var)
        EXP(c).ek_f = raw.(EXP(c).ek_f_var);
    else
        EXP(c).ek_f = [raw.(EXP(c).e1_var)(:)'; ...
                       raw.(EXP(c).e2_var)(:)'];
    end

    % --- kT ---------------------------------------------------------------
    if ~isempty(EXP(c).kT_var) && isfield(raw, EXP(c).kT_var)
        EXP(c).kT = raw.(EXP(c).kT_var)(:);
    else
        EXP(c).kT = (0 : size(EXP(c).Yk, 2) - 1)';
    end

    % --- SS targets -------------------------------------------------------
    if EXP(c).has_ss_targets && ~isempty(EXP(c).xs_vars)
        xs_rows = {};
        for vi = 1:length(EXP(c).xs_vars)
            v = raw.(EXP(c).xs_vars{vi});
            if min(size(v)) > 1
                xs_rows{end+1} = v;       %#ok
            else
                xs_rows{end+1} = v(:)';   %#ok
            end
        end
        xs_raw = vertcat(xs_rows{:});

        us_rows = {};
        for vi = 1:length(EXP(c).us_vars)
            v = raw.(EXP(c).us_vars{vi});
            if min(size(v)) > 1
                us_rows{end+1} = v;       %#ok
            else
                us_rows{end+1} = v(:)';   %#ok
            end
        end
        us_raw = vertcat(us_rows{:});

        N_ss = size(xs_raw, 2);
        if EXP(c).xs_absolute
            EXP(c).Xs_k = xs_raw;
        else
            EXP(c).Xs_k = xs_raw + repmat(TCL.Xs(:), 1, N_ss);
        end
        if EXP(c).us_absolute
            EXP(c).Us_k = us_raw;
        else
            EXP(c).Us_k = us_raw + repmat(Us, 1, N_ss);
        end
    else
        EXP(c).Xs_k = [];
        EXP(c).Us_k = [];
    end

    EXP(c).N = length(EXP(c).kT);
    fprintf('     %d samples  (%.1f min)\n', EXP(c).N, EXP(c).N*Ts/60);
end

% Comparison window: shortest across all sim & exp datasets
Ns = min([SIM.N, EXP.N]);
kT = (0:Ns-1)';
fprintf('\nComparison window  Ns = %d samples  (%.1f min)\n\n', Ns, Ns*Ts/60);

% Reference for plotting (use SIM controller 1 as the common reference)
Rk_f_ref = SIM(1).Rk_f(:, 1:Ns);

%%  SECTION 3 ── PERFORMANCE INDICES

for c = 1:nC
    for src = 1:2   % 1 = SIM, 2 = EXP
        if src == 1, D = SIM(c); else, D = EXP(c); end
        Yk_c = D.Yk(:, 1:Ns);
        Uk_c = D.Uk(:, 1:Ns);
        Rk_c = D.Rk_f(:, 1:Ns);
        sse  = zeros(1,2);  ssmv = zeros(1,2);  ssdel = zeros(1,2);
        for i = 1:2
            sse(i)  = sum((Yk_c(i,:) - Rk_c(i,:)).^2);
            ssmv(i) = sum((Uk_c(i,:) - Us(i)).^2);
            ssdel(i)= sum(diff(Uk_c(i,:)).^2);
        end
        if src == 1
            SIM(c).SSE = sse;  SIM(c).SSMV = ssmv;  SIM(c).SSdelMV = ssdel;
        else
            EXP(c).SSE = sse;  EXP(c).SSMV = ssmv;  EXP(c).SSdelMV = ssdel;
        end
    end
end

% ── Print side-by-side table ─────────────────────────────────────────────
idx_rows  = {'SSE (Loop 1)', 'SSE (Loop 2)', ...
             'SSMV (Loop 1)', 'SSMV (Loop 2)', ...
             'SSdeltaMV (L1)', 'SSdeltaMV (L2)'};
idx_fields= {'SSE','SSE','SSMV','SSMV','SSdelMV','SSdelMV'};
idx_loops = [1 2 1 2 1 2];

col_w = 14;
sep   = repmat('=', 1, 22 + col_w*2*nC);
sep2  = repmat('-', 1, 22 + col_w*2*nC);

fprintf('%s\n', sep);
fprintf('  SIMULATION vs EXPERIMENT – PERFORMANCE INDICES  (Ns=%d, %.1f min)\n', Ns, Ns*Ts/60);
fprintf('%s\n', sep);

% Header row 1: controller names
hdr1 = sprintf('  %-20s', 'Index');
for c = 1:nC
    hdr1 = [hdr1 sprintf('| %-27s', sprintf('── %s ──', SIM(c).label))]; %#ok
end
fprintf('%s\n', hdr1);

% Header row 2: Sim / Exp columns
hdr2 = sprintf('  %-20s', '');
for c = 1:nC
    hdr2 = [hdr2 sprintf('| %-13s %-13s ', 'Sim', 'Exp')]; %#ok
end
fprintf('%s\n%s\n', hdr2, sep2);

for r = 1:6
    line = sprintf('  %-20s', idx_rows{r});
    for c = 1:nC
        sv = SIM(c).(idx_fields{r})(idx_loops(r));
        ev = EXP(c).(idx_fields{r})(idx_loops(r));
        line = [line sprintf('| %13.3f %13.3f ', sv, ev)]; %#ok
    end
    fprintf('%s\n', line);
end
fprintf('%s\n\n', sep);

%%  SECTION 4 ── PLOTS

% Linestyle convention:
%   Simulation  →  solid   ( '-'  )
%   Experiment  →  dashed  ( '--' )
% Each controller keeps its colour.

% ── Figure 1: Output Y_i(k) – Sim vs Exp, one sub-plot per controller ────
figure('Name','SimExp – Output Profiles','NumberTitle','off', ...
       'Position',[30 50 1200 800])

for c = 1:nC
    col = CTRL_COLORS{c};
    for i = 1:2
        subplot(nC, 2, (c-1)*2 + i), hold on
        % Simulation
        plot(kT, SIM(c).Yk(i,1:Ns), '-', ...
             'Color', col, 'LineWidth', 2.0, 'DisplayName', 'Sim');
        % Experiment
        plot(kT, EXP(c).Yk(i,1:Ns), '--', ...
             'Color', col*0.65, 'LineWidth', 1.6, 'DisplayName', 'Exp');
        % Reference
        plot(kT, Rk_f_ref(i,:), 'k:', 'LineWidth', 1.2, 'DisplayName', 'Ref');
        xline(step_locs(1),'k:','HandleVisibility','off')
        xline(step_locs(2),'k:','HandleVisibility','off')
        grid on
        ylabel(sprintf('T_%d  (°C)', i), 'FontSize', 10)
        if i == 1
            title(sprintf('%s  –  Output Y_i(k)', SIM(c).label), 'FontSize', 11)
        end
        if c == nC, xlabel('Sample  k', 'FontSize', 10); end
        legend('Location','Best','FontSize',8)
    end
end
sgtitle('Output Profiles Y_i(k) – Simulation (solid) vs Experiment (dashed)', ...
        'FontSize', 13)

% ── Figure 2: Heater U_i(k) – Sim vs Exp ────────────────────────────────
figure('Name','SimExp – Heater Inputs','NumberTitle','off', ...
       'Position',[60 50 1200 800])

for c = 1:nC
    col = CTRL_COLORS{c};
    for i = 1:2
        subplot(nC, 2, (c-1)*2 + i), hold on
        stairs(kT, SIM(c).Uk(i,1:Ns)', '-',  'Color', col,      'LineWidth', 2.0, 'DisplayName', 'Sim');
        stairs(kT, EXP(c).Uk(i,1:Ns)', '--', 'Color', col*0.65, 'LineWidth', 1.6, 'DisplayName', 'Exp');
        yline(5,  'k--', 'LineWidth', 0.8, 'HandleVisibility','off')
        yline(80, 'k--', 'LineWidth', 0.8, 'HandleVisibility','off')
        xline(step_locs(1),'k:','HandleVisibility','off')
        xline(step_locs(2),'k:','HandleVisibility','off')
        grid on
        ylim([0 105])
        ylabel(sprintf('U_%d(k)  (%% heater)', i), 'FontSize', 10)
        if i == 1
            title(sprintf('%s  –  Input U_i(k)', SIM(c).label), 'FontSize', 11)
        end
        if c == nC, xlabel('Sample  k', 'FontSize', 10); end
        legend('Location','Best','FontSize',8)
    end
end
sgtitle('Heater Inputs U_i(k) – Simulation (solid) vs Experiment (dashed)', ...
        'FontSize', 13)

% ── Figure 3: Filtered error e_{f,i}(k) – Sim vs Exp ────────────────────
figure('Name','SimExp – Filtered Error','NumberTitle','off', ...
       'Position',[90 50 1200 800])

for c = 1:nC
    col = CTRL_COLORS{c};
    for i = 1:2
        subplot(nC, 2, (c-1)*2 + i), hold on
        stairs(kT, SIM(c).ek_f(i,1:Ns)', '-',  'Color', col,      'LineWidth', 2.0, 'DisplayName', 'Sim');
        stairs(kT, EXP(c).ek_f(i,1:Ns)', '--', 'Color', col*0.65, 'LineWidth', 1.6, 'DisplayName', 'Exp');
        yline(0, 'k--', 'LineWidth', 0.8, 'HandleVisibility','off')
        xline(step_locs(1),'k:','HandleVisibility','off')
        xline(step_locs(2),'k:','HandleVisibility','off')
        grid on
        ylabel(sprintf('e_{f,%d}(k)  (°C)', i), 'FontSize', 10)
        if i == 1
            title(sprintf('%s  –  Error e_{f,i}(k)', SIM(c).label), 'FontSize', 11)
        end
        if c == nC, xlabel('Sample  k', 'FontSize', 10); end
        legend('Location','Best','FontSize',8)
    end
end
sgtitle('Filtered Errors e_{f,i}(k) – Simulation (solid) vs Experiment (dashed)', ...
        'FontSize', 13)

% ── Figure 4: All-controller overlay – Outputs  (Sim + Exp together) ─────
figure('Name','SimExp – Output Overlay (all controllers)','NumberTitle','off', ...
       'Position',[120 50 1000 700])

for i = 1:2
    subplot(2,1,i), hold on
    for c = 1:nC
        col = CTRL_COLORS{c};
        lbl = SIM(c).label;
        plot(kT, SIM(c).Yk(i,1:Ns), '-',  'Color', col,      'LineWidth', 1.8, ...
             'DisplayName', sprintf('Sim %s', lbl));
        plot(kT, EXP(c).Yk(i,1:Ns), '--', 'Color', col*0.65, 'LineWidth', 1.4, ...
             'DisplayName', sprintf('Exp %s', lbl));
    end
    plot(kT, Rk_f_ref(i,:), 'k:', 'LineWidth', 1.5, 'DisplayName', 'Reference')
    xline(step_locs(1),'k:','Step 1','LabelVerticalAlignment','bottom', ...
          'FontSize',8,'HandleVisibility','off')
    xline(step_locs(2),'k:','Step 2','LabelVerticalAlignment','bottom', ...
          'FontSize',8,'HandleVisibility','off')
    grid on
    ylabel(sprintf('T_%d  (°C)', i), 'FontSize', 11)
    if i == 1
        title('Output Y_i(k) – All Controllers: Sim (solid) vs Exp (dashed)', 'FontSize', 12)
    end
    if i == 2, xlabel('Sample  k', 'FontSize', 11); end
    legend('Location','Best','FontSize',8, 'NumColumns', 2)
end

% ── Figure 5: All-controller overlay – Inputs ────────────────────────────
figure('Name','SimExp – Input Overlay (all controllers)','NumberTitle','off', ...
       'Position',[150 50 1000 700])

for i = 1:2
    subplot(2,1,i), hold on
    for c = 1:nC
        col = CTRL_COLORS{c};
        lbl = SIM(c).label;
        stairs(kT, SIM(c).Uk(i,1:Ns)', '-',  'Color', col,      'LineWidth', 1.8, ...
               'DisplayName', sprintf('Sim %s', lbl));
        stairs(kT, EXP(c).Uk(i,1:Ns)', '--', 'Color', col*0.65, 'LineWidth', 1.4, ...
               'DisplayName', sprintf('Exp %s', lbl));
    end
    yline(5,  'k--', 'LineWidth', 0.9, 'HandleVisibility','off')
    yline(80, 'k--', 'LineWidth', 0.9, 'HandleVisibility','off')
    xline(step_locs(1),'k:','HandleVisibility','off')
    xline(step_locs(2),'k:','HandleVisibility','off')
    grid on
    ylim([0 105])
    ylabel(sprintf('U_%d(k)  (%% heater)', i), 'FontSize', 11)
    if i == 1, title('Heater Input U_i(k) – All Controllers: Sim (solid) vs Exp (dashed)', 'FontSize', 12); end
    if i == 2, xlabel('Sample  k', 'FontSize', 11); end
    legend('Location','Best','FontSize',8, 'NumColumns', 2)
end

% ── Figure 6: SS targets (Experiment, MPC controllers only) ──────────────
ss_idx = find([EXP.has_ss_targets]);

if ~isempty(ss_idx)
    n_ss = length(ss_idx);
    figure('Name','Exp SS Targets xs(k) and us(k)','NumberTitle','off', ...
           'Position',[200 50 1100 300*n_ss])

    for row = 1:n_ss
        c = ss_idx(row);
        if isempty(EXP(c).Xs_k), continue; end

        Ns_c   = min(Ns, EXP(c).N);
        Xs_c   = EXP(c).Xs_k(:, 1:Ns_c);
        Us_c   = EXP(c).Us_k(:, 1:Ns_c);
        kT_c   = EXP(c).kT(1:Ns_c);
        n_st_c = size(Xs_c, 1);
        col    = CTRL_COLORS{c};

        subplot(n_ss, 2, 2*row-1), hold on
        for s = 1:n_st_c
            plot(kT_c, Xs_c(s,:), 'LineWidth', 1.5)
        end
        xline(step_locs(1),'k:'); xline(step_locs(2),'k:')
        grid on
        ylabel('x_{s,i}(k)  (°C)', 'FontSize', 10)
        title(sprintf('%s [Exp] – SS Target States', EXP(c).label), 'FontSize', 11)
        lx = arrayfun(@(s) sprintf('x_{s,%d}(k)',s), 1:n_st_c, 'UniformOutput', false);
        legend(lx{:}, 'Location','Best','FontSize',9)
        if row == n_ss, xlabel('Sample  k', 'FontSize', 10); end

        subplot(n_ss, 2, 2*row), hold on
        stairs(kT_c, Us_c(1,:)', 'Color', col,      'LineWidth', 1.5)
        stairs(kT_c, Us_c(2,:)', 'Color', col*0.65, 'LineWidth', 1.5)
        xline(step_locs(1),'k:'); xline(step_locs(2),'k:')
        grid on
        ylabel('u_{s,i}(k)  (% heater)', 'FontSize', 10)
        title(sprintf('%s [Exp] – SS Target Inputs', EXP(c).label), 'FontSize', 11)
        legend('u_{s,1}(k)','u_{s,2}(k)', 'Location','Best','FontSize',9)
        if row == n_ss, xlabel('Sample  k', 'FontSize', 10); end
    end
end

% ── Figure 7: Bar-chart – Sim vs Exp per index per controller ────────────
figure('Name','Performance Index Bar Chart – Sim vs Exp','NumberTitle','off', ...
       'Position',[250 50 1200 600])

idx_labels = {'SSE','SSMV','SS\DeltaMV'};
idx_fields2= {'SSE','SSMV','SSdelMV'};
group = categorical(idx_labels);
group = reordercats(group, idx_labels);

for sp = 1:2          % sp=1: Loop 1, sp=2: Loop 2
    subplot(1,2,sp), hold on

    % Build data matrix: rows = 3 indices, cols = 2*nC (Sim1,Exp1,Sim2,Exp2,…)
    data = zeros(3, 2*nC);
    bar_labels = cell(1, 2*nC);
    bar_colors = zeros(2*nC, 3);

    for c = 1:nC
        col_s = 2*c - 1;   % Sim column
        col_e = 2*c;        % Exp column
        for r = 1:3
            data(r, col_s) = SIM(c).(idx_fields2{r})(sp);
            data(r, col_e) = EXP(c).(idx_fields2{r})(sp);
        end
        bar_labels{col_s} = sprintf('Sim\n%s', SIM(c).label);
        bar_labels{col_e} = sprintf('Exp\n%s', EXP(c).label);
        bar_colors(col_s,:) = CTRL_COLORS{c};
        bar_colors(col_e,:) = CTRL_COLORS{c} * 0.6;
    end

    b = bar(group, data', 'grouped');
    for col_b = 1:2*nC
        b(col_b).FaceColor = bar_colors(col_b,:);
        b(col_b).DisplayName = bar_labels{col_b};
    end
    grid on
    ylabel('Index Value', 'FontSize', 10)
    if sp == 1
        title('Loop 1  (T_1 / H_1)', 'FontSize', 11)
    else
        title('Loop 2  (T_2 / H_2)', 'FontSize', 11)
    end
    legend('Location','Best','FontSize',8)
end
sgtitle('Performance Indices – Simulation vs Experiment', 'FontSize', 13)

%%  SECTION 5 ── SAVE

results.Ns = Ns;
results.kT = kT;
results.Ts = Ts;

for c = 1:nC
    fn = matlab.lang.makeValidName(SIM(c).label);

    results.(fn).label = SIM(c).label;

    % Simulation
    results.(fn).Sim.Yk      = SIM(c).Yk(:,1:Ns);
    results.(fn).Sim.Uk      = SIM(c).Uk(:,1:Ns);
    results.(fn).Sim.Rk_f    = SIM(c).Rk_f(:,1:Ns);
    results.(fn).Sim.ek_f    = SIM(c).ek_f(:,1:Ns);
    results.(fn).Sim.SSE     = SIM(c).SSE;
    results.(fn).Sim.SSMV    = SIM(c).SSMV;
    results.(fn).Sim.SSdelMV = SIM(c).SSdelMV;

    % Experiment
    results.(fn).Exp.Yk      = EXP(c).Yk(:,1:Ns);
    results.(fn).Exp.Uk      = EXP(c).Uk(:,1:Ns);
    results.(fn).Exp.Rk_f    = EXP(c).Rk_f(:,1:Ns);
    results.(fn).Exp.ek_f    = EXP(c).ek_f(:,1:Ns);
    results.(fn).Exp.SSE     = EXP(c).SSE;
    results.(fn).Exp.SSMV    = EXP(c).SSMV;
    results.(fn).Exp.SSdelMV = EXP(c).SSdelMV;

    if EXP(c).has_ss_targets && ~isempty(EXP(c).Xs_k)
        Ns_c = min(Ns, EXP(c).N);
        results.(fn).Exp.Xs_k = EXP(c).Xs_k(:,1:Ns_c);
        results.(fn).Exp.Us_k = EXP(c).Us_k(:,1:Ns_c);
    end
end

save TCL_Compare_SimExpResults  results
fprintf('\nResults saved to  TCL_Compare_SimExpResults.mat\n');
fprintf('All done.\n');