
clear all
close all
clc


%% ____  0. OPERATING POINT _____

load TCL_MechModel_Parameters.mat   % provides TCL struct (Ys, Us, Samp_T, Xs)
Ys = TCL.Ys(:);
Us = TCL.Us(:);
Ts = TCL.Samp_T;   % 4 s

k_step1   = 101;
k_step2   = 301;
step_locs = [k_step1-1, k_step2-1];

%  SECTION 1 -- CONTROLLER / DATASET REGISTRY

% -- Controller 1: Decentralised PI ---------------------------------

CONTROLLERS(1).label          = 'PI';
CONTROLLERS(1).file           = 'TCL_PI_Servo_ExpResults.mat';
CONTROLLERS(1).Yk_var         = '';       % PATH B: assemble from scalar vecs
CONTROLLERS(1).Uk_var         = '';
CONTROLLERS(1).Rk_f_var       = '';
CONTROLLERS(1).ek_f_var       = '';
CONTROLLERS(1).t1_var         = 't1s';
CONTROLLERS(1).t2_var         = 't2s';
CONTROLLERS(1).h1_var         = 'h1s';
CONTROLLERS(1).h2_var         = 'h2s';
CONTROLLERS(1).R1_var         = 'R1s';
CONTROLLERS(1).R2_var         = 'R2s';
CONTROLLERS(1).e1_var         = 'e1s';
CONTROLLERS(1).e2_var         = 'e2s';
CONTROLLERS(1).kT_var         = '';       % not saved -> auto-generate
CONTROLLERS(1).has_ss_targets = false;
CONTROLLERS(1).xs_vars        = {};
CONTROLLERS(1).us_vars        = {};
CONTROLLERS(1).xs_absolute    = false;
CONTROLLERS(1).us_absolute    = false;
CONTROLLERS(1).color          = [0.15 0.35 0.75];   % blue

% -- Controller 2: Black-Box MPC ----------------------------------------
CONTROLLERS(2).label          = 'MPC\_BB';
CONTROLLERS(2).file           = 'TCL_MPC_BB_ExpResults.mat';
CONTROLLERS(2).Yk_var         = '';       % PATH B
CONTROLLERS(2).Uk_var         = '';
CONTROLLERS(2).Rk_f_var       = '';
CONTROLLERS(2).ek_f_var       = '';
CONTROLLERS(2).t1_var         = 't1s';
CONTROLLERS(2).t2_var         = 't2s';
CONTROLLERS(2).h1_var         = 'h1s';
CONTROLLERS(2).h2_var         = 'h2s';
CONTROLLERS(2).R1_var         = 'R1s';
CONTROLLERS(2).R2_var         = 'R2s';
CONTROLLERS(2).e1_var         = 'e1s';
CONTROLLERS(2).e2_var         = 'e2s';
CONTROLLERS(2).kT_var         = 'kT';    % row vector 0:N-1
CONTROLLERS(2).has_ss_targets = true;
CONTROLLERS(2).xs_vars        = {'xs1s', 'xs2s'};  % perturbation scalars
CONTROLLERS(2).us_vars        = {'us1s', 'us2s'};  % perturbation scalars
CONTROLLERS(2).xs_absolute    = false;   % add TCL.Xs to convert to absolute
CONTROLLERS(2).us_absolute    = false;   % add Us to convert to absolute
CONTROLLERS(2).color          = [0.80 0.15 0.15];   % red

% -- Controller 3: Mechanistic MPC ---------------------------------------
CONTROLLERS(3).label          = 'MPC\_Mech';
CONTROLLERS(3).file           = 'TCL_MPC_Mech_ExpResults.mat';
CONTROLLERS(3).Yk_var         = 'Yk';       % PATH B
CONTROLLERS(3).Uk_var         = 'Uk';
CONTROLLERS(3).Rk_f_var       = 'Rk_f';
CONTROLLERS(3).ek_f_var       = 'ek_f';
CONTROLLERS(3).t1_var         = 't_sec';
CONTROLLERS(3).t2_var         = 't_sec';
CONTROLLERS(3).h1_var         = '';
CONTROLLERS(3).h2_var         = '';
CONTROLLERS(3).R1_var         = '';
CONTROLLERS(3).R2_var         = '';
CONTROLLERS(3).e1_var         = '';
CONTROLLERS(3).e2_var         = '';
CONTROLLERS(3).kT_var         = 'kT';
CONTROLLERS(3).has_ss_targets = true;
CONTROLLERS(3).xs_vars        = {'Xs_k'};  % perturbation scalars
CONTROLLERS(3).us_vars        = {'Us_k'};  % perturbation scalars
CONTROLLERS(3).xs_absolute    = false;   % add TCL.Xs to convert to absolute
CONTROLLERS(3).us_absolute    = false;   % add Us to convert to absolute
CONTROLLERS(3).color          = [0.10 0.60 0.10];   % green

% -- Controller 4: Output-Tracking MPC -----------------------------------

CONTROLLERS(4).label          = 'MPC\_OT';
CONTROLLERS(4).file           = 'TCL_MPC_OT_ExpResults.mat';
CONTROLLERS(4).Yk_var         = 'Yk_ot';      % PATH A: 2xN matrix
CONTROLLERS(4).Uk_var         = 'Uk_ot';
CONTROLLERS(4).Rk_f_var       = 'Rk_f_ot';
CONTROLLERS(4).ek_f_var       = 'ek_f_ot';
CONTROLLERS(4).t1_var         = '';            % not used (PATH A)
CONTROLLERS(4).t2_var         = '';
CONTROLLERS(4).h1_var         = '';
CONTROLLERS(4).h2_var         = '';
CONTROLLERS(4).R1_var         = '';
CONTROLLERS(4).R2_var         = '';
CONTROLLERS(4).e1_var         = '';
CONTROLLERS(4).e2_var         = '';
CONTROLLERS(4).kT_var         = 'kT';         % column vector 0:N-1
CONTROLLERS(4).has_ss_targets = true;
CONTROLLERS(4).xs_vars        = {'Xs_k_ot'};  % single 2xN matrix, already absolute
CONTROLLERS(4).us_vars        = {'Us_k_ot'};
CONTROLLERS(4).xs_absolute    = true;         % no offset needed
CONTROLLERS(4).us_absolute    = true;
CONTROLLERS(4).color          = [0.80 0.40 0.00];   % orange



%%  SECTION 2 -- LOAD DATA

nC = length(CONTROLLERS);

for c = 1:nC
    fprintf('Loading: %s  (%s) ...\n', CONTROLLERS(c).label, CONTROLLERS(c).file);
    raw = load(CONTROLLERS(c).file);

    % -- Output  Yk  (2 x N) ---------------------------------------------
    if ~isempty(CONTROLLERS(c).Yk_var) && isfield(raw, CONTROLLERS(c).Yk_var)
        % PATH A: variable is already a 2xN matrix
        CONTROLLERS(c).Yk = raw.(CONTROLLERS(c).Yk_var);
    else
        % PATH B: assemble from two separate row/column vector
        fprintf('%d\n',c);
        CONTROLLERS(c).Yk = [raw.(CONTROLLERS(c).t1_var)(:)'; ...
                              raw.(CONTROLLERS(c).t2_var)(:)'];
    end

    % -- Input  Uk  (2 x N) ----------------------------------------------
    if ~isempty(CONTROLLERS(c).Uk_var) && isfield(raw, CONTROLLERS(c).Uk_var)
        CONTROLLERS(c).Uk = raw.(CONTROLLERS(c).Uk_var);
    else
        CONTROLLERS(c).Uk = [raw.(CONTROLLERS(c).h1_var)(:)'; ...
                              raw.(CONTROLLERS(c).h2_var)(:)'];
    end

    % -- Reference  Rk_f  (2 x N) ----------------------------------------
    if ~isempty(CONTROLLERS(c).Rk_f_var) && isfield(raw, CONTROLLERS(c).Rk_f_var)
        CONTROLLERS(c).Rk_f = raw.(CONTROLLERS(c).Rk_f_var);
    else
        CONTROLLERS(c).Rk_f = [raw.(CONTROLLERS(c).R1_var)(:)'; ...
                                raw.(CONTROLLERS(c).R2_var)(:)'];
    end

    % -- Filtered error  ek_f  (2 x N) -----------------------------------
    if ~isempty(CONTROLLERS(c).ek_f_var) && isfield(raw, CONTROLLERS(c).ek_f_var)
        CONTROLLERS(c).ek_f = raw.(CONTROLLERS(c).ek_f_var);
    else
        CONTROLLERS(c).ek_f = [raw.(CONTROLLERS(c).e1_var)(:)'; ...
                                raw.(CONTROLLERS(c).e2_var)(:)'];
    end

    % -- Sample index  kT -------------------------------------------------
    if ~isempty(CONTROLLERS(c).kT_var) && isfield(raw, CONTROLLERS(c).kT_var)
        CONTROLLERS(c).kT = raw.(CONTROLLERS(c).kT_var)(:);  % force column
    else
        CONTROLLERS(c).kT = (0 : size(CONTROLLERS(c).Yk, 2) - 1)';
    end

    % -- Optional SS targets ----------------------------------------------
    if CONTROLLERS(c).has_ss_targets && ~isempty(CONTROLLERS(c).xs_vars)

        % xs: each entry is either a scalar vector or a pre-stacked matrix
        xs_rows = {};
        for vi = 1:length(CONTROLLERS(c).xs_vars)
            v = raw.(CONTROLLERS(c).xs_vars{vi});
            if min(size(v)) > 1
                xs_rows{end+1} = v;      % already 2xN (or nxN) matrix %#ok
            else
                xs_rows{end+1} = v(:)';  % scalar vector -> row %#ok
            end
        end
        xs_raw = vertcat(xs_rows{:});    % n_states x N

        % us: same logic
        us_rows = {};
        for vi = 1:length(CONTROLLERS(c).us_vars)
            v = raw.(CONTROLLERS(c).us_vars{vi});
            if min(size(v)) > 1
                us_rows{end+1} = v;      %#ok
            else
                us_rows{end+1} = v(:)';  %#ok
            end
        end
        us_raw = vertcat(us_rows{:});    % 2 x N

        N_ss = size(xs_raw, 2);

        % Apply operating-point offset if saved in perturbation form
        if CONTROLLERS(c).xs_absolute
            CONTROLLERS(c).Xs_k = xs_raw;
        else
            CONTROLLERS(c).Xs_k = xs_raw + repmat(TCL.Xs(:), 1, N_ss);
        end

        if CONTROLLERS(c).us_absolute
            CONTROLLERS(c).Us_k = us_raw;
        else
            CONTROLLERS(c).Us_k = us_raw + repmat(Us, 1, N_ss);
        end

    else
        CONTROLLERS(c).Xs_k = [];
        CONTROLLERS(c).Us_k = [];
    end

    CONTROLLERS(c).N = length(CONTROLLERS(c).kT);
    fprintf('   %d samples  (%.1f min)\n', CONTROLLERS(c).N, CONTROLLERS(c).N*Ts/60);
end

% Comparison window = shortest dataset
Ns = min([CONTROLLERS.N]);
kT = (0:Ns-1)';
fprintf('\nComparison window  Ns = %d samples  (%.1f min)\n\n', Ns, Ns*Ts/60);

% Reference signal for plotting (use first controller's reference)
Rk_f_ref = CONTROLLERS(1).Rk_f(:, 1:Ns);

%%  SECTION 3 -- PERFORMANCE INDICES

for c = 1:nC
    Yk_c = CONTROLLERS(c).Yk(:, 1:Ns);
    Uk_c = CONTROLLERS(c).Uk(:, 1:Ns);
    Rk_c = CONTROLLERS(c).Rk_f(:, 1:Ns);

    for i = 1:2
        CONTROLLERS(c).SSE(i)     = sum((Yk_c(i,:) - Rk_c(i,:)).^2);
        CONTROLLERS(c).SSMV(i)    = sum((Uk_c(i,:) - Us(i)).^2);
        CONTROLLERS(c).SSdelMV(i) = sum(diff(Uk_c(i,:)).^2);
    end
end

% -- Print formatted table ------------------------------------------------
col_w = 16;
sep   = repmat('=', 1, 20 + col_w*nC);
sep2  = repmat('-', 1, 20 + col_w*nC);
fprintf('%s\n', sep);
fprintf('  EXPERIMENT PERFORMANCE INDEX COMPARISON  (Ns=%d, %.1f min)\n', Ns, Ns*Ts/60);
fprintf('%s\n', sep);
hdr = sprintf('  %-18s', 'Index');
for c = 1:nC, hdr = [hdr sprintf('| %-14s', CONTROLLERS(c).label)]; end %#ok
fprintf('%s\n%s\n', hdr, sep2);

rows   = {'SSE (Loop 1)', 'SSE (Loop 2)', ...
          'SSMV (Loop 1)', 'SSMV (Loop 2)', ...
          'SSdeltaMV (L1)', 'SSdeltaMV (L2)'};
fields = {'SSE','SSE','SSMV','SSMV','SSdelMV','SSdelMV'};
loops  = [1 2 1 2 1 2];

for r = 1:6
    line = sprintf('  %-18s', rows{r});
    for c = 1:nC
        line = [line sprintf('| %14.4f ', CONTROLLERS(c).(fields{r})(loops(r)))]; %#ok
    end
    fprintf('%s\n', line);
end
fprintf('%s\n\n', sep);

%%  SECTION 4 -- PLOTS

line_styles = {'-','-','-','-','-'};   % cycles for up to 5 controllers

%% -- Figure 1: Yi(k) vs k ------------------------------------------------
figure('Name','Exp Comparison - Output Profiles','NumberTitle','off', ...
       'Position',[50 50 1000 700])

for i = 1:2
    subplot(2,1,i), hold on
    for c = 1:nC
        plot(kT, CONTROLLERS(c).Yk(i,1:Ns), ...
             'Color',     CONTROLLERS(c).color, ...
             'LineStyle', line_styles{min(c,5)}, ...
             'LineWidth', 1.8);
    end
    plot(kT, Rk_f_ref(i,:), 'k--', 'LineWidth', 1.2)
    xline(step_locs(1),'k:','Step 1','LabelVerticalAlignment','bottom','FontSize',8)
    xline(step_locs(2),'k:','Step 2','LabelVerticalAlignment','bottom','FontSize',8)
    grid on
    ylabel(sprintf('T_%d  (deg C)', i), 'FontSize', 11)
    if i == 1
        title('Output Y_i(k) and Reference R_i(k) - Experiment', 'FontSize', 12)
    end
    if i == 2, xlabel('Sample  k', 'FontSize', 11); end
    leg = cellfun(@(lb) sprintf('Y_%d^{%s}(k)', i, lb), {CONTROLLERS.label}, ...
                  'UniformOutput', false);
    leg{end+1} = sprintf('R_%d(k)', i);
    legend(leg{:}, 'Location','Best','FontSize',9)
end

%% -- Figure 2: Ui(k) vs k ------------------------------------------------
figure('Name','Exp Comparison - Heater Inputs','NumberTitle','off', ...
       'Position',[100 50 1000 700])

for i = 1:2
    subplot(2,1,i), hold on
    for c = 1:nC
        stairs(kT, CONTROLLERS(c).Uk(i,1:Ns)', ...
               'Color',     CONTROLLERS(c).color, ...
               'LineStyle', line_styles{min(c,5)}, ...
               'LineWidth', 1.8);
    end
    yline(5,  'k--', 'LineWidth', 0.9)
    yline(80, 'k--', 'LineWidth', 0.9)
    xline(step_locs(1),'k:')
    xline(step_locs(2),'k:')
    grid on
    ylim([0 105])
    ylabel(sprintf('U_%d(k)  (%% heater)', i), 'FontSize', 11)
    if i == 1, title('Heater Input U_i(k) - Experiment', 'FontSize', 12); end
    if i == 2, xlabel('Sample  k', 'FontSize', 11); end
    leg = cellfun(@(lb) sprintf('U_%d^{%s}(k)', i, lb), {CONTROLLERS.label}, ...
                  'UniformOutput', false);
    leg{end+1} = 'Safety bounds';
    legend(leg{:}, 'Location','Best','FontSize',9)
end

%% -- Figure 3: ef(k) vs k ------------------------------------------------
figure('Name','Exp Comparison - Filtered Error','NumberTitle','off', ...
       'Position',[150 50 1000 700])

for i = 1:2
    subplot(2,1,i), hold on
    for c = 1:nC
        stairs(kT, CONTROLLERS(c).ek_f(i,1:Ns)', ...
               'Color',     CONTROLLERS(c).color, ...
               'LineStyle', line_styles{min(c,5)}, ...
               'LineWidth', 1.8);
    end
    yline(0, 'k--', 'LineWidth', 0.8)
    xline(step_locs(1),'k:')
    xline(step_locs(2),'k:')
    grid on
    ylabel(sprintf('e_{f,%d}(k)  (deg C)', i), 'FontSize', 11)
    if i == 1, title('Filtered Error e_{f,i}(k) - Experiment', 'FontSize', 12); end
    if i == 2, xlabel('Sample  k', 'FontSize', 11); end
    legend({CONTROLLERS.label}, 'Location','Best','FontSize',9)
end

%% -- Figure 4: xs(k) and us(k) -- one row per controller with SS targets --
ss_idx = find([CONTROLLERS.has_ss_targets]);

if ~isempty(ss_idx)
    n_ss = length(ss_idx);
    figure('Name','Exp SS Targets xs(k) and us(k)','NumberTitle','off', ...
           'Position',[200 50 1000 300*n_ss])

    for row = 1:n_ss
        c = ss_idx(row);
        if isempty(CONTROLLERS(c).Xs_k), continue; end

        Ns_c   = min(Ns, CONTROLLERS(c).N);
        Xs_c   = CONTROLLERS(c).Xs_k(:, 1:Ns_c);
        Us_c   = CONTROLLERS(c).Us_k(:, 1:Ns_c);
        kT_c   = CONTROLLERS(c).kT(1:Ns_c);
        n_st_c = size(Xs_c, 1);

        subplot(n_ss, 2, 2*row-1), hold on
        for s = 1:n_st_c
            plot(kT_c, Xs_c(s,:), 'LineWidth', 1.5)
        end
        xline(step_locs(1),'k:')
        xline(step_locs(2),'k:')
        grid on
        ylabel('x_{s,i}(k)  (deg C)', 'FontSize', 10)
        title(sprintf('%s - SS Target States x_s(k)', CONTROLLERS(c).label), 'FontSize', 11)
        lx = arrayfun(@(s) sprintf('x_{s,%d}(k)', s), 1:n_st_c, 'UniformOutput', false);
        legend(lx{:}, 'Location','Best','FontSize',9)
        if row == n_ss, xlabel('Sample  k', 'FontSize', 10); end

        subplot(n_ss, 2, 2*row), hold on
        stairs(kT_c, Us_c(1,:)', 'LineWidth', 1.5)
        stairs(kT_c, Us_c(2,:)', 'LineWidth', 1.5)
        xline(step_locs(1),'k:')
        xline(step_locs(2),'k:')
        grid on
        ylabel('u_{s,i}(k)  (% heater)', 'FontSize', 10)
        title(sprintf('%s - SS Target Inputs u_s(k)', CONTROLLERS(c).label), 'FontSize', 11)
        legend('u_{s,1}(k)','u_{s,2}(k)', 'Location','Best','FontSize',9)
        if row == n_ss, xlabel('Sample  k', 'FontSize', 10); end
    end
end

%% -- Figure 5: Performance index bar chart --------------------------------
figure('Name','Exp Performance Index Comparison','NumberTitle','off', ...
       'Position',[250 50 1000 500])

data_loop1 = zeros(3, nC);
data_loop2 = zeros(3, nC);
for c = 1:nC
    data_loop1(:,c) = [CONTROLLERS(c).SSE(1);
                       CONTROLLERS(c).SSMV(1);
                       CONTROLLERS(c).SSdelMV(1)];
    data_loop2(:,c) = [CONTROLLERS(c).SSE(2);
                       CONTROLLERS(c).SSMV(2);
                       CONTROLLERS(c).SSdelMV(2)];
end

labels = {'SSE','SSMV','SS\DeltaMV'};
group  = categorical(labels);
group  = reordercats(group, labels);
colors = reshape([CONTROLLERS.color], 3, nC)';

for sp = 1:2
    subplot(1,2,sp)
    data   = data_loop1;
    ltitle = 'Loop 1  (T_1 / H_1)';
    if sp == 2
        data   = data_loop2;
        ltitle = 'Loop 2  (T_2 / H_2)';
    end
    b = bar(group, data', 'grouped');
    for c = 1:nC, b(c).FaceColor = colors(c,:); end
    grid on
    ylabel('Index Value', 'FontSize', 10)
    title(ltitle, 'FontSize', 11)
    legend({CONTROLLERS.label}, 'Location','Best','FontSize',9)
end
sgtitle('Performance Indices - Experiment Comparison', 'FontSize', 13)


%%  SECTION 5 -- SAVE

exp_compare.Ns = Ns;
exp_compare.kT = kT;
exp_compare.Ts = Ts;

for c = 1:nC
    fn = matlab.lang.makeValidName(CONTROLLERS(c).label);
    exp_compare.(fn).label   = CONTROLLERS(c).label;
    exp_compare.(fn).Yk      = CONTROLLERS(c).Yk(:, 1:Ns);
    exp_compare.(fn).Uk      = CONTROLLERS(c).Uk(:, 1:Ns);
    exp_compare.(fn).Rk_f    = CONTROLLERS(c).Rk_f(:, 1:Ns);
    exp_compare.(fn).ek_f    = CONTROLLERS(c).ek_f(:, 1:Ns);
    exp_compare.(fn).SSE     = CONTROLLERS(c).SSE;
    exp_compare.(fn).SSMV    = CONTROLLERS(c).SSMV;
    exp_compare.(fn).SSdelMV = CONTROLLERS(c).SSdelMV;
    if CONTROLLERS(c).has_ss_targets && ~isempty(CONTROLLERS(c).Xs_k)
        Ns_c = min(Ns, CONTROLLERS(c).N);
        exp_compare.(fn).Xs_k = CONTROLLERS(c).Xs_k(:, 1:Ns_c);
        exp_compare.(fn).Us_k = CONTROLLERS(c).Us_k(:, 1:Ns_c);
    end
end

save TCL_Compare_ExpResults  exp_compare
fprintf('Results saved to  TCL_Compare_ExpResults.mat\n');
fprintf('All done.\n');
