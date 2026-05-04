function [pop, GD, SP] = dqn_main(vel_profile)
%DQN_MAIN  Deep Q-Network (DQN) train controller — Pareto approximation.
%
% Trains a DQN agent to control the train step-by-step in a position-step
% environment.  To produce a Pareto front (time vs energy), the agent is
% trained with a scalarised reward under multiple weight values:
%       r_step = -(w * Δt + (1-w) * ΔE)
%
% Architecture: 4→64→64→5 MLP with ReLU hidden layers.
% Training: epsilon-greedy DQN with experience replay buffer.
%
% Spec from E3:
%   Episodes     = 500 per weight
%   epsilon      : 1.0 → 0.1 (linear decay)
%   Replay buffer: 10 000
%   Mini-batch   : 64
%   Discount γ   : 0.99
%   Target net   : hard update every 50 episodes
%
% Actions (5 discrete):
%   1 = full brake   (−a_brake_max)
%   2 = half brake   (−a_brake_max/2)
%   3 = coast        (only gravity + resistance)
%   4 = half traction(+a_trac_max/2)
%   5 = full traction(+a_trac_max)
%
% Requires globals: vel_profile, gradient, inertial_mass, Davis,
%   Max_tractive_power, max_accel_trac, max_accel_brake, dimension, pop_size

    global vel_profile gradient dimension
    global inertial_mass Davis Max_tractive_power max_accel_trac max_accel_brake
    global pop_size iterations

    % ---- environment constants ----
    DS          = 10;        % position step (m)
    BRAKE_REGEN = 0.6;
    N_ACTIONS   = 5;
    PENALTY     = 1e3;       % large penalty for time violation

    % ---- DQN hyper-parameters ----
    N_EPS    = 500;          % episodes per weight
    BUF_MAX  = 10000;        % replay buffer size
    BATCH    = 64;           % mini-batch size
    GAMMA    = 0.99;         % discount factor
    LR       = 1e-3;         % learning rate
    TGT_UPD  = 50;           % target network update interval (episodes)
    EPS_START = 1.0;
    EPS_END   = 0.1;

    % ---- route discretisation ----
    S_max     = max(vel_profile(:,1)) * 1000;   % m
    K         = floor(S_max / DS);
    s_grid    = (0:K) * DS;

    vp_dist   = vel_profile(:,1) * 1000;
    vp_speed  = vel_profile(:,2) / 3.6;         % m/s
    vlim_at   = @(s) interp1(vp_dist, vp_speed, s, 'previous', vp_speed(end));

    if ~isempty(gradient) && size(gradient,1) > 1
        gr_dist = gradient(:,1) * 1000;
        gr_val  = gradient(:,2);
        grad_at = @(s) interp1(gr_dist, gr_val*9.81/1000, s, 'linear', 0);
    else
        grad_at = @(s) 0;
    end

    % ---- weight sweep ----
    w_vals = [0, 0.05, 0.1, 0.2, 0.3, 0.5, 0.7, 0.8, 0.9, 0.95, 1.0];
    n_w    = length(w_vals);

    T_out = zeros(n_w, 1);
    E_out = zeros(n_w, 1);

    fprintf('[DQN] Starting %d weight sweeps × %d episodes each\n', n_w, N_EPS);

    for wi = 1:n_w
        w = w_vals(wi);

        % ---- initialise network weights ----
        net  = init_net(4, 64, N_ACTIONS);   % online network
        net_tgt = net;                        % target network

        % ---- replay buffer ----
        buf.s  = zeros(BUF_MAX, 4);
        buf.a  = zeros(BUF_MAX, 1, 'int32');
        buf.r  = zeros(BUF_MAX, 1);
        buf.s2 = zeros(BUF_MAX, 4);
        buf.d  = zeros(BUF_MAX, 1);   % done flag
        buf.n  = 0;   % current fill
        buf.i  = 0;   % write pointer

        best_T = Inf; best_E = Inf;
        eps_now = EPS_START;

        for ep = 1:N_EPS
            eps_now = EPS_START + (EPS_END - EPS_START) * (ep-1) / (N_EPS-1);

            % ---- reset episode ----
            k    = 0;        % position index (0..K)
            v    = 0;        % speed (m/s)
            T_ep = 0;
            E_ep = 0;
            done = false;

            % Estimate T_target for state normalisation (rough)
            T_ref = S_max / (max(vp_speed) * 0.7 + eps);

            while ~done
                s_norm = k / K;
                v_norm = v / (max(vp_speed) + eps);
                rd_norm = (K - k) / K;
                rt_norm = max(0, (T_ref - T_ep)) / T_ref;
                state = [s_norm, v_norm, rd_norm, rt_norm];

                % --- epsilon-greedy action ---
                if rand < eps_now
                    a = randi(N_ACTIONS);
                else
                    Q = forward_net(net, state);
                    [~, a] = max(Q);
                end

                % --- environment step ---
                pos_cur  = k * DS;
                vl       = vlim_at(pos_cur);
                g_acc    = grad_at(pos_cur);
                resist   = (Davis(1) + Davis(2)*v + Davis(3)*v^2) / inertial_mass;

                switch a
                    case 1, a_applied = -(max_accel_brake + g_acc - resist);
                    case 2, a_applied = -(max_accel_brake/2 + g_acc - resist);
                    case 3, a_applied = g_acc - resist;        % coast
                    case 4, a_applied = max_accel_trac/2 - resist + g_acc;
                    case 5, a_applied = max_accel_trac   - resist + g_acc;
                end

                v_next_sq = v^2 + 2 * a_applied * DS;
                v_next    = sqrt(max(0, v_next_sq));
                v_next    = min(v_next, vl);

                v_avg = (v + v_next) / 2;
                dT = DS / max(v_avg, 0.5);

                net_a = (v_next^2 - v^2) / (2*DS) + resist - g_acc;
                if net_a > 0
                    P_t = min(net_a * inertial_mass * 1000 * v_avg, Max_tractive_power);
                    dE  = P_t * dT / 3.6e6;
                else
                    dE = -BRAKE_REGEN * abs(net_a) * inertial_mass * 1000 * DS / 3.6e6;
                    dE = max(dE, -1e-6);
                end

                T_ep = T_ep + dT;
                E_ep = E_ep + dE;
                k    = k + 1;
                v    = v_next;

                % --- done? ---
                done  = (k >= K);
                % small penalty if speed at terminal not near zero
                r_term = 0;
                if done
                    r_term = -10 * (v / (max(vp_speed)+eps));
                end
                reward = -(w * dT + (1-w) * max(0,dE)) + r_term;

                % ---- store transition ----
                s_norm2 = k / K;
                v_norm2 = v / (max(vp_speed)+eps);
                rd_norm2 = (K - k) / K;
                rt_norm2 = max(0,(T_ref-T_ep)) / T_ref;
                next_state = [s_norm2, v_norm2, rd_norm2, rt_norm2];

                buf = store_buf(buf, state, a, reward, next_state, done, BUF_MAX);

                % ---- training step ----
                if buf.n >= BATCH
                    [s_b, a_b, r_b, s2_b, d_b] = sample_buf(buf, BATCH);
                    net = train_step(net, net_tgt, s_b, a_b, r_b, s2_b, d_b, GAMMA, LR);
                end
            end

            % ---- target network update ----
            if mod(ep, TGT_UPD) == 0
                net_tgt = net;
            end

            % ---- track best feasible solution ----
            if v < 2.0  % arrived with low speed (approximately at rest)
                if E_ep < best_E
                    best_T = T_ep; best_E = max(0, E_ep);
                end
            end

            if mod(ep,100)==0 || ep==1
                fprintf('[DQN] w=%.2f | ep %3d/%d | eps=%.3f | T=%.1fs E=%.4f kWh\n', ...
                    w, ep, N_EPS, eps_now, T_ep, max(0,E_ep));
            end
        end

        if isinf(best_T)
            % Fallback: use last episode result
            best_T = T_ep; best_E = max(0, E_ep);
        end
        T_out(wi) = best_T;
        E_out(wi) = best_E;
        fprintf('[DQN] w=%.2f done | T=%.1f s | E=%.4f kWh\n', w, best_T, best_E);
    end

    % ---- assemble population matrix ----
    valid = ~isnan(T_out) & ~isnan(E_out);
    T_v = T_out(valid); E_v = E_out(valid);
    n_v = sum(valid);

    if n_v == 0
        warning('dqn_main: no valid solutions produced.');
        pop = zeros(0, dimension+4);
        GD = NaN; SP = NaN;
        return;
    end

    F_mat   = [T_v, E_v];
    rank_v  = pareto_rank(F_mat);
    crowd_v = pareto_crowd(F_mat);
    pop = [zeros(n_v, dimension), T_v, E_v, rank_v, crowd_v];

    GD = NaN(1, max(1,iterations));
    SP = NaN(1, max(1,iterations));
    fprintf('[DQN] done | %d solutions on Pareto approximation\n', n_v);
end

% =========================================================================
%  NEURAL NETWORK HELPERS  (no toolbox required)
% =========================================================================

function net = init_net(n_in, n_hid, n_out)
    s = 0.1;
    net.W1 = s * randn(n_hid, n_in);  net.b1 = zeros(n_hid, 1);
    net.W2 = s * randn(n_hid, n_hid); net.b2 = zeros(n_hid, 1);
    net.W3 = s * randn(n_out, n_hid); net.b3 = zeros(n_out, 1);
    % Adam moments
    net.mW1=0; net.vW1=0; net.mb1=0; net.vb1=0;
    net.mW2=0; net.vW2=0; net.mb2=0; net.vb2=0;
    net.mW3=0; net.vW3=0; net.mb3=0; net.vb3=0;
    net.t=0;  % Adam step counter
end

function Q = forward_net(net, x)
    x  = x(:);
    h1 = relu(net.W1 * x + net.b1);
    h2 = relu(net.W2 * h1 + net.b2);
    Q  = net.W3 * h2 + net.b3;   % linear output (N_ACTIONS × 1)
    Q  = Q(:)';                   % row vector
end

function [h1, h2, Q] = forward_full(net, x)
    h1 = relu(net.W1 * x + net.b1);
    h2 = relu(net.W2 * h1 + net.b2);
    Q  = net.W3 * h2 + net.b3;
end

function y = relu(x)
    y = max(0, x);
end

function net = train_step(net, net_tgt, s_b, a_b, r_b, s2_b, d_b, gamma, lr)
% One gradient step on a mini-batch of transitions.
    B = size(s_b, 1);
    beta1 = 0.9; beta2 = 0.999; adam_eps = 1e-8;
    net.t = net.t + 1;

    % Compute targets using target network
    targets = zeros(B, 1);
    for i = 1:B
        Q_next = forward_net(net_tgt, s2_b(i,:));
        if d_b(i)
            targets(i) = r_b(i);
        else
            targets(i) = r_b(i) + gamma * max(Q_next);
        end
    end

    % Gradient accumulation
    gW1 = zeros(size(net.W1)); gb1 = zeros(size(net.b1));
    gW2 = zeros(size(net.W2)); gb2 = zeros(size(net.b2));
    gW3 = zeros(size(net.W3)); gb3 = zeros(size(net.b3));

    for i = 1:B
        x = s_b(i,:)';
        [h1, h2, Q] = forward_full(net, x);
        a = a_b(i);
        delta = Q(a) - targets(i);      % scalar error

        % Output layer gradient
        dL_dQ3 = zeros(size(Q)); dL_dQ3(a) = 2*delta/B;
        gW3 = gW3 + dL_dQ3 * h2';
        gb3 = gb3 + dL_dQ3;

        % Hidden layer 2
        dL_dh2 = net.W3' * dL_dQ3;
        dL_dh2 = dL_dh2 .* (h2 > 0);  % ReLU derivative
        gW2 = gW2 + dL_dh2 * h1';
        gb2 = gb2 + dL_dh2;

        % Hidden layer 1
        dL_dh1 = net.W2' * dL_dh2;
        dL_dh1 = dL_dh1 .* (h1 > 0);
        gW1 = gW1 + dL_dh1 * x';
        gb1 = gb1 + dL_dh1;
    end

    % Adam update helper (inline)
    t = net.t;
    [net.W1, net.mW1, net.vW1] = adam_upd(net.W1, gW1, net.mW1, net.vW1, lr, beta1, beta2, adam_eps, t);
    [net.b1, net.mb1, net.vb1] = adam_upd(net.b1, gb1, net.mb1, net.vb1, lr, beta1, beta2, adam_eps, t);
    [net.W2, net.mW2, net.vW2] = adam_upd(net.W2, gW2, net.mW2, net.vW2, lr, beta1, beta2, adam_eps, t);
    [net.b2, net.mb2, net.vb2] = adam_upd(net.b2, gb2, net.mb2, net.vb2, lr, beta1, beta2, adam_eps, t);
    [net.W3, net.mW3, net.vW3] = adam_upd(net.W3, gW3, net.mW3, net.vW3, lr, beta1, beta2, adam_eps, t);
    [net.b3, net.mb3, net.vb3] = adam_upd(net.b3, gb3, net.mb3, net.vb3, lr, beta1, beta2, adam_eps, t);
end

function [p, m, v] = adam_upd(p, g, m, v, lr, b1, b2, ep, t)
    m = b1*m + (1-b1)*g;
    v = b2*v + (1-b2)*g.^2;
    mh = m / (1-b1^t);
    vh = v / (1-b2^t);
    p = p - lr * mh ./ (sqrt(vh) + ep);
end

% ---- replay buffer ----
function buf = store_buf(buf, s, a, r, s2, d, buf_max)
    buf.i = mod(buf.i, buf_max) + 1;
    buf.s(buf.i,:)  = s;
    buf.a(buf.i)    = a;
    buf.r(buf.i)    = r;
    buf.s2(buf.i,:) = s2;
    buf.d(buf.i)    = d;
    buf.n = min(buf.n + 1, buf_max);
end

function [s_b, a_b, r_b, s2_b, d_b] = sample_buf(buf, batch)
    idx  = randperm(buf.n, min(batch, buf.n));
    s_b  = buf.s(idx,:);
    a_b  = buf.a(idx);
    r_b  = buf.r(idx);
    s2_b = buf.s2(idx,:);
    d_b  = buf.d(idx);
end

% ---- Pareto rank / crowd ----
function rank = pareto_rank(F)
    n = size(F,1); rank = ones(n,1);
    for i = 1:n
        for j = 1:n
            if i~=j && all(F(j,:)<=F(i,:)) && any(F(j,:)<F(i,:))
                rank(i) = rank(i)+1;
            end
        end
    end
end

function cd = pareto_crowd(F)
    n = size(F,1);
    if n<=2, cd=Inf(n,1); return; end
    cd = zeros(n,1);
    for m = 1:size(F,2)
        [~,idx]=sort(F(:,m)); cd(idx(1))=Inf; cd(idx(end))=Inf;
        fr=F(idx(end),m)-F(idx(1),m)+eps;
        for i=2:n-1
            cd(idx(i))=cd(idx(i))+(F(idx(i+1),m)-F(idx(i-1),m))/fr;
        end
    end
end
