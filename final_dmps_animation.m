%  Animation της προσομοίωσης DMP (3 κύκλοι)
% Αυτοτελές script: εκπαιδεύει τα DMP, παράγει την τροχιά 3 κύκλων (με τα
% δ/θδ) και κάνει animation, δείχνοντας το Τμήμα Β μετατοπισμένο/στραμμένο
% σε κάθε κύκλο. 
% =========================================================================
clear; clc; close all;

% Επιλογή: true => DMP και για τον προσανατολισμό (Ερώτημα 3)
%          false => DMP μόνο για τη θέση     (Ερώτημα 2)
USE_ORIENTATION_DMP = true;

ANIM_STEP = 10;        % βήμα δειγμάτων ανά frame (μεγαλύτερο = πιο γρήγορο)

% Ρομπότ 
[scara, q0] = create_robot();

%  Παράμετροι & σημεία 
dt = 0.001;
P0     = [0.707, -0.293, 1.85, pi/4];
P1     = [0.19,  -1.00,  1.63,  0   ];   % σύλληψη @t=2.3s: x = -0.5 + 0.3*2.3
P2     = [0.34,  -1.00,  1.63,  0   ];   % τέλος σύλληψης @t=2.8s
P_VIA1 = [0.7,   -1.00,  1.85,  0   ];
P3     = [0.08,   1.00,  1.85,  0   ];
P4     = [0.08,   1.00,  1.15,  pi/2];

alpha_z = 25; beta_z = alpha_z/4; alpha_x = 1.0; n_bfs = 40;
[c, h_bfs] = init_bfs(n_bfs);

%  Εκπαίδευση DMP Θέσης (Τμήμα 4)
tau_A = 3.4; t_A = (0:dt:tau_A)';
w_dmp_A = zeros(n_bfs, 2);
for j = 1:2
    Y = quintic_full(P_VIA1(j), P3(j), 0, 0, 0, 0, 0, tau_A, t_A);
    dY = [0; diff(Y)/dt]; ddY = [0; diff(dY)/dt];
    w_dmp_A(:,j) = train_dmp(Y, dY, ddY, t_A, tau_A, alpha_z, beta_z, alpha_x, c, h_bfs);
end

%  Εκπαίδευση DMP Προσανατολισμού (Τμήμα 5β) 
tau_B = 0.95; t_B = (0:dt:tau_B)';
Phi = quintic_full(P3(4), P4(4), 0, 0, 0, 0, 0, tau_B, t_B);
dPhi = [0; diff(Phi)/dt]; ddPhi = [0; diff(dPhi)/dt];
w_dmp_Phi = train_dmp(Phi, dPhi, ddPhi, t_B, tau_B, alpha_z, beta_z, alpha_x, c, h_bfs);

%  Σενάρια (όπως η εκφώνηση) 
deltas = [ 0.0,  0.0,  0.0;
          -0.2,  0.2,  0.0;
           0.4, -0.15, 0.0 ];
theta_deltas = [0, 20*(pi/180), -25*(pi/180)];

if USE_ORIENTATION_DMP
    w_phi_use = w_dmp_Phi;  th_use = theta_deltas;   % Ερώτημα 3
else
    w_phi_use = [];         th_use = zeros(1,3);     % Ερώτημα 2
end

%  Τροχιά 3 κύκλων 
[cart_traj, q_traj] = simulate_3_cycles(deltas, th_use, P0, P1, P2, P_VIA1, ...
    P3, P4, dt, tau_A, tau_B, w_dmp_A, w_phi_use, alpha_z, beta_z, alpha_x, c, h_bfs);
N = size(cart_traj,1);
t_traj = (0:N-1)' * dt;

% =========================================================================
% ANIMATION (3 κύκλοι)
% =========================================================================
figure('Name','DMP Assembly Animation (3 Cycles)','Position',[100,100,950,720]);
view(3); grid on; hold on;
axis([-1.5 3.5 -2.5 2.5 0 3]);
xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)');

% Ταινιόδρομοι (επιφάνειες)
[xa,ya] = meshgrid(-1.5:0.5:3.0, -1.5:0.1:-0.5);
surf(xa,ya,1.5*ones(size(xa)),'FaceColor',[0.7 0.8 0.9],'EdgeAlpha',0.1);
[xb,yb] = meshgrid(-1.5:0.5:3.0, 0.5:0.1:1.5);
surf(xb,yb,1.0*ones(size(xb)),'FaceColor',[0.9 0.8 0.7],'EdgeAlpha',0.1);

h_partA = []; h_partB = []; h_gripper = []; h_target = [];

for i = 1:ANIM_STEP:N
    t = t_traj(i);
    cyc = min(floor(t/12) + 1, 3);          % τρέχων κύκλος 1..3
    t_cycle = mod(t, 12.0);
    dx = deltas(cyc,1); dy = deltas(cyc,2); thd = th_use(cyc);

    T_tcp = scara.fkine(q_traj(i,:));

    %  Τμήμα Β (πράσινο): κινείται στον ταινιόδρομο Β και κάθεται στη
    %     μετατοπισμένη/στραμμένη θέση (δ, θδ) που εντόπισε η κάμερα 
    if t_cycle <= 2.4
        x_B = (2.0 + dx) - 0.8*t_cycle;     % -0.8 m/s για 2.4s
    else
        x_B = 0.08 + dx;
    end
    posB = [x_B, 1.0 + dy, 1.0];
    rotB = t2r(trotz(thd));
    if ~isempty(h_partB), delete(h_partB); end
    h_partB = draw_rotated_box(posB, [0.3,0.3,0.1], [0.2 0.8 0.2], rotB);

    % Στόχος συναρμολόγησης (δείκτης)
    if ~isempty(h_target), delete(h_target); end
    h_target = plot3(0.08+dx, 1.0+dy, 1.1, 'kx', 'MarkerSize',10, 'LineWidth',1.5);

    % --- Τμήμα Α (κόκκινο): ταινιόδρομος Α -> μεταφορά -> τοποθέτηση ---
    if t_cycle <= 2.8
        posA = [-0.5 + 0.3*t_cycle, -1, 1.5]; rotA = eye(3);
    elseif t_cycle <= 8.75
        posA = T_tcp.t' - [0,0,0.13]; rotA = t2r(trotz(cart_traj(i,4)));
    else
        posA = [0.08+dx, 1.0+dy, 1.02]; rotA = t2r(trotz(pi/2 + thd));
    end
    if ~isempty(h_partA), delete(h_partA); end
    h_partA = draw_rotated_box(posA, [0.25,0.25,0.18], [0.8 0.2 0.2], rotA);

    %  Αρπάγη 
    if ~isempty(h_gripper), delete(h_gripper); end
    s_grip = gripper_open_frac(t_cycle);
    h_gripper = draw_gripper(T_tcp.T, s_grip);

    title(sprintf('DMP — Κύκλος %d/3   |   \\delta=[%.2f, %.2f] m,  \\theta_\\delta=%.0f°   |   t=%.1fs', ...
        cyc, dx, dy, thd*180/pi, t));

    scara.plot(q_traj(i,:), 'workspace',[-1.5 3.5 -2.5 2.5 0 3], 'delay',0);
    drawnow;
end

% =========================================================================
% ΤΟΠΙΚΕΣ ΣΥΝΑΡΤΗΣΕΙΣ
% =========================================================================
function [cart_traj, q_traj, sim_X_all, sim_Y_all, sim_Phi_all] = simulate_3_cycles(...
    deltas, theta_deltas, P0, P1, P2, P_VIA1, P3, P4, dt, ...
    tau_A, tau_B, w_dmp_A, w_dmp_Phi, alpha_z, beta_z, alpha_x, c, h_bfs)

    cart_traj = [];
    sim_X_all   = zeros(round(tau_A/dt)+1, 3);
    sim_Y_all   = zeros(round(tau_A/dt)+1, 3);
    sim_Phi_all = zeros(round(tau_B/dt)+1, 3);

    for cycle = 1:3
        g_P3    = P3; g_P3(1:3) = P3(1:3) + deltas(cycle,:);
        g_P4    = P4; g_P4(1:3) = P4(1:3) + deltas(cycle,:);
        g_P4(4) = P4(4) + theta_deltas(cycle);

        %  Τμήμα 1: P0 → P1 (0-2.3s)
        t1 = (0:dt:2.3)'; n1 = length(t1);
        tr1 = zeros(n1, 4);
        tr1(:,1) = quintic_full(P0(1), P1(1), 0, 0.3, 0, 0, 0, 2.3, t1);
        tr1(:,2) = quintic_full(P0(2), P1(2), 0, 0.0, 0, 0, 0, 2.3, t1);
        tr1(:,3) = P0(3)*ones(n1,1);
        idx_z = t1 >= 1.3;
        tr1(idx_z,3) = quintic_full(P0(3), P1(3), 0, 0, 0, 0, 1.3, 2.3, t1(idx_z));
        tr1(:,4) = quintic_full(P0(4), P1(4), 0, 0.0, 0, 0, 0, 2.3, t1);

        %  Τμήμα 2: P1 → P2 (2.3-2.8s) σταθερή v=0.3
        t2 = (2.3:dt:2.8)'; n2 = length(t2);
        tr2 = [P1(1)+0.3*(t2-2.3), repmat(P1(2:4), n2, 1)];

        %  Τμήμα 3: P2 → P_VIA1 (2.8→3.6s)
        t3 = (2.8:dt:3.6)'; n3 = length(t3);
        tr3 = zeros(n3, 4);
        tr3(:,1) = quintic_full(P2(1), P_VIA1(1), 0.3, 0, 0, 0, 2.8, 3.6, t3);
        tr3(:,2) = P2(2)*ones(n3,1);
        tr3(:,3) = quintic_full(P2(3), P_VIA1(3), 0, 0, 0, 0, 2.8, 3.6, t3);
        tr3(:,4) = P2(4)*ones(n3,1);

        %  Τμήμα 4: DMP Θέσης P_VIA1 → g_P3 (3.6→7.0s, tau_A=3.4)
        [sim_X, ~, ~] = rollout_dmp(P_VIA1(1), g_P3(1), w_dmp_A(:,1), tau_A, dt, alpha_z, beta_z, alpha_x, c, h_bfs);
        [sim_Y, ~, ~] = rollout_dmp(P_VIA1(2), g_P3(2), w_dmp_A(:,2), tau_A, dt, alpha_z, beta_z, alpha_x, c, h_bfs);
        sim_X = sim_X + (g_P3(1)-sim_X(end))*linspace(0,1,numel(sim_X))';
        sim_Y = sim_Y + (g_P3(2)-sim_Y(end))*linspace(0,1,numel(sim_Y))';
        sim_X_all(:,cycle) = sim_X;
        sim_Y_all(:,cycle) = sim_Y;
        n4 = length(sim_X);
        tr4 = [sim_X, sim_Y, P_VIA1(3)*ones(n4,1), P_VIA1(4)*ones(n4,1)];

        %  Τμήμα 5: Κατάβαση + Προσανατολισμός (7.0→8.75s, σε δύο φάσεις)
        z_rim = 1.21;

        %  Τμήμα 5α: ευθεία κάθοδος ως το χείλος (7.0→7.8s)
        t5a = (7.0:dt:7.8)'; n5a = length(t5a);
        tr5a_Z = quintic_full(g_P3(3), z_rim, 0, 0, 0, 0, 7.0, 7.8, t5a);
        tr5a = [g_P3(1)*ones(n5a,1), g_P3(2)*ones(n5a,1), tr5a_Z, g_P3(4)*ones(n5a,1)];

        %  Τμήμα 5β: βίδωμα (7.8→8.75s) — τελική κάθοδος + στροφή
        t5b = (7.8:dt:8.75)'; n5b = length(t5b);
        tr5b_Z = quintic_full(z_rim, g_P4(3), 0, 0, 0, 0, 7.8, 8.75, t5b);
        if isempty(w_dmp_Phi)
            sim_Phi = quintic_full(P3(4), g_P4(4), 0, 0, 0, 0, 7.8, 8.75, t5b);   % Ερώτημα 2
        else
            [sim_Phi, ~, ~] = rollout_dmp(g_P3(4), g_P4(4), w_dmp_Phi, tau_B, dt, alpha_z, beta_z, alpha_x, c, h_bfs);
            sim_Phi = sim_Phi + (g_P4(4)-sim_Phi(end))*linspace(0,1,numel(sim_Phi))';   % Ερώτημα 3
        end
        sim_Phi_all(:,cycle) = sim_Phi;
        tr5b = [g_P4(1)*ones(n5b,1), g_P4(2)*ones(n5b,1), tr5b_Z, sim_Phi];

        %  Τμήμα 6: Pause (8.75-9.25s)
        t6 = (8.75:dt:9.25)'; n6 = length(t6);
        tr6 = repmat(g_P4, n6, 1);

        %  Τμήμα 7: Επιστροφή με fly-by (9.25→12s)
        P_RETRACT = [g_P4(1), g_P4(2), 1.85, g_P4(4)];
        P_VIA_RET = [1.00, 0.50, 1.85, 3*pi/8];
        t_7b = 11.15 - 10.05; t_7c = 12.0 - 11.15;
        v_via = zeros(1,4);
        for j = 1:4
            V_prev = (P_VIA_RET(j) - P_RETRACT(j)) / t_7b;
            V_next = (P0(j)        - P_VIA_RET(j)) / t_7c;
            if sign(V_prev) == sign(V_next)
                v_via(j) = sign(V_prev) * sqrt(abs(V_prev * V_next));
            else
                v_via(j) = 0;
            end
        end
        t7a = (9.25:dt:10.05)';  tr7a = zeros(length(t7a),4);
        t7b = (10.05:dt:11.15)'; tr7b = zeros(length(t7b),4);
        t7c = (11.15:dt:12.0)';  tr7c = zeros(length(t7c),4);
        for j = 1:4
            tr7a(:,j) = quintic_full(g_P4(j),     P_RETRACT(j), 0,        0,        0, 0,  9.25, 10.05, t7a);
            tr7b(:,j) = quintic_full(P_RETRACT(j), P_VIA_RET(j), 0,        v_via(j), 0, 0, 10.05, 11.15, t7b);
            tr7c(:,j) = quintic_full(P_VIA_RET(j), P0(j),        v_via(j), 0,        0, 0, 11.15, 12.0,  t7c);
        end

        cycle_traj = [tr1(1:end-1,:); tr2(1:end-1,:); tr3(1:end-1,:);
                      tr4(1:end-1,:); tr5a(1:end-1,:); tr5b(1:end-1,:); tr6(1:end-1,:);
                      tr7a(1:end-1,:); tr7b(1:end-1,:); tr7c];

        cart_traj = [cart_traj; cycle_traj(1:end-1,:)];
    end
    cart_traj = [cart_traj; P0];

    %  Αντίστροφη Κινηματική (L1=L2=1)
    N = size(cart_traj,1); q_traj = zeros(N,4);
    for i = 1:N
        x = cart_traj(i,1); y = cart_traj(i,2);
        z = cart_traj(i,3); phi = cart_traj(i,4);
        D  = (x^2 + y^2 - 2) / 2;
        q2 = atan2(sqrt(max(0,1-D^2)), D);
        q1 = atan2(y,x) - atan2(sin(q2), 1+cos(q2));
        q_traj(i,:) = [q1, q2, 1.95-z, q1+q2-phi];
    end
end

function w = train_dmp(y, yd, ydd, t, tau, alpha_z, beta_z, alpha_x, c, h)
    x = exp(-alpha_x * (t / tau));
    y0 = y(1); g = y(end);
    f_target = (tau^2)*ydd - alpha_z*(beta_z*(g - y) - tau*yd);
    scale = (g - y0); if abs(scale) < 1e-4, scale = 1; end
    w = zeros(length(c), 1);
    for i = 1:length(c)
        psi = exp(-h(i) * (x - c(i)).^2);
        w(i) = sum(x .* psi .* f_target) / (sum((x.^2) .* psi) * scale + 1e-10);
    end
end

function [Y, dY, ddY] = rollout_dmp(y0, g, w, tau, dt, alpha_z, beta_z, alpha_x, c, h)
    N = round(tau/dt) + 1;
    Y = zeros(N,1); dY = zeros(N,1); ddY = zeros(N,1);
    y = y0; z = 0; x = 1;
    scale = (g - y0); if abs(scale) < 1e-4, scale = 1; end
    for i = 1:N
        Y(i) = y; dY(i) = z/tau;
        psi = exp(-h .* (x - c).^2);
        f = (sum(w .* psi) / (sum(psi) + 1e-10)) * x * scale;
        dz = alpha_z*(beta_z*(g - y) - z) + f;
        ddY(i) = dz/tau^2;
        y = y + z*(dt/tau);
        z = z + dz*(dt/tau);
        x = x + (-alpha_x*x)*(dt/tau);
    end
end

function [c, h] = init_bfs(n_bfs)
    c = exp(-1.0 * linspace(0, 1, n_bfs))';
    h = 1 ./ (diff(c).^2); h = [h; h(end)];
end

function q_out = quintic_full(q0, qf, v0, vf, a0, af, t0, tf, t_arr)
    H = tf - t0;
    if H <= 0, q_out = q0*ones(length(t_arr),1); return; end
    k0 = q0; k1 = v0; k2 = a0/2;
    A = qf - q0 - v0*H - (a0/2)*H^2;
    B = vf - v0 - a0*H;
    C = af - a0;
    k3 = (10*A - 4*B*H + 0.5*C*H^2) / H^3;
    k4 = (-15*A + 7*B*H - C*H^2)    / H^4;
    k5 = (6*A  - 3*B*H + 0.5*C*H^2) / H^5;
    q_out = zeros(length(t_arr),1);
    for k = 1:length(t_arr)
        tau = max(0, min(t_arr(k)-t0, H));
        q_out(k) = k0 + k1*tau + k2*tau^2 + k3*tau^3 + k4*tau^4 + k5*tau^5;
    end
end

function h = draw_rotated_box(base_center, dims, color, R)
    dx=dims(1)/2; dy=dims(2)/2; dz=dims(3);
    v=[-dx -dy 0; dx -dy 0; dx dy 0; -dx dy 0;
       -dx -dy dz; dx -dy dz; dx dy dz; -dx dy dz]';
    v = R*v + base_center(:); v = v';
    f = [1 2 6 5; 2 3 7 6; 3 4 8 7; 4 1 5 8; 1 2 3 4; 5 6 7 8];
    h = patch('Vertices',v,'Faces',f,'FaceColor',color,'FaceAlpha',0.9);
end

function h = draw_gripper(T, s)
    if nargin < 2, s = 0; end
    R = T(1:3,1:3); pos = T(1:3,4);
    xf = 0.16; z_top = -0.05; z_wrist = -0.20; L = 0.25;
    th = s * (pi/2);
    tipR = [ xf + L*sin(th), 0, z_top + L*cos(th)];
    tipL = [-xf - L*sin(th), 0, z_top + L*cos(th)];
    P = [   0      0      z_wrist;     0      0      z_top;
          NaN    NaN    NaN;
           xf      0      z_top;     -xf      0      z_top;
          NaN    NaN    NaN;
           xf      0      z_top;   tipR(1) tipR(2) tipR(3);
          NaN    NaN    NaN;
          -xf      0      z_top;   tipL(1) tipL(2) tipL(3) ];
    W = (R*P' + pos)';
    h = plot3(W(:,1), W(:,2), W(:,3), 'k-', 'LineWidth', 4);
end

function s = gripper_open_frac(t_cycle)
%  - Προσέγγιση (t<2.3): ανοιχτή · Σύλληψη (2.3–2.8): κλείνει
%  - Μεταφορά/Συναρμολόγηση (2.8–8.75): κλειστή
%  - Απελευθέρωση (8.75–9.25): ανοίγει · Επιστροφή (t>9.25): ανοιχτή
    if     t_cycle < 2.3,  s = 1;
    elseif t_cycle < 2.8,  s = 1 - (t_cycle - 2.3)/0.5;
    elseif t_cycle < 8.75, s = 0;
    elseif t_cycle < 9.25, s = (t_cycle - 8.75)/0.5;
    else,                  s = 1;
    end
end