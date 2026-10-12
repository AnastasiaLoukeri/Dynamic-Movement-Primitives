% ΜΕΡΟΣ 2 — DMP με Χωρική Κλιμάκωση (Ερωτήματα 1, 2, 3)

clear; clc; close all;

%  Βασικές Παράμετροι & Σημεία (Μέρος 1) 
dt = 0.001;
P0     = [0.707, -0.293, 1.85, pi/4];
P1     = [0.19,  -1.00,  1.63,  0   ];   % σύλληψη @t=2.3s: x = -0.5 + 0.3*2.3
P2     = [0.34,  -1.00,  1.63,  0   ];   % τέλος συλληψης @t=2.8s: x = -0.5 + 0.3*2.8
P_VIA1 = [0.7,   -1.00,  1.85,  0   ];   
P3     = [0.08,   1.00,  1.85,  0   ];
P4     = [0.08,   1.00,  1.15,  pi/2];

% Παράμετροι DMP
alpha_z = 25; beta_z = alpha_z / 4; alpha_x = 1.0; n_bfs = 40;
[c, h_bfs] = init_bfs(n_bfs);


% ΕΡΩΤΗΜΑ 1: Εκπαίδευση DMP Θέσης (Τμήμα 4: P_VIA1→P3, διάρκεια 3.4s)
%            ΚΑΙ DMP Προσανατολισμού (Τμήμα 5β: P3(4)→P4(4), διάρκεια 0.95s)
% Επιθυμητές τροχιές: QUINTIC (ίδιες με Μέρος 1)

tau_A = 3.4; t_A = (0:dt:tau_A)';
w_dmp_A = zeros(n_bfs, 2);
Y_des_A = zeros(length(t_A), 2);

for j = 1:2
    % quintic
    Y_des_A(:,j) = quintic_full(P_VIA1(j), P3(j), 0, 0, 0, 0, 0, tau_A, t_A);
    dY_des  = [0; diff(Y_des_A(:,j))/dt];
    ddY_des = [0; diff(dY_des)/dt];
    w_dmp_A(:,j) = train_dmp(Y_des_A(:,j), dY_des, ddY_des, t_A, tau_A, ...
                              alpha_z, beta_z, alpha_x, c, h_bfs);
end

% DMP Προσανατολισμού (Τμήμα 5: γωνία φ)
tau_B = 0.95; t_B = (0:dt:tau_B)';
Phi_des_B = quintic_full(P3(4), P4(4), 0, 0, 0, 0, 0, tau_B, t_B);
dPhi_des  = [0; diff(Phi_des_B)/dt];
ddPhi_des = [0; diff(dPhi_des)/dt];
w_dmp_Phi = train_dmp(Phi_des_B, dPhi_des, ddPhi_des, t_B, tau_B, ...
                      alpha_z, beta_z, alpha_x, c, h_bfs);

% Σενάρια
deltas = [
     0.0,   0.0,  0.0;
    -0.2,   0.2,  0.0;
     0.4,  -0.15, 0.0
];
theta_deltas = [0, 20*(pi/180), -25*(pi/180)];


% ΕΡΩΤΗΜΑ 2: 3 Κύκλοι — DMP ΜΟΝΟ για Θέση, χωρίς DMP Προσανατολισμού

fprintf(' Ερώτημα 2 \n');
[cart_Q2, q_Q2, DMP_X_Q2, DMP_Y_Q2, ~] = simulate_3_cycles(...
    deltas, zeros(1,3), P0, P1, P2, P_VIA1, P3, P4, dt, ...
    tau_A, tau_B, w_dmp_A, [], alpha_z, beta_z, alpha_x, c, h_bfs);

t_Q2    = ((0:size(cart_Q2,1)-1)*dt)';
v_Q2    = diff(cart_Q2(:,1:3))/dt;
tv_Q2   = t_Q2(1:end-1);
qd_Q2   = diff(q_Q2)/dt;

% Σχήμα Q2-1: DMP Σύγκριση (Επιθυμητή - Παραγόμενη)
figure('Name','Q2: DMP Θέσης — Σύγκριση','Position',[50,50,900,500]);
labels = {'\delta_1=0', '\delta_2=[-0.2,0.2]', '\delta_3=[0.4,-0.15]'};
subplot(2,1,1); hold on; grid on;
plot(t_A, Y_des_A(:,1), '--k', 'LineWidth', 2, 'DisplayName', 'Επιθυμητή');
colors = {'r','b','g'};
for k=1:3, plot(t_A, DMP_X_Q2(:,k), colors{k}, 'LineWidth', 1.5, 'DisplayName', labels{k}); end
title('X: Επιθυμητή vs Παραγόμενη (DMP)'); xlabel('Χρόνος (s)'); ylabel('x (m)');
legend('Location','best');

subplot(2,1,2); hold on; grid on;
plot(t_A, Y_des_A(:,2), '--k', 'LineWidth', 2, 'DisplayName', 'Επιθυμητή');
for k=1:3, plot(t_A, DMP_Y_Q2(:,k), colors{k}, 'LineWidth', 1.5, 'DisplayName', labels{k}); end
title('Y: Επιθυμητή vs Παραγόμενη (DMP)'); xlabel('Χρόνος (s)'); ylabel('y (m)');
legend('Location','best');

%  Σχήμα Q2-2: Καρτεσιανή Θέση TCP 
figure('Name','Q2: Καρτεσιανή Θέση TCP','Position',[100,100,900,400]);
plot(t_Q2, cart_Q2(:,1),'r', t_Q2, cart_Q2(:,2),'g', t_Q2, cart_Q2(:,3),'b','LineWidth',1.5);
hold on; xline(12,'--k','LineWidth',1.2); xline(24,'--k','LineWidth',1.2);
title('Θέση TCP — 3 Κύκλοι (Ερώτημα 2)'); xlabel('Χρόνος (s)'); ylabel('Θέση (m)');
legend('x','y','z','Location','best'); grid on;

%  Σχήμα Q2-3: Καρτεσιανή Ταχύτητα TCP 
figure('Name','Q2: Καρτεσιανή Ταχύτητα TCP','Position',[150,150,900,400]);
plot(tv_Q2, v_Q2(:,1),'r', tv_Q2, v_Q2(:,2),'g', tv_Q2, v_Q2(:,3),'b','LineWidth',1.5);
hold on; xline(12,'--k','LineWidth',1.2); xline(24,'--k','LineWidth',1.2);
title('Ταχύτητα TCP — 3 Κύκλοι (Ερώτημα 2)'); xlabel('Χρόνος (s)'); ylabel('Ταχύτητα (m/s)');
legend('Vx','Vy','Vz','Location','best'); grid on; ylim([-2 2]);

%  Σχήμα Q2-4: Θέσεις Αρθρώσεων 
figure('Name','Q2: Θέσεις Αρθρώσεων','Position',[200,200,900,400]);
plot(t_Q2, q_Q2,'LineWidth',1.5);
hold on; xline(12,'--k','LineWidth',1.2); xline(24,'--k','LineWidth',1.2);
title('Θέσεις Αρθρώσεων q(t) — 3 Κύκλοι (Ερώτημα 2)'); xlabel('Χρόνος (s)'); ylabel('Θέση (rad ή m)');
legend('q_1','q_2','q_3','q_4','Location','best'); grid on;

%  Σχήμα Q2-5: Ταχύτητες Αρθρώσεων 
figure('Name','Q2: Ταχύτητες Αρθρώσεων','Position',[250,250,900,400]);
plot(tv_Q2, qd_Q2,'LineWidth',1.5);
hold on; xline(12,'--k','LineWidth',1.2); xline(24,'--k','LineWidth',1.2);
title('Ταχύτητες Αρθρώσεων dq/dt — 3 Κύκλοι (Ερώτημα 2)'); xlabel('Χρόνος (s)'); ylabel('Ταχύτητα (rad/s ή m/s)');
legend('dq_1/dt','dq_2/dt','dq_3/dt','dq_4/dt','Location','best'); grid on;


% ΕΡΩΤΗΜΑ 3: 3 Κύκλοι — DMP για Θέση ΚΑΙ Προσανατολισμό

fprintf(' Ερώτημα 3 \n');
[cart_Q3, q_Q3, DMP_X_Q3, DMP_Y_Q3, DMP_Phi_Q3] = simulate_3_cycles(...
    deltas, theta_deltas, P0, P1, P2, P_VIA1, P3, P4, dt, ...
    tau_A, tau_B, w_dmp_A, w_dmp_Phi, alpha_z, beta_z, alpha_x, c, h_bfs);

t_Q3    = ((0:size(cart_Q3,1)-1)*dt)';
v_Q3    = diff(cart_Q3(:,1:3))/dt;
tv_Q3   = t_Q3(1:end-1);
qd_Q3   = diff(q_Q3)/dt;

%  Σχήμα Q3-1: DMP Θέσης Σύγκριση 
figure('Name','Q3: DMP Θέσης — Σύγκριση','Position',[50,50,900,500]);
labels3 = {'\theta_{\delta1}=0°', '\theta_{\delta2}=20°', '\theta_{\delta3}=-25°'};
subplot(2,1,1); hold on; grid on;
plot(t_A, Y_des_A(:,1), '--k', 'LineWidth', 2, 'DisplayName', 'Επιθυμητή');
for k=1:3, plot(t_A, DMP_X_Q3(:,k), colors{k}, 'LineWidth', 1.5, 'DisplayName', labels3{k}); end
title('X: Επιθυμητή vs Παραγόμενη (DMP)'); xlabel('Χρόνος (s)'); ylabel('x (m)');
legend('Location','best');

subplot(2,1,2); hold on; grid on;
plot(t_A, Y_des_A(:,2), '--k', 'LineWidth', 2, 'DisplayName', 'Επιθυμητή');
for k=1:3, plot(t_A, DMP_Y_Q3(:,k), colors{k}, 'LineWidth', 1.5, 'DisplayName', labels3{k}); end
title('Y: Επιθυμητή vs Παραγόμενη (DMP)'); xlabel('Χρόνος (s)'); ylabel('y (m)');
legend('Location','best');

%  Σχήμα Q3-2: DMP Προσανατολισμού Σύγκριση 
figure('Name','Q3: DMP Προσανατολισμού — Σύγκριση','Position',[100,100,900,400]);
hold on; grid on;
plot(t_B, Phi_des_B, '--k', 'LineWidth', 2, 'DisplayName', 'Επιθυμητή');
for k=1:3, plot(t_B, DMP_Phi_Q3(:,k), colors{k}, 'LineWidth', 1.5, 'DisplayName', labels3{k}); end
title('φ: Επιθυμητή vs Παραγόμενη (DMP Προσανατολισμού)');
xlabel('Χρόνος (s)'); ylabel('φ (rad)'); legend('Location','best');

%  Σχήμα Q3-3: Καρτεσιανή Θέση TCP 
figure('Name','Q3: Καρτεσιανή Θέση TCP','Position',[150,150,900,400]);
plot(t_Q3, cart_Q3(:,1),'r', t_Q3, cart_Q3(:,2),'g', t_Q3, cart_Q3(:,3),'b','LineWidth',1.5);
hold on; xline(12,'--k','LineWidth',1.2); xline(24,'--k','LineWidth',1.2);
title('Θέση TCP — 3 Κύκλοι (Ερώτημα 3)'); xlabel('Χρόνος (s)'); ylabel('Θέση (m)');
legend('x','y','z','Location','best'); grid on;

%  Σχήμα Q3-4: Καρτεσιανή Ταχύτητα TCP 
figure('Name','Q3: Καρτεσιανή Ταχύτητα TCP','Position',[200,200,900,400]);
plot(tv_Q3, v_Q3(:,1),'r', tv_Q3, v_Q3(:,2),'g', tv_Q3, v_Q3(:,3),'b','LineWidth',1.5);
hold on; xline(12,'--k','LineWidth',1.2); xline(24,'--k','LineWidth',1.2);
title('Ταχύτητα TCP — 3 Κύκλοι (Ερώτημα 3)'); xlabel('Χρόνος (s)'); ylabel('Ταχύτητα (m/s)');
legend('Vx','Vy','Vz','Location','best'); grid on; ylim([-2 2]);

%  Σχήμα Q3-5: Θέσεις Αρθρώσεων 
figure('Name','Q3: Θέσεις Αρθρώσεων','Position',[250,250,900,400]);
plot(t_Q3, q_Q3,'LineWidth',1.5);
hold on; xline(12,'--k','LineWidth',1.2); xline(24,'--k','LineWidth',1.2);
title('Θέσεις Αρθρώσεων q(t) — 3 Κύκλοι (Ερώτημα 3)'); xlabel('Χρόνος (s)'); ylabel('Θέση (rad ή m)');
legend('q_1','q_2','q_3','q_4','Location','best'); grid on;

%  Σχήμα Q3-6: Ταχύτητες Αρθρώσεων
figure('Name','Q3: Ταχύτητες Αρθρώσεων','Position',[300,300,900,400]);
plot(tv_Q3, qd_Q3,'LineWidth',1.5);
hold on; xline(12,'--k','LineWidth',1.2); xline(24,'--k','LineWidth',1.2);
title('Ταχύτητες Αρθρώσεων dq/dt — 3 Κύκλοι (Ερώτημα 3)'); xlabel('Χρόνος (s)'); ylabel('Ταχύτητα (rad/s ή m/s)');
legend('dq_1/dt','dq_2/dt','dq_3/dt','dq_4/dt','Location','best'); grid on;

% =========================================================================
%   Προσανατολισμός άκρου φ (γωνία TCP) — 3 Κύκλοι
%  φ = cart(:,4). Δείχνει το βίδωμα 0 -> (π/2 + θδ) σε κάθε κύκλο.
% =========================================================================

% --- Ερώτημα 2 (θδ=0): ίδιο βίδωμα 0->π/2 κάθε κύκλο ---
figure('Name','Q2: Προσανατολισμός TCP \phi','Position',[300,300,900,350]);
plot(t_Q2, cart_Q2(:,4),'m','LineWidth',1.6);
hold on; xline(12,'--k','LineWidth',1.2); xline(24,'--k','LineWidth',1.2);
yline(pi/2,':k','\phi = 90°','LabelHorizontalAlignment','left');
title('Προσανατολισμός TCP \phi — 3 Κύκλοι (Ερώτημα 2)');
xlabel('Χρόνος (s)'); ylabel('\phi (rad)'); grid on; ylim([-0.2 2.2]);

% --- Ερώτημα 3 (θδ = 0, 20°, -25°): βίδωμα 0 -> 90°/110°/65° ---
phi_goals = pi/2 + theta_deltas;            % στόχοι φ ανά κύκλο (rad)
figure('Name','Q3: Προσανατολισμός TCP \phi','Position',[350,350,900,350]);
plot(t_Q3, cart_Q3(:,4),'m','LineWidth',1.6);
hold on; xline(12,'--k','LineWidth',1.2); xline(24,'--k','LineWidth',1.2);
for k = 1:3
    yline(phi_goals(k),':', ...
        sprintf('\\phi_{goal,%d} = %.0f°', k, rad2deg(phi_goals(k))), ...
        'LabelHorizontalAlignment','left');
end
title('Προσανατολισμός TCP \phi — 3 Κύκλοι (Ερώτημα 3)');
xlabel('Χρόνος (s)'); ylabel('\phi (rad)'); grid on; ylim([-0.2 2.2]);
% ΚΥΡΙΑ ΣΥΝΑΡΤΗΣΗ ΠΡΟΣΟΜΟΙΩΣΗΣ 3 ΚΥΚΛΩΝ

function [cart_traj, q_traj, sim_X_all, sim_Y_all, sim_Phi_all] = simulate_3_cycles(...
    deltas, theta_deltas, P0, P1, P2, P_VIA1, P3, P4, dt, ...
    tau_A, tau_B, w_dmp_A, w_dmp_Phi, alpha_z, beta_z, alpha_x, c, h_bfs)

    cart_traj = [];
    sim_X_all   = zeros(round(tau_A/dt)+1, 3);
    sim_Y_all   = zeros(round(tau_A/dt)+1, 3);
    sim_Phi_all = zeros(round(tau_B/dt)+1, 3);

    for cycle = 1:3
        % Στόχοι με offset
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
        z_rim = 1.21;     % TCP με μύτη πείρου στο χείλος υποδοχής (1.08 + 0.13)

        %  Τμήμα 5α: ευθεία κάθοδος ως το χείλος (7.0→7.8s), φ = g_P3(4)
        t5a = (7.0:dt:7.8)'; n5a = length(t5a);
        tr5a_Z = quintic_full(g_P3(3), z_rim, 0, 0, 0, 0, 7.0, 7.8, t5a);
        tr5a = [g_P3(1)*ones(n5a,1), g_P3(2)*ones(n5a,1), tr5a_Z, g_P3(4)*ones(n5a,1)];

        %  Τμήμα 5β: βίδωμα (7.8→8.75s) — τελική κάθοδος 0.06m + στροφή 90° 
        t5b = (7.8:dt:8.75)'; n5b = length(t5b);
        tr5b_Z = quintic_full(z_rim, g_P4(3), 0, 0, 0, 0, 7.8, 8.75, t5b);
        if isempty(w_dmp_Phi)
            % Ερώτημα 2: κλασικό Quintic για φ
            sim_Phi = quintic_full(P3(4), g_P4(4), 0, 0, 0, 0, 7.8, 8.75, t5b);
        else
            % Ερώτημα 3: DMP για φ
            [sim_Phi, ~, ~] = rollout_dmp(g_P3(4), g_P4(4), w_dmp_Phi, tau_B, dt, alpha_z, beta_z, alpha_x, c, h_bfs);
            sim_Phi = sim_Phi + (g_P4(4)-sim_Phi(end))*linspace(0,1,numel(sim_Phi))';
        end
        sim_Phi_all(:,cycle) = sim_Phi;
        tr5b = [g_P4(1)*ones(n5b,1), g_P4(2)*ones(n5b,1), tr5b_Z, sim_Phi];

        %  Τμήμα 6: Pause (8.75-9.25s) 
        t6 = (8.75:dt:9.25)'; n6 = length(t6);
        tr6 = repmat(g_P4, n6, 1);

        %  Τμήμα 7: Επιστροφή με fly-by (9.25→12s) 
        P_RETRACT = [g_P4(1), g_P4(2), 1.85, g_P4(4)];
        P_VIA_RET = [1.00, 0.50, 1.85, 3*pi/8];

        %  Υπολογισμός v_via για σωστό fly-by
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
            tr7a(:,j) = quintic_full(g_P4(j),     P_RETRACT(j),  0,        0,        0, 0,  9.25,  10.05, t7a);
            tr7b(:,j) = quintic_full(P_RETRACT(j), P_VIA_RET(j), 0,        v_via(j), 0, 0, 10.05,  11.15, t7b);
            tr7c(:,j) = quintic_full(P_VIA_RET(j), P0(j),        v_via(j), 0,        0, 0, 11.15,  12.0, t7c);
        end

        cycle_traj = [tr1(1:end-1,:); tr2(1:end-1,:); tr3(1:end-1,:);
                      tr4(1:end-1,:); tr5a(1:end-1,:); tr5b(1:end-1,:); tr6(1:end-1,:);
                      tr7a(1:end-1,:); tr7b(1:end-1,:); tr7c];
        cart_traj = [cart_traj; cycle_traj(1:end-1,:)];
    end
    cart_traj = [cart_traj; P0];

    %  Αντίστροφη Κινηματική 
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
        y = y + z*(dt/tau);        % Euler: πρώτα y με παλιό z
        z = z + dz*(dt/tau);       % μετά update z
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