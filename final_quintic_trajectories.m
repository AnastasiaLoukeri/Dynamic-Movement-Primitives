% =========================================================================
% Πλήρως Καρτεσιανός Σχεδιασμός με C2 Συνέχεια — Quintic Παντού
% =========================================================================
clear; clc; close all;
set(0, 'defaultfigurecolor', [1 1 1])

%  Αρχικοποίηση Ρομπότ 
[scara, q0] = create_robot();
q0 = q0(:)';
dt = 0.001;

% Σημεία
P0     = [0.707, -0.293, 1.85, pi/4];
P1     = [0.19,  -1.00,  1.63,  0   ];   % σύλληψη @t=2.3s: x = -0.5 + 0.3*2.3
P2     = [0.34,  -1.00,  1.63,  0   ];   % τέλος συλληψης @t=2.8s: x = -0.5 + 0.3*2.8
P_VIA1 = [0.7, -1.00, 1.85, 0];   
P3     = [0.08,   1.00,  1.85,  0   ];
P4     = [0.08,   1.00,  1.15,  pi/2];
P5     = [0.08,   1.00,  1.15,  pi/2];

% =========================================================================
% ΖΗΤΟΥΜΕΝΟ 1: Ομογενείς Μετασχηματισμοί των Ενδιάμεσων Πόζων (πλαίσιο W)
% Για κάθε πόζα [x y z phi] υπολογίζεται q μέσω αντίστροφης κινηματικής
% και T = scara.fkine(q), ώστε ο μετασχηματισμός να συμπίπτει ακριβώς με
% την πόζα που υλοποιεί το προσομοιωμένο ρομπότ.
% =========================================================================
P_RETRACT = [P5(1), P5(2), 1.85, P5(4)];
P_VIA_RET = [1.00,  0.50,  1.85, 3*pi/8];

pose_names = {'P0','P1','P2','P_VIA1','P3','P4','P5','P_RETRACT','P_VIA_RET'};
pose_list  = [P0; P1; P2; P_VIA1; P3; P4; P5; P_RETRACT; P_VIA_RET];

fprintf('\n ΟΜΟΓΕΝΕΙΣ ΜΕΤΑΣΧΗΜΑΤΙΣΜΟΙ ΕΝΔΙΑΜΕΣΩΝ ΠΟΖΩΝ \n');
for k = 1:size(pose_list,1)
    p  = pose_list(k,:);
    qk = ik_scara(p);
    Tk = scara.fkine(qk);
    Tm = Tk.T;
    fprintf('\n%-9s  (x=%.3f, y=%.3f, z=%.3f, phi=%.4f rad)\n', ...
            pose_names{k}, p(1), p(2), p(3), p(4));
    fprintf('  q = [%.4f  %.4f  %.4f  %.4f]\n', qk(1), qk(2), qk(3), qk(4));
    for r = 1:4
        fprintf('   % .4f  % .4f  % .4f  % .4f\n', Tm(r,1), Tm(r,2), Tm(r,3), Tm(r,4));
    end
end


cart_traj_1cycle = [];

% =========================================================================
% ΤΜΗΜΑ 1: P0 → P1 (0 → 2.3s)
% Προσέγγιση Part A. Vx(tf)=0.3 m/s για συνέχεια στη σύλληψη.
% Z: παραμένει στο 1.85 έως t=1.3s, κατεβαίνει μόνο στο τελευταίο 1s.
% =========================================================================
t_sec = (0:dt:2.3)';
traj  = zeros(length(t_sec), 4);
traj(:,1) = quintic_full(P0(1), P1(1), 0, 0.3, 0, 0, 0,   2.3, t_sec);
traj(:,2) = quintic_full(P0(2), P1(2), 0, 0.0, 0, 0, 0,   2.3, t_sec);
traj(:,4) = quintic_full(P0(4), P1(4), 0, 0.0, 0, 0, 0,   2.3, t_sec);

idx_z = t_sec >= 1.3;
traj(:,3) = P0(3) * ones(length(t_sec), 1);
traj(idx_z, 3) = quintic_full(P0(3), P1(3), 0, 0, 0, 0, 1.3, 2.3, t_sec(idx_z));
cart_traj_1cycle = [cart_traj_1cycle; traj(1:end-1,:)];

% =========================================================================
% ΤΜΗΜΑ 2: P1 → P2 (2.3 → 2.8s)  —  Σταθερή ταχύτητα (Σύλληψη)
% Ο ρομπότ ακολουθεί τον ταινιόδρομο Α με v=0.3 m/s
% =========================================================================
t_sec = (2.3:dt:2.8)';
traj  = zeros(length(t_sec), 4);
traj(:,1) = P1(1) + 0.3*(t_sec - 2.3);
for j = 2:4, traj(:,j) = P1(j) * ones(length(t_sec), 1); end
cart_traj_1cycle = [cart_traj_1cycle; traj(1:end-1,:)];

% =========================================================================
% ΤΜΗΜΑ 3: P2 → P_VIA1 (2.8 → 3.6s)  —  Ανύψωση Z, επιβράδυνση X [QUINTIC]
% Vx(t0)=0.3 m/s (συνέχεια από Τμήμα 2), Vx(tf)=0
% =========================================================================


t_sec = (2.8:dt:3.6)';
traj  = zeros(length(t_sec), 4);
traj(:,1) = quintic_full(P2(1), P_VIA1(1), 0.3, 0, 0, 0, 2.8, 3.6, t_sec);
traj(:,2) = P2(2) * ones(length(t_sec), 1);
traj(:,3) = quintic_full(P2(3), P_VIA1(3), 0,   0, 0, 0, 2.8, 3.6, t_sec);
traj(:,4) = P2(4) * ones(length(t_sec), 1);
cart_traj_1cycle = [cart_traj_1cycle; traj(1:end-1,:)];

% =========================================================================
% ΤΜΗΜΑ 4: P_VIA1 → P3 (3.6 → 7.0s)  —  Οριζόντια Μεταφορά [QUINTIC]
% C2 συνέχεια, μηδενικές ταχύτητες/επιταχύνσεις στα άκρα
% =========================================================================



t_sec = (3.6:dt:7.0)';
traj  = zeros(length(t_sec), 4);
for j = 1:4
    traj(:,j) = quintic_full(P_VIA1(j), P3(j), 0, 0, 0, 0, 3.6, 7.0, t_sec);
end
cart_traj_1cycle = [cart_traj_1cycle; traj(1:end-1,:)];

% =========================================================================
% ΤΜΗΜΑ 5: P3 → P4 (7.0 → 8.75s)  —  Συναρμολόγηση (5α κάθοδος 7.0-7.8 + 5β βίδωμα 7.8-8.75)
% =========================================================================


z_rim = 1.21;     % TCP με μύτη πείρου στο χείλος υποδοχής (1.08 + 0.13)

% 5α: ευθεία κάθοδος ως το χείλος [7.4, 8.0], φ=0
t_sec = (7.0:dt:7.8)';
traj  = zeros(length(t_sec), 4);
traj(:,1) = P3(1);  traj(:,2) = P3(2);
traj(:,3) = quintic_full(P3(3), z_rim, 0,0,0,0, 7.0, 7.8, t_sec);
traj(:,4) = P3(4);                 % φ = 0
cart_traj_1cycle = [cart_traj_1cycle; traj(1:end-1,:)];

% 5β: βίδωμα [8.0, 9.0] — τελική κάθοδος 0.05m + στροφή 90°
t_sec = (7.8:dt:8.75)';
traj  = zeros(length(t_sec), 4);
traj(:,1) = P4(1);  traj(:,2) = P4(2);
traj(:,3) = quintic_full(z_rim, P4(3), 0,0,0,0, 7.8, 8.75, t_sec);
traj(:,4) = quintic_full(P3(4), P4(4), 0,0,0,0, 7.8, 8.75, t_sec);
cart_traj_1cycle = [cart_traj_1cycle; traj(1:end-1,:)];

% =========================================================================
% ΤΜΗΜΑ 6: P4 → P5 (8.75 → 9.25s)  —  Pause (Άνοιγμα αρπάγης)
% =========================================================================
t_sec = (8.75:dt:9.25)';
traj  = repmat(P4, length(t_sec), 1);
cart_traj_1cycle = [cart_traj_1cycle; traj(1:end-1,:)];

% =========================================================================
% ΤΜΗΜΑ 7: P5 → P0 (9.25 → 12.0s)  —  Επιστροφή με "Fly-by" Via-Point
% =========================================================================
P_RETRACT = [P5(1), P5(2), 1.85, P5(4)];
P_VIA_RET = [1.00,  0.50,  1.85, 3*pi/8];

t_7b = 11.15 - 10.05;
t_7c = 12.0 - 11.15;
v_via = zeros(1, 4);
for j = 1:4
    V_prev = (P_VIA_RET(j) - P_RETRACT(j)) / t_7b;
    V_next = (P0(j)        - P_VIA_RET(j)) / t_7c;
    if sign(V_prev) == sign(V_next)
        v_via(j) = sign(V_prev) * sqrt(abs(V_prev * V_next));
    else
        v_via(j) = 0;
    end
end

% Φάση 7α (9.25 - 10.05s): Ανύψωση Z [QUINTIC]
t_sec = (9.25:dt:10.05)';
traj  = zeros(length(t_sec), 4);
for j = 1:4
    traj(:,j) = quintic_full(P5(j), P_RETRACT(j), 0, 0, 0, 0, 9.25, 10.05, t_sec);
end
cart_traj_1cycle = [cart_traj_1cycle; traj(1:end-1,:)];

% Φάση 7β (10.05 - 11.15s): P_RETRACT → P_VIA_RET [QUINTIC]
t_sec = (10.05:dt:11.15)';
traj  = zeros(length(t_sec), 4);
for j = 1:4
    traj(:,j) = quintic_full(P_RETRACT(j), P_VIA_RET(j), 0, v_via(j), 0, 0, 10.05, 11.15, t_sec);
end
cart_traj_1cycle = [cart_traj_1cycle; traj(1:end-1,:)];

% Φάση 7γ (11.15 - 12.0s): P_VIA_RET → P0 [QUINTIC]
t_sec = (11.15:dt:12.0)';
traj  = zeros(length(t_sec), 4);
for j = 1:4
    traj(:,j) = quintic_full(P_VIA_RET(j), P0(j), v_via(j), 0, 0, 0, 11.15, 12.0, t_sec);
end
cart_traj_1cycle = [cart_traj_1cycle; traj];

% =========================================================================
% ΛΟΥΠΑ 2 ΚΥΚΛΩΝ (24s)
% =========================================================================
cart_traj = [cart_traj_1cycle(1:end-1, :); cart_traj_1cycle];
N      = size(cart_traj, 1);
t_traj = ((0:N-1) * dt)';

% =========================================================================
% ΑΝΤΙΣΤΡΟΦΗ ΚΙΝΗΜΑΤΙΚΗ
% =========================================================================
q_traj = zeros(N, scara.n);
L1 = 1; L2 = 1;
for i = 1:N
    x = cart_traj(i,1); y = cart_traj(i,2);
    z = cart_traj(i,3); phi = cart_traj(i,4);

    D  = (x^2 + y^2 - L1^2 - L2^2) / (2 * L1 * L2);
    q2 = atan2(sqrt(1-D^2), D);   % Elbow-up
    q1 = atan2(y, x) - atan2(L2*sin(q2), L1 + L2*cos(q2));
    q3 = 1.95 - z;
    q4 = q1 + q2 - phi;

    q_traj(i,:) = [q1, q2, q3, q4];
end

% =========================================================================
% ANIMATION (2 Κύκλοι)
% =========================================================================
figure('Name','SCARA Assembly (2 Cycles)','Position',[100,100,900,700]);
view(3); grid on; hold on;
axis([-1.5 3.5 -2.5 2.5 0 3]);
xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)');

[xa,ya] = meshgrid(-1.5:0.5:3.0, -1.5:0.1:-0.5);
surf(xa,ya,1.5*ones(size(xa)),'FaceColor',[0.7 0.8 0.9],'EdgeAlpha',0.1);
[xb,yb] = meshgrid(-1.5:0.5:3.0, 0.5:0.1:1.5);
surf(xb,yb,1.0*ones(size(xb)),'FaceColor',[0.9 0.8 0.7],'EdgeAlpha',0.1);

h_partA=[]; h_partB=[]; h_gripper=[];

for i = 1:15:N
    t = t_traj(i);
    T_tcp   = scara.fkine(q_traj(i,:));
    t_cycle = mod(t, 12.0);

    if t_cycle <= 2.8
        pos_A = [-0.5+0.3*t_cycle, -1, 1.5]; rot_A = eye(3);
    elseif t_cycle <= 8.75
        pos_A = T_tcp.t' - [0,0,0.13]; rot_A = t2r(trotz(cart_traj(i,4)));
    else
        pos_A = [0.08, 1, 1.02]; rot_A = t2r(trotz(pi/2));
    end

    if ~isempty(h_partA), delete(h_partA); end
    h_partA = draw_rotated_box(pos_A,[0.25,0.25,0.18],[0.8 0.2 0.2],rot_A);

    x_B = (t_cycle<=2.4)*(2.0-0.8*t_cycle) + (t_cycle>2.4)*0.08;
    if ~isempty(h_partB), delete(h_partB); end
    h_partB = draw_rotated_box([x_B,1,1.0],[0.3,0.3,0.1],[0.2 0.8 0.2],eye(3));

    if ~isempty(h_gripper), delete(h_gripper); end
    s_grip = gripper_open_frac(t_cycle);
    h_gripper = draw_gripper(T_tcp.T, s_grip);

    scara.plot(q_traj(i,:),'workspace',[-1.5 3.5 -2.5 2.5 0 3],'delay',0);
    drawnow;
end

% =========================================================================
% ΔΙΑΓΡΑΜΜΑΤΑ (24s)
% =========================================================================
t_cv    = t_traj(1:end-1);
v_cart  = diff(cart_traj(:,1:3)) / dt;
qd_traj = diff(q_traj) / dt;

t_ca   = t_cv(1:end-1);
a_cart = diff(v_cart) / dt;
qd_acc = diff(qd_traj) / dt;

% seg_times = [0, 3.0, 3.5, 5.0, 8.0, 9.0, 9.5, 12.0];
seg_times = [0, 2.3, 2.8, 3.6, 7.0, 7.8, 8.75, 9.25, 10.05, 11.15, 12.0];

% Σχήμα 1: Καρτεσιανή Θέση TCP
figure('Name','Καρτεσιανή Θέση TCP','Position',[50,50,900,500]);
subplot(2,1,1);
plot(t_traj, cart_traj(:,1),'r', t_traj, cart_traj(:,2),'g', ...
     t_traj, cart_traj(:,3),'b','LineWidth',1.5);
for ts = seg_times, xline(ts,'--k','LineWidth',0.8); end
xline(12.0,'-k','LineWidth',1.5,'Label','Τέλος 1ου Κύκλου');
title('Θέση TCP (2 Κύκλοι)','FontSize',13);
xlabel('Χρόνος (s)'); ylabel('Θέση (m)');
legend('x','y','z','Location','best'); grid on;

subplot(2,1,2);
plot(t_traj, cart_traj(:,4),'m','LineWidth',1.5);
for ts = seg_times, xline(ts,'--k','LineWidth',0.8); end
xline(12.0,'-k','LineWidth',1.5);
title('Γωνία Εργαλείου φ (2 Κύκλοι)','FontSize',13);
xlabel('Χρόνος (s)'); ylabel('φ (rad)'); grid on;

% Σχήμα 2: Γραμμική Ταχύτητα TCP
figure('Name','Γραμμική Ταχύτητα TCP','Position',[100,100,900,400]);
plot(t_cv, v_cart(:,1),'r', t_cv, v_cart(:,2),'g', ...
     t_cv, v_cart(:,3),'b','LineWidth',1.8);
hold on;
for ts = seg_times, xline(ts,'--k','LineWidth',0.8); end
xline(12.0,'-k','LineWidth',1.5);
title('Γραμμική Ταχύτητα TCP (2 Κύκλοι)','FontSize',13);
xlabel('Χρόνος (s)'); ylabel('Ταχύτητα (m/s)');
legend('Vx','Vy','Vz','Location','best'); grid on; ylim([-2 2]);

% Σχήμα 3: Ταχύτητες Αρθρώσεων
figure('Name','Ταχύτητες Αρθρώσεων','Position',[150,150,900,400]);
plot(t_cv, qd_traj,'LineWidth',1.5);
for ts = seg_times, xline(ts,'--k','LineWidth',0.8); end
xline(12.0,'-k','LineWidth',1.5,'Label','Τέλος 1ου Κύκλου');
title('Ταχύτητες Αρθρώσεων dq/dt (2 Κύκλοι)','FontSize',13);
xlabel('Χρόνος (s)'); ylabel('Ταχύτητα (rad/s ή m/s)');
legend('dq_1/dt','dq_2/dt','dq_3/dt','dq_4/dt','Location','best'); grid on;

% Σχήμα 4: Γραμμική Επιτάχυνση TCP
figure('Name','Γραμμική Επιτάχυνση TCP','Position',[200,200,900,400]);
plot(t_ca, a_cart(:,1),'r', t_ca, a_cart(:,2),'g', ...
     t_ca, a_cart(:,3),'b','LineWidth',1.5);
hold on;
for ts = seg_times, xline(ts,'--k','LineWidth',0.8); end
xline(12.0,'-k','LineWidth',1.5);
title('Γραμμική Επιτάχυνση TCP (2 Κύκλοι)','FontSize',13);
xlabel('Χρόνος (s)'); ylabel('Επιτάχυνση (m/s^2)');
legend('Ax','Ay','Az','Location','best'); grid on; ylim([-15 15]);

% Σχήμα 5: Επιτάχυνση Αρθρώσεων
figure('Name','Επιτάχυνση Αρθρώσεων','Position',[250,250,900,400]);
plot(t_ca, qd_acc,'LineWidth',1.5);
for ts = seg_times, xline(ts,'--k','LineWidth',0.8); end
xline(12.0,'-k','LineWidth',1.5);
title('Επιταχύνσεις Αρθρώσεων d^2q/dt^2 (2 Κύκλοι)','FontSize',13);
xlabel('Χρόνος (s)'); ylabel('Επιτάχυνση');
legend('ddq_1','ddq_2','ddq_3','ddq_4','Location','best'); grid on; ylim([-20 20]);

% Σχήμα 6: Θέσεις Αρθρώσεων
figure('Name','Θέσεις Αρθρώσεων','Position',[300,300,900,400]);
plot(t_traj, q_traj,'LineWidth',1.5);
hold on;
for ts = seg_times, xline(ts,'--k','LineWidth',0.8); end
xline(12.0,'-k','LineWidth',1.5,'Label','Τέλος 1ου Κύκλου');
title('Θέσεις Αρθρώσεων q(t) (2 Κύκλοι)','FontSize',13);
xlabel('Χρόνος (s)'); ylabel('Θέση (rad ή m)');
legend('q_1','q_2','q_3','q_4','Location','best'); grid on;

% =========================================================================
%  ΣΥΝΑΡΤΗΣΕΙΣ
% =========================================================================

function q_out = quintic_full(q0, qf, v0, vf, a0, af, t0, tf, t_arr)
% Πολυώνυμο 5ου βαθμού — C2 συνέχεια (θέση + ταχύτητα + επιτάχυνση).
    h = tf - t0;
    if h <= 0; q_out = q0 * ones(length(t_arr),1); return; end

    k0 = q0;
    k1 = v0;
    k2 = a0 / 2;

    A = qf - q0 - v0*h - (a0/2)*h^2;
    B = vf - v0 - a0*h;
    C = af - a0;

    k3 = (10*A - 4*B*h + 0.5*C*h^2) / h^3;
    k4 = (-15*A + 7*B*h - C*h^2)    / h^4;
    k5 = (6*A  - 3*B*h + 0.5*C*h^2) / h^5;

    q_out = zeros(length(t_arr), 1);
    for k = 1:length(t_arr)
        tau = max(0, min(t_arr(k) - t0, h));
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

function q = ik_scara(p)
% Αντίστροφη κινηματική SCARA (elbow-up) — ίδια με τον βρόχο IK του σεναρίου.
    x = p(1); y = p(2); z = p(3); phi = p(4);
    D  = (x^2 + y^2 - 2) / 2;
    q2 = atan2(sqrt(max(0,1-D^2)), D);
    q1 = atan2(y,x) - atan2(sin(q2), 1+cos(q2));
    q  = [q1, q2, 1.95 - z, q1 + q2 - phi];
end

function h = draw_gripper(T, s)
% Γκρίπερ με ανοιγόκλεισμα. s: 1 = πλήρως ανοιχτή (δάχτυλα παράλληλα στην
% παλάμη, Σχήμα 3), 0 = κλειστή (δάχτυλα κάτω, πιάνουν το αντικείμενο).
% Συντεταγμένες στο πλαίσιο εργαλείου: +z προς το αντικείμενο, -z προς τον καρπό.
    if nargin < 2, s = 0; end
    R = T(1:3,1:3); pos = T(1:3,4);
    xf = 0.16; z_top = -0.05; z_wrist = -0.20; L = 0.25;
    th = s * (pi/2);                 % 0 = κλειστή (κάτω), pi/2 = ανοιχτή (οριζόντια)
    tipR = [ xf + L*sin(th), 0, z_top + L*cos(th)];   % άκρο δεξιού δακτύλου
    tipL = [-xf - L*sin(th), 0, z_top + L*cos(th)];   % άκρο αριστερού δακτύλου
    P = [   0      0      z_wrist;     0      0      z_top;     % στέλεχος → παλάμη
          NaN    NaN    NaN;
           xf      0      z_top;     -xf      0      z_top;     % παλάμη
          NaN    NaN    NaN;
           xf      0      z_top;   tipR(1) tipR(2) tipR(3);     % δεξί δάχτυλο
          NaN    NaN    NaN;
          -xf      0      z_top;   tipL(1) tipL(2) tipL(3) ];   % αριστερό δάχτυλο
    W = (R*P' + pos)';
    h = plot3(W(:,1), W(:,2), W(:,3), 'k-', 'LineWidth', 4);
end

function s = gripper_open_frac(t_cycle)
% Κατάσταση αρπάγης μέσα στον κύκλο 12s (1 = ανοιχτή, 0 = κλειστή):
%  - Προσέγγιση (t<2.3): ανοιχτή
%  - Σύλληψη   (2.3–2.8): κλείνει σε 0.5s
%  - Μεταφορά/Συναρμολόγηση (2.8–8.75): κλειστή
%  - Απελευθέρωση (8.75–9.25): ανοίγει σε 0.5s
%  - Επιστροφή (t>9.25): ανοιχτή
    if     t_cycle < 2.3,  s = 1;
    elseif t_cycle < 2.8,  s = 1 - (t_cycle - 2.3)/0.5;
    elseif t_cycle < 8.75, s = 0;
    elseif t_cycle < 9.25, s = (t_cycle - 8.75)/0.5;
    else,                  s = 1;
    end
end