%% ���C���t�@�C��:
% ��ԕϐ��Ƒ���ϐ��𗣎U�����ē��I�v��@(discretized DP)������.

clear;
close all;
format short;

%% *** �J���u���[�V���� ***
beta  = 0.96; % �������q
gamma = 1.0;  % ���ΓI�댯���x(�َ��_�Ԃ̑�ւ̒e�͐��̋t��)
alpha = 0.40; % ���{���z��
delta = 1.00; % �Œ莑�{����(0.08)

% *** ���U���p�̃p�����[�^ ***
nk   = 10001;  % �O���b�h�̐�
kmax = 0.5;   % ���{�O���b�h�̍ő�l
%kmax = 10.0; % ���{�O���b�h�̍ő�l(�Œ莑�{����=0.08�̏ꍇ�Ɏg�p)
kmin = 0.05;  % ���{�O���b�h�̍ŏ��l (0�ɂ���Ɛ��Y���o���Ȃ��Ȃ�)
%========================

% *** �����̊ ***
it = 1;          % ���[�v�E�J�E���^�[
maxit = 1000;    % �J��Ԃ��v�Z�̍ő�l
tol  = 1.0e-005; % ���e�덷(STEP 2)
dif1 = 1;        % ���l�֐��̌J��Ԃ��덷
dif2 = 1.0;      % ����֐��̌J��Ԃ��덷
count = 1;
%==================

%% �v�Z�J�n

tic

disp('')
disp('-+- Solve a neoclassical growth model -+-');

%% STEP 1(a): �O���b�h����

kgrid = linspace(kmin, kmax, nk)';

%% STEP 1(b): ���l�֐��E����֐��̏����l��ݒ�

vfcn  = zeros(nk, 1);
pfcn  = zeros(nk, 1);
Tvfcn = zeros(nk, 1);
Tpfcn = zeros(nk, 1);
vkp   = zeros(nk, nk);
val_tmp = zeros(nk, 4);

%% STEP 3: ���p�֐��̑g�ݍ��킹

% ���p�֐��̏����l (���0�ȉ��ɂȂ�g�ݍ��킹�ɂ̓y�i���e�B)
util = -10000.0*ones(nk, nk);

% ������l�ɂȂ�(k,k')�̑g�ݍ��킹�ɂ��Č��p���v�Z
for i = 1:nk
    %  �����鑀��ϐ�k'�ɂ���:
    for j = 1:nk
        wealth = kgrid(i).^alpha + (1.0-delta).*kgrid(i);
        cons = wealth - kgrid(j);
        if cons > 0
           util(j,i) = CRRA(cons, gamma);
        end
    end
end

%% STEP 4: ���l�֐����J��Ԃ��v�Z

while it < maxit && dif1 > tol

    % �x���}��������: V(k;k')
    for i = 1:nk
        vkp(:,i) = util(:,i) + beta.*vfcn;
    end
    
    % �œK��: �ek�ɂ���V(k;k')���ő�ɂ���k'��T��
    [Tvfcn, ploc] = max(vkp);
    Tvfcn = Tvfcn';
    Tpfcn = kgrid(ploc);
    
    % �J��Ԃ��v�Z�덷���m�F
    dif1 = max(abs((Tvfcn-vfcn)./vfcn));
    dif2 = max(abs((Tpfcn-pfcn)./pfcn));

    % ���l�֐��E����֐����A�b�v�f�[�g
    vfcn = Tvfcn;
    pfcn = Tpfcn;
    fprintf('iteration index: %i, iteration diff of value: %d, iteration diff of policy: %d \n', it, dif1, dif2);

    % �����r���̌J��Ԃ��v�Z�덷��ۑ�
    % �r���o�߂�}������ړI�Ȃ̂ŁA�ʏ�͕s�v(�ނ���x���Ȃ�̂ŏ����ׂ�)
    v_conv(it) = dif1;
    p_conv(it) = dif2;

    % ���������l�֐��̎�����}������ړI�ŕۑ�(�{���͕s�v)
    if it==1 || it==3 || it==5
       val_tmp(:, count) = vfcn;
       count = count + 1;
    end
    
    it = it + 1;

end

val_tmp(:, 4) = vfcn;

toc

%% �v�Z���ʂ��R�}���h�E�B���h�E�ɕ\��

disp('-+- Parameter values -+-');
disp('');
fprintf('beta=%5.2f, gamma=%5.2f, alpha=%5.2f, delta=%5.2f \n', beta, gamma, alpha, delta);
disp(''); 
fprintf('kmin=%5.2f, kmax=%5.2f, #grid=%i \n', kmin, kmax, nk);
disp('');

%% ��͓I��

AA = (1.0-beta).^(-1) * (log(1.0-alpha*beta) + ((alpha*beta)/(1.0-alpha*beta))*log(alpha*beta));
BB = alpha/(1.0-alpha*beta);
v_true = AA + BB*log(kgrid);
p_true = beta*alpha*(kgrid.^alpha);

%% �I�C���[����������덷�𑪒�(�񎩓���)

nkk = 21;
kgrid2 = linspace(kmin, kmax, nkk)';
pfcn0 = zeros(nkk, 1);
j = 1;

pfcn0(1, 1) = pfcn(1, 1);
pfcn0(2, 1) = pfcn(501, 1);
pfcn0(3, 1) = pfcn(1001, 1);
pfcn0(4, 1) = pfcn(1501, 1);
pfcn0(5, 1) = pfcn(2001, 1);
pfcn0(6, 1) = pfcn(2501, 1);
pfcn0(7, 1) = pfcn(3001, 1);
pfcn0(8, 1) = pfcn(3501, 1);
pfcn0(9, 1) = pfcn(4001, 1);
pfcn0(10, 1) = pfcn(4501, 1);
pfcn0(11, 1) = pfcn(5001, 1);
pfcn0(12, 1) = pfcn(5501, 1);
pfcn0(13, 1) = pfcn(6001, 1);
pfcn0(14, 1) = pfcn(6501, 1);
pfcn0(15, 1) = pfcn(7001, 1);
pfcn0(16, 1) = pfcn(7501, 1);
pfcn0(17, 1) = pfcn(8001, 1);
pfcn0(18, 1) = pfcn(8501, 1);
pfcn0(19, 1) = pfcn(9001, 1);
pfcn0(20, 1) = pfcn(9501, 1);
pfcn0(21, 1) = pfcn(10001, 1);

cons = kgrid2.^alpha + (1.-delta).*kgrid2 - pfcn0(:,1);
LHS  = mu_CRRA(cons, gamma);
kp   = pfcn0(:,1);
kpp  = interp1(kgrid2, pfcn0(:,1), kp);
cons = kp.^alpha + (1.-delta).*kp - kpp;
rent = alpha.*kp.^(alpha-1.0) - delta;
RHS  = beta.*(1.+rent).*mu_CRRA(cons, gamma);
err  = RHS./LHS-1.0;

csvwrite("err_ddp.csv", err);

%% �}��`��

figure;
plot(kgrid, vfcn, '-', 'linewidth', 3); hold('on');
plot(kgrid, v_true, '--', 'linewidth', 3); hold('off');
%title('���l�֐�', 'fontsize', 16);
xlabel('���{�ۗL�ʁFk', 'Fontsize', 16);
ylabel('���l�֐��FV(k)', 'Fontsize', 16);
xlim([0,kmax]);
legend('�ߎ���', '��͓I��', 'Location', 'SouthEast');
grid on;
set(gca,'Fontsize', 16);
saveas (gcf,'Fig3_dndp1.eps','epsc2');

figure;
plot(kgrid, pfcn, '-', 'linewidth', 3); hold('on');
plot(kgrid, p_true, '--', 'linewidth', 3);
plot(kgrid, kgrid, ':', 'linewidth', 2); hold('off');
%title('����֐�', 'fontsize', 16);
xlabel('�����̎��{�ۗL�ʁFk', 'Fontsize', 16);
ylabel("�����̎��{�ۗL�ʁFk'", 'Fontsize', 16);
xlim([0,kmax]);
legend('�ߎ���', '��͓I��', '45�x��', 'Location', 'NorthWest');
grid on;
set(gca,'Fontsize', 16);
saveas (gcf,'Fig3_dndp2.eps','epsc2');

[r, l] = size(v_conv);
iter = linspace(1, l, l)';

figure;
plot(iter, v_conv, '-', 'linewidth', 3);
%title('���l�֐��̎���', 'fontsize', 16);
xlabel('�v�Z��', 'Fontsize', 16);
ylabel('�J��Ԃ��v�Z�덷', 'Fontsize', 16);
grid on;
set(gca,'Fontsize', 16);
saveas (gcf,'Fig3_dndp3.eps','epsc2');

figure;
plot(iter, p_conv, '-', 'linewidth', 3);
%title('����֐��̎���', 'fontsize', 16);
xlabel('�v�Z��', 'Fontsize', 16);
ylabel('�J��Ԃ��v�Z�덷', 'Fontsize', 16);
grid on;
set(gca,'Fontsize', 16);
saveas (gcf,'Fig3_dndp4.eps','epsc2');

figure;
plot(iter, v_conv, '-', 'linewidth', 2); hold('on');
plot(iter, p_conv, ':', 'linewidth', 2); hold('off');
%title('���l�֐��E����֐��̎���', 'fontsize', 16);
xlabel('�v�Z��', 'Fontsize', 16);
ylabel('�J��Ԃ��v�Z�덷', 'Fontsize', 16);
ylim([0,0.1]);
legend('���l�֐�', '����֐�', 'Location', 'NorthEast');
grid on;
set(gca,'Fontsize', 16);
saveas (gcf,'Fig3_dndp5.eps','epsc2');

figure;
plot(kgrid, val_tmp(:, 1), '-', 'linewidth', 3); hold('on');
plot(kgrid, val_tmp(:, 2), '--', 'linewidth', 3);
plot(kgrid, val_tmp(:, 3), '-.', 'linewidth', 3);
plot(kgrid, val_tmp(:, 4), ':', 'linewidth', 3); hold('off');
title('���l�֐��̎���', 'fontsize', 16);
xlabel('���{�ۗL�ʁFk', 'Fontsize', 16);
ylabel('���l�֐��FV(k)', 'Fontsize', 16);
xlim([0,kmax]);
legend('it=1', 'it=3', 'it=5', '����', 'Location', 'East');
grid on;
set(gca,'Fontsize', 16);
saveas (gcf,'Fig3_dndp6.eps','epsc2');



%% �����̐}

figure;
plot(kgrid, vfcn, '-', 'color', 'black', 'linewidth', 3); hold('on');
plot(kgrid, v_true, '--', 'color', 'black', 'linewidth', 3); hold('off');
%title('���l�֐�', 'fontsize', 16);
xlabel('���{�ۗL�ʁFk', 'Fontsize', 16);
ylabel('���l�֐��FV(k)', 'Fontsize', 16);
xlim([0,kmax]);
legend('�ߎ���', '��͓I��', 'Location', 'SouthEast');
grid on;
set(gca,'Fontsize', 16);
saveas (gcf,'Fig3_dndp1_bk.eps','epsc2');

figure;
plot(kgrid, pfcn, '-', 'color', 'black', 'linewidth', 3); hold('on');
plot(kgrid, p_true, '--', 'color', 'black', 'linewidth', 3);
plot(kgrid, kgrid, ':', 'color', 'black', 'linewidth', 2); hold('off');
%title('����֐�', 'fontsize', 16);
xlabel('�����̎��{�ۗL�ʁFk', 'Fontsize', 16);
ylabel("�����̎��{�ۗL�ʁFk'", 'Fontsize', 16);
xlim([0,kmax]);
legend('�ߎ���', '��͓I��', '45�x��', 'Location', 'NorthWest');
grid on;
set(gca,'Fontsize', 16);
saveas (gcf,'Fig3_dndp2_bk.eps','epsc2');

[r, l] = size(v_conv);
iter = linspace(1, l, l)';

figure;
plot(iter, v_conv, '-', 'color', 'black', 'linewidth', 3);
%title('���l�֐��̎���', 'fontsize', 16);
xlabel('�v�Z��', 'Fontsize', 16);
ylabel('�J��Ԃ��v�Z�덷', 'Fontsize', 16);
grid on;
set(gca,'Fontsize', 16);
saveas (gcf,'Fig3_dndp3_bk.eps','epsc2');

figure;
plot(iter, p_conv, '-', 'color', 'black', 'linewidth', 3);
%title('����֐��̎���', 'fontsize', 16);
xlabel('�v�Z��', 'Fontsize', 16);
ylabel('�J��Ԃ��v�Z�덷', 'Fontsize', 16);
grid on;
set(gca,'Fontsize', 16);
saveas (gcf,'Fig3_dndp4_bk.eps','epsc2');

figure;
plot(iter, v_conv, '-', 'color', 'black', 'linewidth', 2); hold('on');
plot(iter, p_conv, ':', 'color', 'black', 'linewidth', 2); hold('off');
%title('���l�֐��E����֐��̎���', 'fontsize', 16);
xlabel('�v�Z��', 'Fontsize', 16);
ylabel('�J��Ԃ��v�Z�덷', 'Fontsize', 16);
ylim([0,0.1]);
legend('���l�֐�', '����֐�', 'Location', 'NorthEast');
grid on;
set(gca,'Fontsize', 16);
saveas (gcf,'Fig3_dndp5_bk.eps','epsc2');

figure;
plot(kgrid, val_tmp(:, 1), '-', 'color', 'black', 'linewidth', 3); hold('on');
plot(kgrid, val_tmp(:, 2), '--', 'color', 'black', 'linewidth', 3);
plot(kgrid, val_tmp(:, 3), '-.', 'color', 'black', 'linewidth', 3);
plot(kgrid, val_tmp(:, 4), ':', 'color', 'black', 'linewidth', 3); hold('off');
title('���l�֐��̎���', 'fontsize', 16);
xlabel('���{�ۗL�ʁFk', 'Fontsize', 16);
ylabel('���l�֐��FV(k)', 'Fontsize', 16);
xlim([0,kmax]);
legend('it=1', 'it=3', 'it=5', '����', 'Location', 'East');
grid on;
set(gca,'Fontsize', 16);
saveas (gcf,'Fig3_dndp6_bk.eps','epsc2');

return

