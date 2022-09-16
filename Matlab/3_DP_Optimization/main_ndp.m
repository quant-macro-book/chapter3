%% ���C���t�@�C��:
% ��ԕϐ��̂ݗ��U�����đ���ϐ��͘A���I�ɒl�����ꍇ�̓��I�v��@(parametric DP)�̉�@.
% �A���S���Y���̏ڍׂ́AJohnson et al. (1993)���Q��

clear;
clear global;
close all;
format short;

global beta gamma alpha delta capital vfcn kgrid

%% *** �J���u���[�V���� ***
beta  = 0.96; % �������q
gamma = 1.0;  % ���ΓI�댯���x(�َ��_�Ԃ̑�ւ̒e�͐��̋t��)
alpha = 0.40; % ���{���z��
delta = 1.00; % �Œ莑�{����(0.08)

% *** ���U���p�̃p�����[�^ ***
nk   = 21;    % �O���b�h�̐�
kmax = 0.5;   % ���{�O���b�h�̍ő�l
%kmax = 10.0; % ���{�O���b�h�̍ő�l(�Œ莑�{����=0.08�̏ꍇ�Ɏg�p)
kmin = 0.05;  % ���{�O���b�h�̍ŏ��l (0�ɂ���Ɛ��Y���o���Ȃ��Ȃ�)
%========================

% *** �����̊ ***
it = 1;          % ���[�v�E�J�E���^�[
maxit = 1000;    % �J��Ԃ��v�Z�̍ő�l
tol  = 1.0e-005; % ���e�덷(STEP 2)
dif1 = 1.0;      % ���l�֐��̌J��Ԃ��덷
dif2 = 1.0;      % ����֐��̌J��Ԃ��덷
count = 1;
%=================

%% �v�Z�J�n

tic

disp('')
disp('-+- Solve a neoclassical growth model -+-');

%% STEP 1(a): �O���b�h����

kgrid = linspace(kmin, kmax, nk)';
% kgrid = grid_exp1(kmin, kmax, nk)';
% kgrid = grid_exp2(kmin, kmax, nk)';
% kgrid = grid_exp3(kmin, kmax, nk)';

%% STEP 1(b): ���l�֐��E����֐��̏����l�𓖂Đ���

pfcn0 = zeros(nk, 1);
vfcn0 = CRRA(kgrid.^alpha + (1.-delta).*kgrid, gamma);

pfcn1 = zeros(nk, 1);
vfcn1 = zeros(nk, 1);

% ���l�֐��E����֐��̌o�H���L�^(�Ȃ��Ă���)
vpath(:, 1) = vfcn0;
ppath(:, 1) = pfcn0;

%% STEP 4: ���l�֐����J��Ԃ��v�Z

while it < maxit && dif1 > tol

    fprintf('iteration index: %i \n', it);
    fprintf('value function iteration error: %e\n', dif1);
    fprintf('policy function iteration error: %e\n', dif2);

    for i = 1:nk

        % �O���[�o���ϐ���ݒ�
        % fminsearch�Ŏg���֐�(BellmanEq)�ɍœK������ϐ�"�ȊO"�̕ϐ���n��(�O���[�o���ϐ����g��Ȃ����@������͂�)
        capital = kgrid(i);
        vfcn = vfcn0;

        % MATLAB�̍œK���֐�(fminsearch)���g���ăO���b�h��ŉ��l�֐��Ɛ���֐��̒l��T��
        % �����l��0.01
        [pfcn1(i,1), vfcn1(i,1)] = fminsearch(@BellmanEq, 0.01);

    end

    % fminsearch�͍ŏ��l��T���֐��Ȃ̂ŕ����𔽓]������
    vfcn1 = -1*vfcn1;

    % �J��Ԃ��v�Z�덷���m�F
    dif1 = max(abs((vfcn1-vfcn0)./vfcn0));
    dif2 = max(abs((pfcn1-pfcn0)./pfcn0));

    % �����r���̌J��Ԃ��v�Z�덷��ۑ�
    % �r���o�߂�}������ړI�Ȃ̂ŁA�ʏ�͕s�v(�ނ���x���Ȃ�̂ŏ����ׂ�)
    % �v�Z���ɍs��̃T�C�Y���ς���Ă����͖̂]�܂����Ȃ��������Ȃ̂Ŗ{���͔�����ׂ�
    dif(1, it) = dif1;
    dif(2, it) = dif2;

    % ���l�֐��̌o�H���L�^(�Ȃ��Ă����Ȃ�)
    vpath(:, it) = vfcn0;
    ppath(:, it) = pfcn0;

    % ���l�֐��E����֐����A�b�v�f�[�g
    vfcn0 = vfcn1;
    pfcn0 = pfcn1;

    it = it + 1;

end

% �ŏI�I�Ȑ���֐��������Ă������֐����v�Z
cfcn = kgrid.^alpha + (1.-delta).*kgrid - pfcn0(:,1);

%% �v�Z���ʂ��R�}���h�E�B���h�E�ɕ\��

disp('-+- PARAMETER VALUES -+-');
disp('');
fprintf('beta=%5.2f, gamma=%5.2f, alpha=%5.2f, delta=%5.2f \n', beta, gamma, alpha, delta);
disp(''); 
fprintf('kmin=%5.2f, kmax=%5.2f, #grid=%i \n', kmin, kmax, nk);
disp('');

toc

%% ��͓I��

AA = (1.0-beta).^(-1) * (log(1.0-alpha*beta) + ((alpha*beta)/(1.0-alpha*beta))*log(alpha*beta));
BB = alpha/(1.0-alpha*beta);
v_true = AA + BB*log(kgrid);
p_true = beta*alpha*(kgrid.^alpha);

%% �I�C���[����������덷�𑪒�

cons = kgrid.^alpha + (1.-delta).*kgrid - pfcn0(:,1);
LHS  = mu_CRRA(cons, gamma);
kp   = pfcn0(:,1);
kpp  = interp1(kgrid, pfcn0(:,1), kp);
cons = kp.^alpha + (1.-delta).*kp - kpp;
rent = alpha.*kp.^(alpha-1.0) - delta;
RHS  = beta.*(1.+rent).*mu_CRRA(cons, gamma);
err  = RHS./LHS-1.0;

err2 = csvread("err_ddp.csv");

%% �}��`��

figure;
plot(kgrid, vfcn0, '-', 'color', 'black', 'linewidth', 3); hold('on');
plot(kgrid, v_true, '--', 'color', 'black', 'linewidth', 3); hold('off');
%title('���l�֐�', 'fontsize', 16);
xlabel('���{�ۗL�ʁFk', 'Fontsize', 16);
ylabel('���l�֐��FV(k)', 'Fontsize', 16);
xlim([0, kmax]);
legend('�ߎ���', '��͓I��', 'Location', 'SouthEast');
grid on;
set(gca,'Fontsize', 16);
saveas (gcf,'Fig3_pndp1.eps','epsc2');

figure;
plot(kgrid, pfcn0, '-', 'color', 'black', 'linewidth', 3); hold('on');
plot(kgrid, p_true, '--', 'color', 'black', 'linewidth', 3);
plot(kgrid, kgrid, ':', 'color', 'black', 'linewidth', 2); hold('off');
%title('����֐�', 'fontsize', 16);
xlabel('�����̎��{�ۗL�ʁFk', 'Fontsize', 16);
ylabel("�����̎��{�ۗL�ʁFk'", 'Fontsize', 16);
xlim([0, kmax]);
legend('�ߎ���', '��͓I��', '45�x��', 'Location', 'NorthWest');
grid on;
set(gca,'Fontsize', 16);
saveas (gcf,'Fig3_pndp2.eps','epsc2');

figure;
plot(kgrid, cfcn(:,1), '-', 'color', 'black', 'linewidth', 3);
%title('����֐�', 'fontsize', 16);
xlabel('���{�ۗL�ʁFk', 'Fontsize', 16);
ylabel('����Fc', 'Fontsize', 16);
grid on;
set(gca,'Fontsize', 16);
saveas (gcf,'Fig3_pndp3.eps','epsc2');

iter = linspace(1, it-1, it-1);

figure;
plot(iter, dif(1, :), '-', 'color', 'black', 'linewidth', 2); hold('on');
plot(iter, dif(2, :), ':', 'color', 'black', 'linewidth', 2); hold('off');
%title('���l�֐��E����֐��̎���', 'fontsize', 16);
xlabel('�v�Z��', 'Fontsize', 16);
ylabel('�J��Ԃ��v�Z�덷', 'Fontsize', 16);
ylim([0,0.1]);
legend('���l�֐�', '����֐�', 'Location', 'NorthEast');
grid on;
set(gca,'Fontsize', 16);
saveas (gcf,'Fig3_pndp4.eps','epsc2');

figure;
plot(kgrid, err, '-', 'color', 'black', 'linewidth', 3);
%title('�I�C���[�������̌덷', 'fontsize', 16);
xlabel('���{�ۗL�ʁFk', 'Fontsize', 16);
ylabel('�I�C���[�������덷', 'Fontsize', 16);
grid on;
set(gca,'Fontsize', 16);
saveas (gcf,'Fig3_pndp5.eps','epsc2');

figure;
plot(kgrid, err2, '-', 'color', 'black', 'linewidth', 3); hold('on');
plot(kgrid, err, '--', 'color', 'black', 'linewidth', 3); hold('off');
%title('�I�C���[�������̌덷', 'fontsize', 16);
xlabel('���{�ۗL�ʁFk', 'Fontsize', 16);
ylabel('�I�C���[�������덷', 'Fontsize', 16);
ylim([-15e-004,5e-004]);
legend('����ϐ��F���U', '����ϐ��F�A��', 'Location', 'SouthEast');
grid on;
set(gca,'Fontsize', 16);
saveas (gcf,'Fig3_pndp6.eps','epsc2');

return
