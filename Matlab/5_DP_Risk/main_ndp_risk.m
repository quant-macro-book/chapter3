%% ���C���t�@�C��:
% ��ԕϐ��̂ݗ��U�����đ���ϐ��͘A���I�ɒl�����ꍇ�̓��I�v��@(parametric DP)�̉�@.
% �A���S���Y���̏ڍׂ́AJohnson et al. (1993)���Q��

clear;
clear global;
close all;
format short;

global beta gamma alpha delta A tran capital vfcn kgrid

%% *** �J���u���[�V���� ***
beta  = 0.96; % �������q
gamma = 1.0;  % ���ΓI�댯���x(�َ��_�Ԃ̑�ւ̒e�͐��̋t��)
alpha = 0.40; % ���{���z��
delta = 1.00; % �Œ莑�{����(0.08)
tfp   = [1.01, 0.99];
Pi_A  = [1-0.125 0.125; 0.125 1-0.125];

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
disp('-+- Solve a neoclassical growth model with TFP shock -+-');

%% STEP 1(a): �O���b�h����

kgrid = linspace(kmin, kmax, nk)';
% kgrid = grid_exp1(kmin, kmax, nk)';
% kgrid = grid_exp2(kmin, kmax, nk)';
% kgrid = grid_exp3(kmin, kmax, nk)';

%% STEP 1(b): ���l�֐��E����֐��̏����l�𓖂Đ���

% A�̏�Ԃ�good��bad��2��
pfcn0 = zeros(nk, 2);
vfcn0 = zeros(nk, 2);
vfcn0(:, 1) = CRRA(tfp(1)*kgrid.^alpha + (1.-delta).*kgrid, gamma);
vfcn0(:, 2) = CRRA(tfp(2)*kgrid.^alpha + (1.-delta).*kgrid, gamma);

pfcn1 = zeros(nk, 2);
vfcn1 = zeros(nk, 2);

%% STEP 4: ���l�֐����J��Ԃ��v�Z

while it < maxit && dif1(1) > tol

    fprintf('iteration index: %i \n', it);
    fprintf('value function iteration error: %e\n', dif1);
    fprintf('policy function iteration error: %e\n', dif2);

    for z = 1:2
        for i = 1:nk
    
            % �O���[�o���ϐ���ݒ�
            % fminsearch�Ŏg���֐�(BellmanEq)�ɍœK������ϐ�"�ȊO"�̕ϐ���n��(�O���[�o���ϐ����g��Ȃ����@������͂�)
            capital = kgrid(i);
            A = tfp(z);
            tran = Pi_A(z, :);
            vfcn = vfcn0;
    
            % MATLAB�̍œK���֐�(fminsearch)���g���ăO���b�h��ŉ��l�֐��Ɛ���֐��̒l��T��
            % �����l��0.01
            [pfcn1(i,z), vfcn1(i,z)] = fminsearch(@BellmanEq, 0.01);
    
        end
    end

    % fminsearch�͍ŏ��l��T���֐��Ȃ̂ŕ����𔽓]������
    vfcn1 = -1*vfcn1;

    % �J��Ԃ��v�Z�덷���m�F
    dif1 = max(abs((vfcn1-vfcn0)./vfcn0));
    dif2 = max(abs((pfcn1-pfcn0)./pfcn0));

    % ���l�֐��E����֐����A�b�v�f�[�g
    vfcn0 = vfcn1;
    pfcn0 = pfcn1;

    it = it + 1;

end

%% �v�Z���ʂ��R�}���h�E�B���h�E�ɕ\��

disp('-+- PARAMETER VALUES -+-');
disp('');
fprintf('beta=%5.2f, gamma=%5.2f, alpha=%5.2f, delta=%5.2f \n', beta, gamma, alpha, delta);
disp(''); 
fprintf('kmin=%5.2f, kmax=%5.2f, #grid=%i \n', kmin, kmax, nk);
disp('');

toc

%% �}��`��

figure;
plot(kgrid, vfcn0(:, 1), '-', 'linewidth', 3); hold('on');
plot(kgrid, vfcn0(:, 2), '-.', 'linewidth', 3); hold('off');
%title('���l�֐�', 'fontsize', 16);
xlabel('���{�ۗL�ʁFk', 'Fontsize', 16);
ylabel('���l�֐��FV(k)', 'Fontsize', 16);
xlim([0, kmax]);
legend('�D��', '�s��', 'Location', 'NorthWest');
grid on;
set(gca,'Fontsize', 16);
saveas (gcf,'Fig3_pndp_risk1.eps','epsc2');

figure;
plot(kgrid, pfcn0(:, 1), '-', 'linewidth', 3); hold('on');
plot(kgrid, pfcn0(:, 2), '-.', 'linewidth', 3);
plot(kgrid, kgrid, ':', 'linewidth', 2); hold('off');
%title('����֐�', 'fontsize', 16);
xlabel('�����̎��{�ۗL�ʁFk', 'Fontsize', 16);
ylabel("�����̎��{�ۗL�ʁFk'", 'Fontsize', 16);
xlim([0, kmax]);
ylim([0, kmax]);
legend('�D��', '�s��', '45�x��', 'Location', 'NorthWest');
grid on;
set(gca,'Fontsize', 16);
saveas (gcf,'Fig3_pndp_risk2.eps','epsc2');


%% �����̐}

figure;
plot(kgrid, vfcn0(:, 1), '-', 'color', 'black', 'linewidth', 3); hold('on');
plot(kgrid, vfcn0(:, 2), '-.', 'color', 'black', 'linewidth', 3); hold('off');
%title('���l�֐�', 'fontsize', 16);
xlabel('���{�ۗL�ʁFk', 'Fontsize', 16);
ylabel('���l�֐��FV(k)', 'Fontsize', 16);
xlim([0, kmax]);
legend('�D��', '�s��', 'Location', 'NorthWest');
grid on;
set(gca,'Fontsize', 16);
saveas (gcf,'Fig3_pndp_risk1_bk.eps','epsc2');

figure;
plot(kgrid, pfcn0(:, 1), '-', 'color', 'black', 'linewidth', 3); hold('on');
plot(kgrid, pfcn0(:, 2), '-.', 'color', 'black', 'linewidth', 3);
plot(kgrid, kgrid, ':', 'color', 'black', 'linewidth', 2); hold('off');
%title('����֐�', 'fontsize', 16);
xlabel('�����̎��{�ۗL�ʁFk', 'Fontsize', 16);
ylabel("�����̎��{�ۗL�ʁFk'", 'Fontsize', 16);
xlim([0, kmax]);
ylim([0, kmax]);
legend('�D��', '�s��', '45�x��', 'Location', 'NorthWest');
grid on;
set(gca,'Fontsize', 16);
saveas (gcf,'Fig3_pndp_risk2_bk.eps','epsc2');

return
