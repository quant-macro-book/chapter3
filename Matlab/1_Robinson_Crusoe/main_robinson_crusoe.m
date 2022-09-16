%% ���C���t�@�C��:
% �œK��(optimization)�Ɠ��}�@(interpolation)�������ă��r���\���E�N���[�\�[�o�ς�����.

clear;
clear global;
close all;
format short;

global beta gamma alpha capital val_fcn kgrid

%% *** �J���u���[�V���� ***
beta  = 0.96; % �������q
gamma = 1.0;  % ���ΓI�댯���x(�َ��_�Ԃ̑�ւ̒e�͐��̋t��)
alpha = 0.4;  % ���{���z��

% *** ���U���p�̃p�����[�^ ***
nk   = 11;    % �O���b�h�̐�
kmax =  1.0;  % ���{�O���b�h�̍ő�l
kmin =  0.05; % ���{�O���b�h�̍ŏ��l
%==========================

% *** �ϐ����` ***
TT   = 10;            % ���l���ɑ؍݂������
vfcn = zeros(nk, TT); % ���l�֐�
pfcn = zeros(nk, TT); % ����֐�
cfcn = zeros(nk, TT); % ����֐�
%=================

%% �v�Z�J�n

tic

disp(' ');
fprintf('Solve Robinson Crusoe economy with %d periods \n', TT);

%% �O���b�h�|�C���g���v�Z

kgrid = linspace(kmin, kmax, nk)';

%% �ŏI��(�S�Ă�����)

pfcn(:, TT) = 0; % �S�ď����̂Œ��~�̓[��
cfcn(:, TT) = kgrid.^alpha; % ���Y=����
vfcn(:, TT) = CRRA(cfcn(:, TT), gamma); % ����瓾������p

%% ���C�����[�v

for t = TT-1:-1:1

    fprintf('period %d: \n', t);

    for i = 1:nk

        % �O���[�o���ϐ���ݒ�
        % fminsearch�Ŏg���֐�(BellmanEq)�ɍœK������ϐ�"�ȊO"�̕ϐ���n��(�O���[�o���ϐ����g��Ȃ����@������͂�)
        capital = kgrid(i);
        val_fcn = vfcn(:, t+1);

        % MATLAB�̍œK���֐�(fminsearch)���g���ăO���b�h��ŉ��l�֐��Ɛ���֐��̒l��T��
        % �����l��0.01
        [pfcn(i, t), vfcn(i, t)] = fminsearch(@BellmanEq, 0.01);

    end

    % ����֐����v�Z(�K�������v�Z����K�v�͂Ȃ�)
    cfcn(:, t) = kgrid.^alpha - pfcn(:, t);

    % fminsearch�͍ŏ��l��T���֐��Ȃ̂ŕ����𔽓]������
    vfcn(:, t) = -1*vfcn(:, t);

end

%% �v�Z���ʂ��R�}���h�E�B���h�E�ɕ\��

disp('-+- PARAMETER VALUES -+-');
disp('');
fprintf('beta=%5.2f, gamma=%5.2f, alpha=%5.2f \n', beta, gamma, alpha);
disp(''); 
fprintf('kmin=%5.2f, kmax=%5.2f, #grid=%5.2f \n', kmin, kmax, nk);
disp('');

toc

%% ��͓I��

p_true = zeros(nk, TT);

for t = 1:TT
    for i = 1:nk
        p_true(i, t) = alpha*beta*( (1-(alpha*beta)^(TT-t)) / (1-(alpha*beta)^(TT-t+1)) )*(kgrid(i)^alpha);
    end
end

%% �}��`��

figure;
plot(kgrid, vfcn(:,TT), '-', 'color', 'black', 'linewidth', 3); hold('on');
plot(kgrid, vfcn(:,9), '-.', 'color', 'black', 'linewidth', 3);
plot(kgrid, vfcn(:,8), '--', 'color', 'black', 'linewidth', 3);
plot(kgrid, vfcn(:,7), ':', 'color', 'black', 'linewidth', 3);
plot(kgrid, vfcn(:,1), '-', 'color', 'black', 'linewidth', 3); hold('off');
%title('���l�֐�', 'fontsize', 16);
xlabel('���{�ۗL�ʁFk_{t}', 'Fontsize', 16);
ylabel('���l�֐��FV_{t}(k_{t})', 'Fontsize', 16);
legend('10��', '9��', '8��', '7��', '1��', 'Location', 'SouthEast');
grid on;
set(gca,'Fontsize',16);
saveas (gcf,'Fig3_rc1.eps','epsc2');

figure;
plot(kgrid, pfcn(:,TT), '-', 'color', 'black', 'linewidth', 3); hold('on');
plot(kgrid, pfcn(:,9), '-.', 'color', 'black', 'linewidth', 3);
plot(kgrid, pfcn(:,8), '--', 'color', 'black', 'linewidth', 3);
plot(kgrid, pfcn(:,7), ':', 'color', 'black', 'linewidth', 3);
plot(kgrid, pfcn(:,1), '-', 'color', 'black', 'linewidth', 3); hold('off');
%title('����֐�', 'fontsize', 16);
xlabel('t���̎��{�ۗL�ʁFk_{t}', 'Fontsize', 16);
ylabel("t+1���̎��{�ۗL�ʁFk_{t+1}", 'Fontsize', 16);
legend('10��', '9��', '8��', '7��', '1��', 'Location', 'SouthEast');
grid on;
set(gca,'Fontsize',16);
saveas (gcf,'Fig3_rc2.eps','epsc2');

figure;
plot(kgrid, cfcn(:,TT), '-', 'color', 'black', 'linewidth', 3); hold('on');
plot(kgrid, cfcn(:,9), '-.', 'color', 'black', 'linewidth', 3);
plot(kgrid, cfcn(:,8), '--', 'color', 'black', 'linewidth', 3);
plot(kgrid, cfcn(:,7), ':', 'color', 'black', 'linewidth', 3);
plot(kgrid, cfcn(:,1), '-', 'color', 'black', 'linewidth', 3); hold('off');
%title('����֐�', 'fontsize', 16);
xlabel('���{�ۗL�ʁFk_{t}', 'Fontsize', 16);
ylabel('����Fc_{t}', 'Fontsize', 16);
legend('10��', '9��', '8��', '7��', '1��', 'Location', 'SouthEast');
grid on;
set(gca,'Fontsize',16);
saveas (gcf,'Fig3_rc3.eps','epsc2');

figure;
plot(kgrid, pfcn(:,9), '-', 'color', 'black', 'linewidth', 3); hold('on');
plot(kgrid, p_true(:,9), '-.', 'color', 'black', 'linewidth', 3);
plot(kgrid, pfcn(:,1), ':', 'color', 'black', 'linewidth', 3);
plot(kgrid, p_true(:,1), '-', 'color', 'black', 'linewidth', 3); hold('off');
%title('�^�̐���֐��Ƌߎ���', 'fontsize', 16);
xlabel('t���̎��{�ۗL�ʁFk_{t}', 'Fontsize', 16);
ylabel('t+1���̎��{�ۗL�ʁFk_{t+1}', 'Fontsize', 16);
legend('9���̋ߎ���', '9���̉�͓I��', '1���̋ߎ���', '1���̉�͓I��', 'Location', 'SouthEast');
grid on;
set(gca,'Fontsize',16);
saveas (gcf,'Fig3_rc4.eps','epsc2');

return
