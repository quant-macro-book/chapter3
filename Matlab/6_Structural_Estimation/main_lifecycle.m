%% ���C���t�@�C��:
% ��ԕϐ��̂ݗ��U�����đ���ϐ��͘A���I�ɒl�����ꍇ�̓��I�v��@(parametric DP)�̉�@.
% �A���S���Y���̏ڍׂ́AJohnson et al. (1993)���Q��

clear;
clear global;
close all;
format short;

global beta gamma asset ss r tran eta endow surv ret_age age zt na nz vfcn_yng vfcn_old agrid

%% *** �J���u���[�V���� ***
beta  = 0.98; % �������q
gamma = 1.0;  % ���ΓI�댯���x(�َ��_�Ԃ̑�ւ̒e�͐��̋t��)
r     = 0.04; % ���q��
ss    = 0.5;  % ������֗��ŃJ���u���[�g

eta = readmatrix("earnings_profiles.csv");
surv = readmatrix("surv.csv");

nz = 3;
endow = [0.8027, 1.0, 1.2457];
tran = [0.7451 0.2528 0.0021; 0.1360 0.7281 0.1360; 0.0021 0.2528 0.7451];

max_age = 86; %20�΂���105�΂܂�
ret_age = 45; %20�΂���64�΂܂�

% *** ���U���p�̃p�����[�^ ***
na   = 101;  % �O���b�h�̐�
amax = 40.0; % ���Y�O���b�h�̍ő�l
amin = 0.0;  % ���Y�O���b�h�̍ŏ��l 
%========================

%% �v�Z�J�n

tic

disp('')
disp('-+- Solve a life cycle model -+-');

%% STEP 1(a): �O���b�h����

agrid = linspace(amin, amax, na)';

%% STEP 1(b): ���l�֐��E����֐������邽�߂̕ϐ�

pfcn_yng = zeros(na, nz, ret_age);
vfcn_yng = zeros(na, nz, ret_age);

pfcn_old = zeros(na, max_age - ret_age);
vfcn_old = zeros(na, max_age - ret_age);

%% STEP 4: ���l�֐����������Ɍv�Z

% �ŏI���͂��ׂĂ̎������g���؂�
for i = 1:na
    vfcn_old(i, max_age-ret_age) = CRRA((1+r)*agrid(i) + ss, gamma);
    pfcn_old(i, max_age-ret_age) = 0.0;
end

figure;
plot(agrid, vfcn_old(:, max_age-ret_age), '-', 'linewidth', 3);
xlabel('���~', 'Fontsize', 16);
ylabel('���l�֐�', 'Fontsize', 16);
grid on;

% ���ތ�̉��l�֐�
for t = max_age-1:-1:ret_age+1

    fprintf('age index: %i \n', t);

    for i = 1:na
    
        % �O���[�o���ϐ���ݒ�
        % fminsearch�Ŏg���֐�(BellmanEq)�ɍœK������ϐ�"�ȊO"�̕ϐ���n��(�O���[�o���ϐ����g��Ȃ����@������͂�)
        asset = agrid(i);
        age = t;
    
        % MATLAB�̍œK���֐�(fminsearch)���g���ăO���b�h��ŉ��l�֐��Ɛ���֐��̒l��T��
        % �����l��0.01
        [pfcn_old(i, age-ret_age), vfcn_old(i, age-ret_age)] = fminsearch(@bellman_eq_retired, 0.01);
    
    end

    % fminsearch�͍ŏ��l��T���֐��Ȃ̂ŕ����𔽓]������
    vfcn_old(:, age-ret_age) = -1*vfcn_old(:, age-ret_age);

end

figure;
plot(agrid, vfcn_old(:, 1), '-', 'linewidth', 3);
xlabel('���~', 'Fontsize', 16);
ylabel('���l�֐�', 'Fontsize', 16);
grid on;

figure;
plot(agrid, pfcn_old(:, 1), '-', 'linewidth', 3);
xlabel('���~', 'Fontsize', 16);
ylabel('����֐�', 'Fontsize', 16);
grid on;

% �J���҂̉��l�֐�
for t = ret_age:-1:1

    fprintf('age index: %i \n', t);

    for z = 1:nz
        for i = 1:na
        
            % �O���[�o���ϐ���ݒ�
            % fminsearch�Ŏg���֐�(BellmanEq)�ɍœK������ϐ�"�ȊO"�̕ϐ���n��(�O���[�o���ϐ����g��Ȃ����@������͂�)
            asset = agrid(i);
            age = t;
            zt  = z;

            % MATLAB�̍œK���֐�(fminsearch)���g���ăO���b�h��ŉ��l�֐��Ɛ���֐��̒l��T��
            % �����l��0.01
            [pfcn_yng(i, z, age), vfcn_yng(i, z, age)] = fminsearch(@bellman_eq_worker, 0.01);
        
        end
    end

    % fminsearch�͍ŏ��l��T���֐��Ȃ̂ŕ����𔽓]������
    vfcn_yng(:, :, age) = -1*vfcn_yng(:, :, age);

end

toc

%% �}��`��

age = linspace(20, 64, 45);

figure;
plot(age, eta, '-', 'linewidth', 3);
xlabel('�N��', 'Fontsize', 16);
ylabel('�N��Ƃ̘J�����Y��', 'Fontsize', 16);
xlim([20,65]);
ylim([0,1.5]);
grid on;
set(gca,'Fontsize', 16);
saveas (gcf,'Fig3_eta.eps','epsc2');

age = linspace(20, 105, 86);

figure;
plot(age, surv, '-', 'linewidth', 3);
xlabel('�N��', 'Fontsize', 16);
ylabel('�����m��', 'Fontsize', 16);
xlim([20,105]);
ylim([0,1.0]);
grid on;
set(gca,'Fontsize', 16);
saveas (gcf,'Fig3_surv.eps','epsc2');

figure;
plot(agrid, vfcn_yng(:, 1, 1), '-', 'linewidth', 3); hold('on');
plot(agrid, vfcn_yng(:, 2, 1), '-', 'linewidth', 3);
plot(agrid, vfcn_yng(:, 3, 1), '-', 'linewidth', 3); hold('off');
xlabel('���~', 'Fontsize', 16);
ylabel('���l�֐�', 'Fontsize', 16);
grid on;

figure;
plot(agrid, pfcn_old(:, 1, 1), '-', 'linewidth', 3); hold('on');
plot(agrid, pfcn_old(:, 2, 1), '-', 'linewidth', 3);
plot(agrid, pfcn_old(:, 3, 1), '-', 'linewidth', 3); hold('off');
xlabel('���~', 'Fontsize', 16);
ylabel('����֐�', 'Fontsize', 16);
grid on;

return
