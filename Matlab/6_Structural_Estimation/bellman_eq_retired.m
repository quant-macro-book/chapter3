function  value = bellman_eq_retired(aprime)
% Function bellman_eq_retired
%  [value] = bellman_eq_retired(kprime)
%
% �ړI:
% a'����^�����Ƃ��̃x���}����������Ԃ��֐�.
% main_lifecycle.m����Ăяo���Ďg��.
%
% �O���[�o���ϐ�: beta gamma alpha delta A tran capital vfcn kgrid

global beta gamma asset ss r surv ret_age age vfcn_old agrid

%% �x���}��������

% �\�Z����Ə���
wealth = (1+r)*asset + ss;
cons = wealth - aprime;

% ������l�̏ꍇ�A�y�i���e�B��^���Ă��̒l���I�΂�Ȃ��悤�ɂ���
if cons > 0.0
    util = CRRA(cons, gamma);
else
    util = -10000.0;
end

% �����̉��l�֐����X�v���C�����
vnext = interp1(agrid, vfcn_old(:, age-ret_age+1), aprime, 'spline');

% ���l�֐�
value = util + surv(age)*beta*vnext;

%% �g���b�N(1): k'�͐��̒l�������Ȃ��̂ŁA�y�i���e�B��^���Ă��̒l���I�΂�Ȃ��悤�ɂ���
if aprime < 0
    value = -1000000.0;
end

%% �g���b�N(2): "�ŏ���"������̂ŕ����𔽓]
value = -1.0 * value;
 
return
