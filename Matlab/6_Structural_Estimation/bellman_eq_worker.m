function  value = bellman_eq_worker(aprime)
% Function bellman_eq_worker
%  [value] = bellman_eq_worker(kprime)
%
% �ړI:
% a'����^�����Ƃ��̃x���}����������Ԃ��֐�.
% main_lifecycle.m����Ăяo���Ďg��.
%
% �O���[�o���ϐ�: beta gamma alpha delta A tran capital vfcn kgrid

global beta gamma asset ss r tran eta endow surv ret_age age zt na nz vfcn_yng vfcn_old agrid

%% �x���}��������

% �\�Z����Ə���
wealth = (1+r)*asset + eta(age)*endow(zt);
cons = wealth - aprime;

% ������l�̏ꍇ�A�y�i���e�B��^���Ă��̒l���I�΂�Ȃ��悤�ɂ���
if cons > 0.0
    util = CRRA(cons, gamma);
else
    util = -10000.0;
end

% �����̊��Ҍ��p
if age == ret_age
    % ���ޒ��O�ł���Ύ����ɂ�z�̕s�m�����͑��݂��Ȃ��̂Ŋ��Ғl�v�Z�͕s�v
    vnext = interp1(agrid, vfcn_old(:, 1), aprime, 'spline');
    vnext = surv(age)*beta*vnext;
else
    vnext_z = zeros(nz, 1);
    % ���ꂼ��̎�����z�ɑΉ�����value���v�Z
    for zz = 1:nz
        if aprime <= agrid(na)
            vnext_z(zz) = interp1(agrid, vfcn_yng(:, zz, age+1), aprime, 'spline');
        else
            vnext_z(zz) = interp1(agrid, vfcn_yng(:, zz, age+1), aprime, 'linear');
        end
    end
    % �J�ڊm���s��Ŋ��Ғl���v�Z
    vnext = surv(age)*beta*dot(tran(zt, :), vnext_z);
end

% ���l�֐�
value = util + vnext;

%% �g���b�N(1): k'�͐��̒l�������Ȃ��̂ŁA�y�i���e�B��^���Ă��̒l���I�΂�Ȃ��悤�ɂ���
if aprime < 0
    value = -1000000.0;
end

%% �g���b�N(2): "�ŏ���"������̂ŕ����𔽓]
value = -1.0 * value;
 
return
