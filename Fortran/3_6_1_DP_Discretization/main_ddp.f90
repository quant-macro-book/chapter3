program main_ddp
! ���C���t�@�C��:
! ��ԕϐ��Ƒ���ϐ��𗣎U�����ē��I�v��@(discretized DP)������.

    use mod_types, only:dp
    use mod_my_econ_fcn, only:CRRA
    use mod_generate_grid, only:grid_uniform

    implicit none

    ! *** �J���u���[�V���� ***
    real(dp), parameter :: beta  = 0.96 ! �������q
    real(dp), parameter :: gamma = 1.0  ! ���ΓI�댯���x(�َ��_�Ԃ̑�ւ̒e�͐��̋t��)
    real(dp), parameter :: alpha = 0.40 ! ���{���z��
    real(dp), parameter :: delta = 1.00 ! �Œ莑�{����(0.08)

    ! *** ���U���p�̃p�����[�^ ***
    integer,  parameter :: nk = 1001    ! �O���b�h�̐�
    real(dp), parameter :: kmax = 0.5  ! ���{�O���b�h�̍ő�l
    real(dp), parameter :: kmin = 0.05 ! ���{�O���b�h�̍ŏ��l (0�ɂ���Ɛ��Y���o���Ȃ��Ȃ�)
    real(dp), dimension(nk) :: kgrid = 0.0

    ! *** ���[�J���ϐ� ***
    integer :: i, j
    real(dp) :: cons, wealth
    real(dp) :: time_begin, time_end
    real(dp), dimension(nk) :: vfcn, pfcn
    real(dp), dimension(nk) :: Tvfcn, Tpfcn
    real(dp), dimension(nk, nk) :: vkp
    real(dp), dimension(nk, nk) :: util

    ! *** ��͓I�� ***
    real(dp) :: AA, BB
    real(dp), dimension(nk) :: v_true
    real(dp), dimension(nk) :: p_true
    
    ! *** �����̊ ***
    integer :: it = 1 ! ���[�v�E�J�E���^�[
    integer,  parameter :: maxit = 1000   ! �J��Ԃ��v�Z�̍ő�l
    real(dp), parameter :: tol = 1.0e-005 ! ���e�덷(STEP 2)
    real(dp) :: dif1 = 1.0 ! ���l�֐��̌J��Ԃ��v�Z�덷
    real(dp) :: dif2 = 1.0 ! ����֐��̌J��Ԃ��v�Z�덷
    !----------------------------------------------

    ! �v�Z���Ԃ��J�E���g�J�n
    call cpu_time( time_begin )

    write (*,*) ""
    write (*,*) "-+- Solve a neoclassical growth model -+-"
    write (*,*) ""

    ! STEP 1(a): �O���b�h����
    call grid_uniform(kmin, kmax, nk, kgrid)

    ! STEP 1(b): ���l�֐��E����֐��̏����l��ݒ�
    vfcn  = 0.0
    pfcn  = 0.0
    Tvfcn = 0.0
    Tpfcn = 0.0
    vkp   = 0.0

    ! STEP 3: ���p�֐��̑g�ݍ��킹

    ! ���p�֐��̏����l (���0�ȉ��ɂȂ�g�ݍ��킹�ɂ̓y�i���e�B)
    util = -10000.0

    ! ������l�ɂȂ�(k,k')�̑g�ݍ��킹�ɂ��Č��p���v�Z
    do i = 1,nk
        !  �����鑀��ϐ�k'�ɂ���:
        do j = 1,nk
            wealth = kgrid(i)**alpha + (1.0-delta)*kgrid(i)
            cons = wealth - kgrid(j)
            if (cons > 0) then
               util(j,i) = CRRA(cons, gamma)
            end if
        end do
    end do

    ! STEP 4: ���l�֐����J��Ԃ��v�Z
    do it = 1,maxit

        !if (dif1 < tol .or. dif2 < tol) exit
        if (dif1 < tol) exit

        ! �x���}��������: V(k;k')
        do i = 1,nk
            vkp(:,i) = util(:,i) + beta*vfcn
        end do

        ! �œK��: �ek�ɂ���V(k;k')���ő�ɂ���k'��T��
        do i = 1,nk
            Tvfcn(i) = maxval(vkp(:,i))
            Tpfcn(i) = kgrid(maxloc(vkp(:,i),dim=1))
        end do

        ! �J��Ԃ��v�Z�덷���m�F
        dif1 = maxval(abs((Tvfcn-vfcn)/vfcn))
        dif2 = maxval(abs((Tpfcn-pfcn)/pfcn))

        ! ���l�֐��E����֐����A�b�v�f�[�g
        vfcn = Tvfcn
        pfcn = Tpfcn
        write (*,"('iteration index: ', i3, ', iteration diff of value: ', f9.5, ', iteration diff of policy: ', f9.5)") it, dif1, dif2

    end do

    ! �v�Z���Ԃ��J�E���g�I��
    call cpu_time( time_end )
    write (*,*) ""
    write (*,"(' Program finished sucessfully in', f12.9, ' seconds')") time_end - time_begin

    ! �v�Z���ʂ��R�}���h�v�����v�g�ɕ\��
    write (*,*) ""
    write (*,*) "-+- Parameter values -+-"
    write (*,*) ""
    write (*,"('beta=', f5.2, ', gamma=', f5.2, ', alpha=', f5.2, ',delta=', f5.2)") beta, gamma, alpha, delta
    write (*,*) ""
    write (*,"('kmin=', f5.2, ', kmax=', f5.2, ', #grid', i5)") kmin, kmax, nk
    write (*,*) ""

    ! ��͓I��
    AA = (1.0-beta)**(-1) * (log(1.0-alpha*beta) + ((alpha*beta)/(1.0-alpha*beta))*log(alpha*beta))
    BB = alpha/(1.0-alpha*beta)
    v_true = AA + BB*log(kgrid)
    p_true = beta*alpha*(kgrid**alpha)

    ! �O���t�@�C���ɏo��
    open (10, file='result.csv', status='replace')
    do i = 1,nk
        write (10,10) kgrid(i), ',', vfcn(i), ',', pfcn(i), ',', v_true(i), ',', p_true(i)
        10 format(f0.6,a1,f0.6,a1,f0.6,a1,f0.6,a1,f0.6)
    end do
    close (10)

end program