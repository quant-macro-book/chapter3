program main_ddp
! メインファイル:
! 状態変数と操作変数を離散化して動的計画法(discretized DP)を解く.

    use mod_types, only:dp
    use mod_my_econ_fcn, only:CRRA
    use mod_generate_grid, only:grid_uniform

    implicit none

    ! *** カリブレーション ***
    real(dp), parameter :: beta  = 0.96 ! 割引因子
    real(dp), parameter :: gamma = 1.0  ! 相対的危険回避度(異時点間の代替の弾力性の逆数)
    real(dp), parameter :: alpha = 0.40 ! 資本分配率
    real(dp), parameter :: delta = 1.00 ! 固定資本減耗(0.08)

    ! *** 離散化用のパラメータ ***
    integer,  parameter :: nk = 1001    ! グリッドの数
    real(dp), parameter :: kmax = 0.5  ! 資本グリッドの最大値
    real(dp), parameter :: kmin = 0.05 ! 資本グリッドの最小値 (0にすると生産が出来なくなる)
    real(dp), dimension(nk) :: kgrid = 0.0

    ! *** ローカル変数 ***
    integer :: i, j
    real(dp) :: cons, wealth
    real(dp) :: time_begin, time_end
    real(dp), dimension(nk) :: vfcn, pfcn
    real(dp), dimension(nk) :: Tvfcn, Tpfcn
    real(dp), dimension(nk, nk) :: vkp
    real(dp), dimension(nk, nk) :: util

    ! *** 解析的解 ***
    real(dp) :: AA, BB
    real(dp), dimension(nk) :: v_true
    real(dp), dimension(nk) :: p_true
    
    ! *** 収束の基準 ***
    integer :: it = 1 ! ループ・カウンター
    integer,  parameter :: maxit = 1000   ! 繰り返し計算の最大値
    real(dp), parameter :: tol = 1.0e-005 ! 許容誤差(STEP 2)
    real(dp) :: dif1 = 1.0 ! 価値関数の繰り返し計算誤差
    real(dp) :: dif2 = 1.0 ! 政策関数の繰り返し計算誤差
    !----------------------------------------------

    ! 計算時間をカウント開始
    call cpu_time( time_begin )

    write (*,*) ""
    write (*,*) "-+- Solve a neoclassical growth model -+-"
    write (*,*) ""

    ! STEP 1(a): グリッド生成
    call grid_uniform(kmin, kmax, nk, kgrid)

    ! STEP 1(b): 価値関数・政策関数の初期値を設定
    vfcn  = 0.0
    pfcn  = 0.0
    Tvfcn = 0.0
    Tpfcn = 0.0
    vkp   = 0.0

    ! STEP 3: 効用関数の組み合わせ

    ! 効用関数の初期値 (消費が0以下になる組み合わせにはペナルティ)
    util = -10000.0

    ! 消費が正値になる(k,k')の組み合わせについて効用を計算
    do i = 1,nk
        !  あらゆる操作変数k'について:
        do j = 1,nk
            wealth = kgrid(i)**alpha + (1.0-delta)*kgrid(i)
            cons = wealth - kgrid(j)
            if (cons > 0) then
               util(j,i) = CRRA(cons, gamma)
            end if
        end do
    end do

    ! STEP 4: 価値関数を繰り返し計算
    do it = 1,maxit

        !if (dif1 < tol .or. dif2 < tol) exit
        if (dif1 < tol) exit

        ! ベルマン方程式: V(k;k')
        do i = 1,nk
            vkp(:,i) = util(:,i) + beta*vfcn
        end do

        ! 最適化: 各kについてV(k;k')を最大にするk'を探す
        do i = 1,nk
            Tvfcn(i) = maxval(vkp(:,i))
            Tpfcn(i) = kgrid(maxloc(vkp(:,i),dim=1))
        end do

        ! 繰り返し計算誤差を確認
        dif1 = maxval(abs((Tvfcn-vfcn)/vfcn))
        dif2 = maxval(abs((Tpfcn-pfcn)/pfcn))

        ! 価値関数・政策関数をアップデート
        vfcn = Tvfcn
        pfcn = Tpfcn
        write (*,"('iteration index: ', i3, ', iteration diff of value: ', f9.5, ', iteration diff of policy: ', f9.5)") it, dif1, dif2

    end do

    ! 計算時間をカウント終了
    call cpu_time( time_end )
    write (*,*) ""
    write (*,"(' Program finished sucessfully in', f12.9, ' seconds')") time_end - time_begin

    ! 計算結果をコマンドプロンプトに表示
    write (*,*) ""
    write (*,*) "-+- Parameter values -+-"
    write (*,*) ""
    write (*,"('beta=', f5.2, ', gamma=', f5.2, ', alpha=', f5.2, ',delta=', f5.2)") beta, gamma, alpha, delta
    write (*,*) ""
    write (*,"('kmin=', f5.2, ', kmax=', f5.2, ', #grid', i5)") kmin, kmax, nk
    write (*,*) ""

    ! 解析的解
    AA = (1.0-beta)**(-1) * (log(1.0-alpha*beta) + ((alpha*beta)/(1.0-alpha*beta))*log(alpha*beta))
    BB = alpha/(1.0-alpha*beta)
    v_true = AA + BB*log(kgrid)
    p_true = beta*alpha*(kgrid**alpha)

    ! 外部ファイルに出力
    open (10, file='result.csv', status='replace')
    do i = 1,nk
        write (10,10) kgrid(i), ',', vfcn(i), ',', pfcn(i), ',', v_true(i), ',', p_true(i)
        10 format(f0.6,a1,f0.6,a1,f0.6,a1,f0.6,a1,f0.6)
    end do
    close (10)

end program