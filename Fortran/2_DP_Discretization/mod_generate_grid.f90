MODULE mod_generate_grid
!
!  Purpose:
!    Generate meshed grids (evaluation points).
!    The grids are spaced uniformly or exponentially.
!    (1) grid_uniform
!    (2) grid_exp
!    (3) grid_double_exp
!    (4) grid_triple_exp
!    (5) stretch_space
!
!  Record of revisions:
!     Date     Programmer  Description of change
!  ==========  ==========  ======================
!  10/21/2010  T. Yamada   Collect subroutines
!  10/22/2014  T. Yamada   Add stretch_space

    IMPLICIT NONE

    !****** set double precision ******
    integer, parameter :: prec = selected_real_kind(p = 15, r = 307)
    !----------------------------------

    CONTAINS

    SUBROUTINE grid_uniform( mink, maxk, num_grid, grid )
    !
    !  Purpose:
    !    Generate a uniform grid between [maxk, mink].
    !    The number of grids is `num_grid'.
    !
    !  Record of revisions:
    !     Date     Programmer  Description of change
    !  ==========  ==========  =====================
    !  07/19/2003  T. Yamada   Original code
        !****** input ******
        real(prec), intent(in) :: maxk
        real(prec), intent(in) :: mink
        integer, intent(in) :: num_grid
        !****** output ******
        real(prec), intent(out), dimension(num_grid) :: grid
        !****** local variables ******
        integer :: i
        real(prec) :: increment
        !-----------------------------
        increment = (maxk - mink) / dble(num_grid-1)
        do i = 1,num_grid
            grid(i) = (i-1)*increment+mink
        end do
        ! avoid rounding error
        if ( grid(num_grid) /= maxk ) then
            grid(num_grid) = maxk
        end if
    END SUBROUTINE grid_uniform


    SUBROUTINE grid_exp( mink, maxk, num_grid, grid )
    !
    !  Purpose:
    !    Generate exponentially-spaced grids.
    !    "mink" must be 0 or positive.
    !
    !  Record of revisions:
    !     Date     Programmer  Description of change
    !  ==========  ==========  =====================
    !  10/21/2010  T. Yamada   Original code
    !  10/13/2014  T. Yamada   Negative mink allowed
        !****** input ******
        real(prec), intent(in) :: maxk
        real(prec), intent(in) :: mink
        integer, intent(in) :: num_grid
        !****** output ******
        real(prec), intent(out), dimension(num_grid) :: grid
        !****** local variables ******
        real(prec) :: maxd
        real(prec), dimension(num_grid) :: mesh
        !-----------------------------
        maxd = log(maxk+1.0)
        CALL grid_uniform( mink, maxd, num_grid, mesh )
        grid = exp(mesh) - 1.0
    END SUBROUTINE grid_exp


    SUBROUTINE grid_double_exp( mink, maxk, num_grid, grid )
    !
    !  Purpose:
    !    Generate double exponentially-spaced grids.
    !    "mink" must be 0 or positive.
    !
    !  Record of revisions:
    !     Date     Programmer  Description of change
    !  ==========  ==========  =====================
    !  10/21/2010  T. Yamada   Original code
        !****** input ******
        real(prec), intent(in) :: maxk
        real(prec), intent(in) :: mink
        integer, intent(in) :: num_grid
        !****** output ******
        real(prec), intent(out), dimension(num_grid) :: grid
        !****** local variables ******
        real(prec) :: maxd
        real(prec), dimension(num_grid) :: mesh
        !-----------------------------
        maxd = log(log(maxk+1.0)+1.0)
        CALL grid_uniform( mink, maxd, num_grid, mesh )
        grid = exp(exp(mesh)-1.0)-1.0
    END SUBROUTINE grid_double_exp


    SUBROUTINE grid_triple_exp( mink, maxk, num_grid, grid )
    !
    !  Purpose:
    !    Generate triple exponentially-spaced grids.
    !    "mink" must be 0 or positive.
    !
    !  Record of revisions:
    !     Date     Programmer  Description of change
    !  ==========  ==========  =====================
    !  10/21/2010  T. Yamada   Original code
        !****** input ******
        real(prec), intent(in) :: maxk
        real(prec), intent(in) :: mink
        integer, intent(in) :: num_grid
        !****** output ******
        real(prec), intent(out), dimension(num_grid) :: grid
        !****** local variables ******
        real(prec) :: maxd
        real(prec), dimension(num_grid) :: mesh
        !-----------------------------
        maxd = log(log(log(maxk+1.0)+1.0)+1.0)
        CALL grid_uniform( mink, maxd, num_grid, mesh )
        grid = exp(exp(exp(mesh)-1.0)-1.0)-1.0
    END SUBROUTINE grid_triple_exp


    SUBROUTINE stretch_space( grid, num_grid, mink, maxk )
    !
    !  Purpose:
    !    Linear transformation in case of negative mink.
    !
    !  Record of revisions:
    !     Date     Programmer  Description of change
    !  ==========  ==========  =====================
    !  10/22/2014  T. Yamada   Original code
        !****** input ******
        real(prec), intent(in) :: maxk
        real(prec), intent(in) :: mink
        integer, intent(in) :: num_grid
        !****** output ******
        real(prec), intent(inout), dimension(num_grid) :: grid
        !****** local variables ******
        integer :: i
        real(prec), dimension(num_grid) :: temp
        real(prec), dimension(num_grid-1) :: gratio
        !-----------------------------
        do i = 1,num_grid-1
            gratio(i) = (grid(i+1)-grid(i))/(grid(num_grid)-grid(1))
        end do
        temp(1) = mink
        do i = 1,num_grid-1
            temp(i+1) = gratio(i)*(maxk-mink) + temp(i)
        end do
        grid = temp
    END SUBROUTINE stretch_space

END MODULE mod_generate_grid
