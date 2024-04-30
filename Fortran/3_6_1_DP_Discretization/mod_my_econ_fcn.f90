module mod_my_econ_fcn

    implicit none
    contains

    function CRRA( cons, gamma )
    ! CRRA�^���p�֐�
    ! cons: �����
    ! gamma: ���ΓI�댯���x
        use mod_types
        !****** input ******
        real(dp), intent(in) :: cons, gamma
        real(dp) :: CRRA
        !-----------------------------
        if (gamma /= 1.0) then
            CRRA = cons**(1-gamma)/(1-gamma)
        else
            CRRA = log(cons)
        end if
    end function CRRA

end module mod_my_econ_fcn
