
module testing

    implicit none

contains

    logical function is_close(value, expected, reltol, name)
    double precision, intent(in) :: value, expected, reltol
    character(len=*), intent(in) :: name

    double precision relerr

    relerr = dabs((value - expected)/expected)
    is_close = relerr .le. reltol
    if (.not. is_close) then
        print *, 'is_close failed for: ', name
        print *, 'expected value:  ', expected
        print *, 'actual value:    ', value
        print *, 'rel. err:        ', relerr
        print *, 'allowed rel.err: ', reltol
    end if
    end function is_close

end module testing
