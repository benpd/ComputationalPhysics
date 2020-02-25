program contint
    use numtype
    implicit none

    integer :: nint, ifail, i
    real(dp) :: weight(100), abscis(100), r0, t
    complex(dp) :: z, z0, zz, ress
    complex(dp), external :: ff

    nint = 30

    call d01bcf (0, 0._dp, 2*pi, 0._dp, 0._dp, nint, weight, abscis, ifail)

    z0 = 0        !For question, (0, pi/2, pi)
    r0 = 0.01
    ress = 0._dp

    do i = 1, nint
        t = abscis(i)
        z = r0 * exp(iic * t)
        zz = z0 + z
        ress = ress + (ff(zz)*0.5_dp * z * weight(i))
    end do
    
    print '(2f7.3)', ress

end program contint

!_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-

function ff(z)
    use numtype
    implicit none

    complex(dp) :: z, ff

    ff = (1/tan(pi*z)) 

end function ff