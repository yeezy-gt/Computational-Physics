program question_4part2
    implicit none
    integer :: n, i, j, k
    real*8 :: x(6), y, func, int_mc, var, sigma, p(2)
    real*8 :: length, volume, sqrt2, gauss_dev

    length = 5.0d0
    volume = acos(-1.0d0)**3
    sqrt2 = 1.0d0/sqrt(2.0d0)

    open(unit = 10, file = 'importanceint.dat')
    do k = 1, 8
        n = 10**k    
        int_mc = 0.0d0
        var = 0.0d0
        sigma = 0.0d0

        do i = 1, n
            do j = 1,6
                call random_number(p)
                x(j) = gauss_dev(p)*sqrt2
            end do
            int_mc = int_mc + func(x)
            sigma = sigma + func(x)**2
        end do

        int_mc = int_mc/real(n)
        sigma = sigma/real(n)
        var = sigma - int_mc**2

        int_mc = int_mc * volume
        sigma = volume*sqrt(var/real(n))
        write(10, *) n, ' ', int_mc, ' ', sigma
    end do
    close(10)
end program question_4part2

real*8 function func(x) result(ans)
    implicit none
    real*8 :: x(6), xy
    real*8 :: a, b

    a = 0.5d0
    xy = (x(1)-x(4))**2 + (x(2)-x(5))**2 + (x(3)-x(6))**2
    ans = exp(-a*xy)

end function func

real*8 function gauss_dev(x) result(ans)
    implicit none
    real*8 :: fact, sqr, p, x1, x2, x(2)

    sqr = 0.0d0
    do while((sqr >= 1.0d0) .or. (sqr == 0.0d0))
        call random_number(p)
        x1 = 2*p-1.0d0
        call random_number(p)
        x2 = 2*p-1.0d0
        sqr = x1**2 + x2**2
    end do

    ans = x2 * sqrt(-2.0d0*log(sqr)/sqr)
end function gauss_dev

