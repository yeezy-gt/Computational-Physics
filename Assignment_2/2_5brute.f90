program question_4
    implicit none
    integer:: n,i,j, k
    real*8:: x(6), y,func, int_mc,var, sigma, p(6), length, volume

    open(unit = 10, file = "./mcintegral.dat")
    do k = 1, 9
        n = 10 ** k    
        length = 10.0d0
        volume = (2*length)**6

        int_mc = 0.0d0
        var = 0.0d0
        sigma = 0.0d0

        do i = 1, n
            call random_number(p)
            do j = 1,6
                x(j) = (-length) + (2.0d0 * length * p(j))
            end do
            int_mc = int_mc + func(x)
            sigma = sigma + (func(x))**2
        end do
        int_mc = int_mc/real(n)
        sigma = sigma/real(n)
        var = sigma - (int_mc**2)

        int_mc = volume * int_mc
        sigma = volume*sqrt(var/real(n))
        write(10, *) n, ' ', int_mc, ' ', sigma
    end do    
    close(10)
end program question_4

real*8 function func(x) result(ans)
    implicit none
    real*8::x(6), xx, xy, yy, a, b

    a = 1.0d0
    b = 0.5d0
    
    xx = (x(1)**2) + (x(2)**2) + (x(3)**2)
    yy = (x(4)**2) + (x(5)**2) + (x(6)**2)
    xy = (x(1)-x(4))**2 + (x(2)-x(5))**2 + (x(3)-x(6))**2
    ans = exp(-a*xx -a*yy -b*xy )
end function func
