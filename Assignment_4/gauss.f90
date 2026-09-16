program gaussseidel
    implicit none
    real*8, parameter :: f_x = 0.0d0, end_x = 1.0d0
    real*8 :: f_y, end_y
    real*8, parameter :: dx = 0.01d0
    integer, parameter :: nop = int((end_x-f_x)/dx) 
    integer :: i, ll, cond
    real*8 :: x(nop), y(nop), y_old(nop), limit
    real*8, parameter :: dum1 = 1.0d0/(2.0d0 - 10.0d0*dx*dx), dum2 = (1.0d0 - 2.5d0*dx), dum3 = (1.0d0 + 2.5*dx)
    real*8, parameter :: dum4 = -10.0d0*dx*dx
    
    limit = 0.0001d0
    ll = 0
    cond = 0
    open(10, file = 'gauss_seidel.dat', status = 'unknown')
    x(1) = f_x
    x(nop) = end_x
    y(1) = 0.0d0
    y(nop) = 2.0d0

    do i = 2, nop - 1
        x(i) = x(i-1) + dx
        y(i) = (end_y - f_y)*x(i)/(end_x - f_x)
    end do
    
    do
        ll = ll + 1
        if(cond==1) exit
        y_old = y
        do i = 2, nop - 1
            y(i) = dum1*(dum2 * y_old(i+1) + dum3 * y(i-1) + dum4*x(i))
        end do
        cond = 1
        do i = 2, nop - 1
            if(abs(y_old(i) - y(i)) .ge. limit) cond = 0
        end do
    end do
    write(*, *) 'number of iterations to achieve convergence ', ll

    do i = 1, nop
        write(10, *) x(i), y(i)
    end do

end program gaussseidel
