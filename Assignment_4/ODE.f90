program ODE
    implicit none

    ! call euler(0.0d0, 0.0d0, 0.001d0, int(1.55/0.001)+1)
    ! call mod_euler(0.0d0, 0.0d0, 0.001d0, int(1.55/0.001)+1)
    ! call imp_euler(0.0d0, 0.0d0, 0.001d0, int(1.55/0.001)+1)
    ! call rk4(0.0d0, 0.0d0, 0.01d0, int(1.55/0.01))
    ! call diff2sin(0.0d0, 1.999d0, 0.0d0, 0.01d0, 5000)
end program ODE

subroutine euler(x0, y0, dx, niter) 
    implicit none
    real*8, intent(in)::  dx, y0, x0
    integer, intent(in) :: niter
    integer :: i
    real*8 :: dydx, x, y
    
    y = y0
    x = x0
    open(unit = 10, file = './euler.dat')
    do i = 1, niter
        dydx = 1 + y**2
        y = y + (dydx * dx)
        write(10,*) x + i *dx, y
    end do
    close(10)
end subroutine euler

subroutine mod_euler(x0, y0, dx, niter)
    implicit none
    real*8, intent(in)::  dx, y0, x0
    integer, intent(in) :: niter
    integer :: i
    real*8 :: dydx, x, y, ytemp, dydxtemp

    y = y0
    x = x0
    open(unit = 11, file = './mod_euler.dat')
    do i = 1, niter
        dydx = 1 + y**2
        ytemp = y + (dydx * dx/2)
        dydxtemp = 1 + ytemp**2
        y = y + dx * (dydxtemp)
        write(11,*) x + i*dx , y
    end do
    close(11)
    
end subroutine mod_euler

subroutine imp_euler(x0, y0, dx, niter)
    implicit none
    real*8, intent(in)::  dx, y0, x0
    integer, intent(in) :: niter
    integer :: i
    real*8 :: dydx, x, y, ytemp, dydxtemp1, dydxtemp2

    y = y0
    x = x0
    open(unit = 12, file = './imp_euler.dat')
    do i = 1, niter
        dydxtemp1 = 1 + y**2
        ytemp = y + (dydxtemp1 * dx)
        dydxtemp2 = 1 + ytemp**2
        dydx = (dydxtemp1 + dydxtemp2)/2
        y = y + dx * (dydx)
        write(12,*) x + i*dx , y
    end do
    close(12)  
end subroutine imp_euler

subroutine rk4(x0, y0, h, niter)
    implicit none
    real*8, intent(in) :: x0, y0, h
    integer, intent(in) :: niter
    integer :: i
    real*8 :: dydx0, dydx1, dydx2, dydx3, x, y, ytemp1, ytemp2, ytemp3

    x = x0
    y = y0
    open(unit = 20, file = 'rk4.dat')
    write(20, *) x, y
    do i = 1, niter
        dydx0 = 1 + y**2
        ytemp1 = y + (h/2)*dydx0
        dydx1 = 1 + (ytemp1)**2
        ytemp2 = y + (h/2)*dydx1
        dydx2 = 1 + (ytemp2)**2
        ytemp3 = y + h*dydx2
        dydx3 = 1 + ytemp3**2
        x = x + h
        y = y + (h/6)*(dydx0 + (2*dydx1) + (2*dydx2) + dydx3)
        write(20, *) x, y
    end do
    close(20)
end subroutine rk4

subroutine diff2sin(x0, v0, t0, h, niter)
    implicit none
    real*8, intent(in) :: x0, v0, t0, h
    integer, intent(in) :: niter
    real*8 :: x, v, t, xtemp1, xtemp2, xtemp3, dvdt0, dvdt1, dvdt2, dvdt3, dxdt0, dxdt1, dxdt2, dxdt3, vtemp1, vtemp2, vtemp3
    integer :: i

    x = x0
    v = v0
    t = t0
    open(unit = 50, file = "diff2sin.dat")
    write(50, *) x, v, t
    do i = 1, niter
        dvdt0 = -sin(x)
        dxdt0 = v
        xtemp1 = x + (h/2)*dxdt0
        vtemp1 = v + (h/2)*dvdt0

        dvdt1 = -sin(xtemp1)
        dxdt1 = vtemp1
        xtemp2 = x + (h/2)*dxdt1
        vtemp2 = v + (h/2)*dvdt1
        
        dvdt2 = -sin(xtemp2)
        dxdt2 = vtemp2
        xtemp3 = x + (h)*dxdt2
        vtemp3 = v + (h)*dvdt2

        dvdt3 = -sin(xtemp3)
        dxdt3 = vtemp3

        x = x + (h/6)*(dxdt0 + (2*dxdt1) + (2*dxdt2) + dxdt3)
        v = v + (h/6)*(dvdt0 + (2*dvdt1) + (2*dvdt2) + dvdt3)
        t = t + h
        write(50, *) x, v, t
    end do
    close(50)
end subroutine diff2sin