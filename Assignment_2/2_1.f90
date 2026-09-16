program question1
    implicit none
    real*8:: x, trap_integral, trap_integral_sin, trap_integral_gauss, gauss

    ! open(unit = 10, file = "trapdata.dat")
    ! Twrite(10,*) trap_integral(100)
    ! write(10,*) trap_integral(1000)
    ! write(10,*) trap_integral(10000)
    ! write(10,*) trap_integral(100000)
    ! close(10)
    
    x = trap_integral_gauss(100000)
    print *, x

end program question1

function trap_integral(n) result(integral)
    implicit none
    real*8 :: integrand, lowlim, uplim, integral, delx
    integer :: i
    integer, intent(in) :: n
    
    ! lowlim : lower limit of the integral
    ! uplim : upper limit of the integral
    ! flow : integrand value at lower limit
    ! fhi : integrand value at upper limit
    ! n : number of divisions to calculate the trapezoidal integral
    ! integral : value of the integral
    ! delx : width of each division
    
    lowlim = 0.0
    uplim = 1.0

    delx = (uplim - lowlim)/real(n)
    integral = delx * (integrand(lowlim) +integrand(uplim))/2
    do i = 1, (n-1)
        integral = integral + (delx * integrand((lowlim + i*delx)))
    end do

    print *, integral

end function trap_integral

function trap_integral_sin(n) result(integral)
    implicit none
    real*8 :: lowlim, uplim, integral, delx
    integer :: i
    integer, intent(in) :: n
    
    ! lowlim : lower limit of the integral
    ! uplim : upper limit of the integral
    ! flow : integrand value at lower limit
    ! fhi : integrand value at upper limit
    ! n : number of divisions to calculate the trapezoidal integral
    ! integral : value of the integral
    ! delx : width of each division

    lowlim = 0.0d0
    uplim = 2*asin(1.0d0)
    
    delx = (uplim - lowlim)/real(n)
    integral = delx * (sin(lowlim) +sin(uplim))/2
    
    do i = 1, (n-1)
        integral = integral + (delx * sin((lowlim + i*delx)))
    end do

    print *, integral

end function trap_integral_sin

function integrand(x) result(ans)
    implicit none
    real*8, intent(in) :: x
    real*8 :: ans

    ! function to calculate the value of the integral
    ans = 4/(1+ (x**2))
end function integrand

function gaussian(x) result(ans)
    implicit none
    real*8, intent(in) :: x
    real*8 :: ans

    ! function to calculate the value of the gaussian function
    ans = (1/sqrt(4*asin(1.0d0)))*(exp(-(x**2)/2))
end function gaussian

function trap_integral_gauss(n) result(integral)
    implicit none
    real*8 :: gaussian, lowlim, uplim, integral, delx
    integer :: i
    integer, intent(in) :: n
    
    ! lowlim : lower limit of the integral
    ! uplim : upper limit of the integral
    ! flow : integrand value at lower limit
    ! fhi : integrand value at upper limit
    ! n : number of divisions to calculate the trapezoidal integral
    ! integral : value of the integral
    ! delx : width of each division
    
    lowlim = -3.0
    uplim = 3.0

    delx = (uplim - lowlim)/real(n)
    integral = delx * (gaussian(lowlim) +gaussian(uplim))/2
    do i = 1, (n-1)
        integral = integral + (delx * gaussian((lowlim + i*delx)))
    end do

    print *, integral

end function trap_integral_gauss