program question_3
    implicit none
    real*8 :: unirand(10000), gaussrand(10000), exprand(10000),R, unirand2(10000), x
    integer :: i
    real*8::u,v,s,array(20000),t,sigma

    open(unit = 11, file = "gaussrandom.dat")
    call random_number(unirand)
    call random_number(unirand2)

    do i = 1, size(unirand)
        exprand(i) = (log(1 - (unirand(i))))/(-2)
    end do

    sigma=2
    do i=1,10000
    s=2.0
        do while(s>1)
            call random_number(u)
            call random_number(v)
            u=2*u-1
            v=2*v-1
            s=u*u+v*v
        enddo
        t=sigma*sqrt(-(2*log(s))/s)
        array(i)=t*u
        array(10000 + i) = t * v

    enddo
    open(unit=109,file="gaus.dat")
        do i=1,20000
            write(109,*)array(i)
        enddo
    close(109)

    ! do i = 1, size(unirand)/2
    !     R = sqrt((unirand(i)**2) + (unirand2(i)**2))
    !     if ( R<=1 ) then
    !         x = sqrt(-8 * log(R**2)) * (unirand(i)/R)
    !         write(11 , '(f15.6)') x
    !     end if
    ! end do
    ! do i = size(unirand)/2, size(unirand)
    !     R = sqrt((unirand(i)**2) + (unirand2(i)**2))
    !     if ( R<=1 ) then
    !         x = -sqrt(-8 * log(R**2)) * (unirand(i)/R)
    !         write(11 , '(f15.6)') x
    !     end if
    ! end do

    open(unit = 10, file = "exprandom.dat")
    do i = 1, size(exprand)
        write(10 , '(f15.6)') exprand(i)
    end do
    close(10)

end program question_3