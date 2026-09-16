program laplace1
    implicit none
    integer, parameter ::lx = 35, ly = 35
    real*8 :: old_temp(1:lx, 1:ly), temp(1:lx, 1:ly)
    integer :: i, jj, ii
    real*8 :: bound_temp1, bound_temp2, increment_temp
    integer ::  test, counter

    bound_temp1 = 3.7d0
    bound_temp2 = 0.4d0 
    increment_temp = -0.1d0
    old_temp =0.0d0
    do i = 1, ly
        old_temp(1, i) = bound_temp1
        old_temp(lx, i) = bound_temp2
        if (mod(i, 2) ==0) old_temp(1, i) = bound_temp2
    end do

    

    do i = 2, lx-1
        old_temp(i, 1) = old_temp(1,1) + increment_temp * (i-1)
        old_temp(i, ly) = old_temp(1,ly) + increment_temp * (i-1)
    end do

    open(file = "initializelaplace1.dat", unit = 100)
    do i = 1, ly
        write(100,*) old_temp(:, i)
    end do    
    close(100)

    temp = old_temp
    counter = 0
    test = 0
    do
        counter = counter + 1
        test = 0
        do jj = 2, ly-1
            do ii = 2, lx-1
                temp(ii, jj) = 0.25*(old_temp(ii-1, jj) + old_temp(ii+1, jj) + old_temp(ii, jj-1) + old_temp(ii, jj+1))
            end do
        end do
        do jj = 2, ly-1
            do ii = 2, lx-1
                if ((abs(temp(jj,ii)) - old_temp(jj,ii)) > 0.0001d0) then
                    test = 1
                end if
            end do
        end do
        if (test ==0) exit
        old_temp = temp
    end do
    write(*,*) "number of iterations", counter
    open(unit = 101, file = "laplace1.dat")
    do i = 1, ly
        write(101,*) temp(:, i)
    end do 
    close(101)
end program laplace1
