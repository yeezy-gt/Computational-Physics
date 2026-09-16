program question2
    implicit none
    real*8 :: rand(1000000), corr, mean, std, meanosq, cross
    integer :: i, f

    call random_number(rand)
    open(unit = 10, file = "./random.dat")
    do i = 1, size(rand)
        write(10, '(f10.9)') rand(i)        
    end do
    close(10)

    mean = sum(rand)/size(rand)
    meanosq = 0
    do i = 1, size(rand)
        meanosq = meanosq + (rand(i))**2
    end do
    meanosq = meanosq/size(rand)
    
    open(unit = 100, file = "./corr.dat")
    do f = 1, 100
        cross = 0
        do i = 1, size(rand)-f
            cross = cross + rand(i)*rand(i+f)
        end do
        cross = cross/(size(rand)-f)
        corr = (cross - (mean**2))/(meanosq - (mean**2))
        write (100, '(f15.9)') corr
    end do
    close(100)

    std = sqrt(meanosq - (mean**2))

    print *, "the mean is "
    print *, mean    
    print *, "the correlation function is "
    print *, corr
    print *, "the standard deviation is "
    print *, std

end program question2