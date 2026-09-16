program assignment

    implicit none
    integer :: num
    real :: float
    real*8 :: double, mean, mean100, mean10k, mean1m, diff, diff100, diff10k, diff1m
    real*8 :: r(10),r2(10),r3(10,10),r100(100),r10k(10000),r1m(1000000), sum10k(10000), normsum(10000)
    real :: a
    integer :: i, j
    real :: walk

    !q0
    num = 90
    float = 90.0
    double = 90.0
    print *, "Yash Gupta"
    print *, num
    print *, float
    print *, double

    !q1a
    call random_seed(put = [69420, 123456, 762463, 456678, 546768, 789546, 789789, 89789])
    call random_number(r)
    print *, r

    !q1b
    open(unit = 10, file = "test_ran.dat")
    do i = 1, 10
        write(10,'(f11.10)') r(i)        
    end do

    !q1c
    write(10,*) "Changing seed and generating 10 new random numbers"

    !q1d
    call random_number(r2)
    do i = 1, 10
        write(10,'(f11.10)') r2(i)        
    end do
    
    
    call random_number(r3)
    open (unit = 11, file = "test_ran_10_seeds.dat")
    do i = 1, 10
        write(11,'(10f11.10)') r3(:,i)
    end do
    close(11)
    
    !q1e
    write(10,*) "NOW calculating average of 10 random numbers"
    mean = sum(r)/(max(1,size(r)))
    write(10,'(f11.10)') mean

    !q1f
    call random_number(r100)
    call random_number(r10k)
    call random_number(r1m)
    mean100 = sum(r100)/(max(1,size(r100)))
    write(10,'(f11.10)') mean100
    mean10k = sum(r10k)/(max(1,size(r10k)))
    write(10,'(f11.10)') mean10k
    mean1m = sum(r1m)/(max(1,size(r1m)))
    write(10,'(f11.10)') mean1m
    close(10)

    !q1g
    diff = abs(0.50d0 - mean)
    print '(f11.10)', diff
    diff100 = abs(0.50d0 - mean100)
    print '(f11.10)', diff100
    diff10k = abs(0.50d0 - mean10k)
    print '(f11.10)', diff10k    
    diff1m = abs(0.50d0 - mean1m)
    print '(f11.10)', diff1m

    !q1h
    do i = 1, 10000
        call random_number(r10k)
        sum10k(i) = sum(r10k)
    end do

    open (unit = 13, file = "sum_10k.dat")
    do i = 1, size(sum10k)
        write(13,*) sum10k(i)      
    end do
    close(13)

    open (unit = 14, file = "norm_sum_10k.dat") 
    do i = 1, size(sum10k)
        normsum(i) = sum10k(i)/sum(sum10k)        
    end do

    
    do i = 1, size(normsum)
        write(14,*) normsum(i)      
    end do
    close(14)

    !!!!!plotting the distributions

    open (unit = 20, file = "dist_sum_random_walk.dat")
    do i = 1, 10000
        walk = 0.0
        do j = 1, 10000
            call random_number(a)
            if (a > 0.5) then
                walk = walk +1
            else
                walk = walk - 1
            end if
        end do
        write(20, *) walk
    end do
    close(20)
    
    open (unit = 21, file = "dist_sum_random_walk_1m.dat")
    do i = 1, 100000
        walk = 0.0
        do j = 1, 10000
            call random_number(a)
            if (a > 0.5) then
                walk = walk +1
            else
                walk = walk - 1
            end if
        end do
        write(21, *) walk
    end do
    close(21)

    open (unit = 22, file = "dist_sum_random_walk_1m.dat")
    do i = 1, 100000
        walk = 0.0
        do j = 1, 10000
            call random_number(a)
            if (a > 0.5) then
                walk = walk +1
            else
                walk = walk - 1
            end if
        end do
        write(22, *) walk
    end do
    close(22)

    open (unit = 23, file = "dist_sum_random_walk_1m1m.dat")
    do i = 1, 100000
        walk = 0.0
        do j = 1, 100000
            call random_number(a)
            if (a > 0.5) then
                walk = walk +1
            else
                walk = walk - 1
            end if
        end do
        write(23, *) walk
    end do
    close(23)

end program assignment

! subroutine create_hist(data, size, binw, start)
!     implicit none
!     integer :: size
!     real, allocatable :: data(:)
!     integer, allocatable :: frequency(:)
!     real :: binw
!     real, optional :: start
!     real :: pointer
!     integer :: i

!     if(.not. present(start))
!         start = minval(data)
!     end if

!     allocate(data(size))
!     allocate(frequency())
!     pointer = start

!     do i = start, end
        
!     end do
!     do while (pointer>maxval)
        
!     end do

! end subroutine create_hist

