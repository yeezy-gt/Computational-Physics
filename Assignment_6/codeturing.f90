program turing
    implicit none
    integer, parameter ::lx = 60, ly = 60
    real*8 :: old_tempA(1:lx, 1:ly), tempA(1:lx, 1:ly)
    real*8 :: old_tempB(1:lx, 1:ly), tempB(1:lx, 1:ly)
    integer :: i, j, kk, jj, ii
    real*8 :: bound_temp, increment_temp, prefactor, densityA, densityB
    integer :: xp, xn, yp, yn

    character (len = 30) :: filename, filename1
    integer :: iunit, iunitB, ci, cii, test, counter, niter

    !Parameters for simulation
    real*8 :: dx = 1.0d0, dy = 1.0d0, dt = 0.002d0
    real*8 :: toler= 0.0001d0
    real*8 :: Diff_a = 1.0d0, Diff_b = 100.0d0
    real*8 :: alpha = 0.05d0, beta = 10.0d0
    integer :: n_snapshots = 100, tot_niter = 40000

    call random_number(old_tempA)
    call random_number(old_tempB)
    ! old_tempA = 0.0d0
    ! old_tempB = 0.0d0
    ! old_tempA = 0.0d0
    ! old_tempA(lx/4, 3*ly/4) = 1.0d0
    ! old_tempA(3*lx/4, 3*ly/4) = 1.0d0
    ! old_tempA(lx/2, ly/2) = 0.1d0
    ! old_tempA(3*lx/4, ly/4) = 7.0d0

    old_tempA = old_tempA - 0.5d0
    old_tempB = old_tempB - 0.5d0
    tempA = old_tempA
    tempB = old_tempB

    iunit = 71
    iunitB = 72
    ci = 0

    write(filename, '("initialise_", i0, ".dat")') ci
    write(filename1, '("initialiseB_", i0, ".dat")') ci

    open(unit = iunit, file = filename)
    do ii = 1, lx
        do jj = 1, ly
            write(iunit, *) ii, jj, old_tempA(ii, jj)
        end do
    end do
    close(iunit)

    open(unit = iunitB, file = filename1)
    do ii = 1, lx
        do jj = 1, ly
            write(iunitB, *) ii, jj, old_tempB(ii, jj)
        end do
    end do
    close(iunitB)

    test = 0
    counter = 0
    prefactor = 0.5d0*dx*dx*dy*dy/((dx*dx) + (dy*dy))

    do niter = 1, tot_niter
        counter = counter + 1
        test = 0

        do jj = 1, ly
            xp = jj+1
            xn = jj-1

            do ii = 1, lx
                if (jj == 1) xn = ly
                if (jj==ly) xp = 1

                yp = ii + 1
                yn = ii - 1

                if (ii == 1) yn = lx
                if (ii == lx) yp = 1

                tempA(ii,jj) = old_tempA(ii, jj) + dt*Diff_a*(old_tempA(yn, jj) + old_tempA(yp, jj) + old_tempA(ii, xp) + &
                & old_tempA(ii, xn) - 4.0d0*old_tempA(ii, jj))
                tempA(ii,jj) = tempA(ii,jj) + dt*(old_tempA(ii, jj) - (old_tempA(ii, jj))**3 + alpha - old_tempB(ii,jj))
                
                tempB(ii,jj) = old_tempB(ii, jj) + dt*Diff_b*(old_tempB(yn, jj) + old_tempB(yp, jj) + old_tempB(ii, xp) + &
                & old_tempB(ii, xn) - 4.0d0*old_tempB(ii, jj))
                tempB(ii,jj) = tempB(ii,jj) + dt*beta*(old_tempA(ii, jj) - old_tempB(ii,jj))
            end do
        end do
        old_tempA = tempA
        old_tempB = tempB

        if ( mod(niter, n_snapshots) == 0 ) then
            write(filename, '("initialise_", i0, ".dat")') niter
            write(filename1, '("initialiseB_", i0, ".dat")') niter

            densityA = 0.0d0
            densityB = 0.0d0

            open(unit = iunit, file = filename)
            open(unit = iunitB, file = filename1)
            do ii = 1, lx
                do jj = 1, ly
                    write(iunit, *) ii, jj, tempA(ii,jj)
                    write(iunitB, *) ii, jj, tempB(ii,jj)
                    densityA = densityA + tempA(ii,jj)
                    densityB = densityB + tempB(ii,jj)
                end do
            end do
            close(iunit)
            close(iunitB)
            write(*, *) niter, densityA, densityB
        end if
    end do

end program turing