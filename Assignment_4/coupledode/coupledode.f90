program coupledode
    implicit none
    integer, parameter :: nop = 50
    integer :: niter, i, dt1000, kk, kkp1, kkn1
    real*8 :: init_t, init_y(nop), init_vel(nop)
    real*8 :: yim(nop), vim(nop), y_temp, v_temp
    real*8, parameter :: km = 1.0d0, mass = 1.0d0
    real*8 :: dt, t, dtby2, ke, pe, zpos, twopi, r, theta(nop)
    real*8 :: f0(nop), f1(nop), f2(nop), f3(nop), f0v(nop), f1v(nop), f2v(nop), f3v(nop)
    real*8 :: y1(nop), y2(nop), y3(nop), v1(nop), v2(nop), v3(nop), force(nop)
    character(len = 30) :: charac_a, charac_e, charac_m, charac_i

    charac_m = 'x_energy_vl1.0_dt_.dat'
    charac_i = 'kx_coor_dt_'

    niter = 2000
    init_t = 0.0d0
    init_y = 0.0d0
    init_vel = 0.0d0
    dt = 0.02d0
    dt1000 = int(dt*1000.d0)
    dtby2 = dt/2.0d0

    t = init_t
    yim = init_y
    yim(1) = 0.8d0
    yim(26) = 0.8d0
    vim = init_vel

    zpos = 0.0d0
    r = 5.0d0
    twopi = 4.0*asin(1.0d0)
    do i = 1, nop
        theta(i) = dfloat(i)*twopi/(50.0d0)
    end do

    charac_a = charac_m
    ! call addnumtostring(charac_a, dt1000)
    open(unit = 80, file = charac_a, form = 'formatted')
    open(unit = 81, file = 'posi.xyz', form = 'formatted')
    do i = 1, niter

        t = t + dt
        do kk = 1, nop
            kkp1 = kk + 1
            kkn1 = kk - 1
            if (kk == 1) kkn1 = nop
            if (kk == nop) kkp1 = 1
            
            f0(kk) = vim(kk)
            f0v(kk) = km*(yim(kkp1) + yim(kkn1) - 2.0d0*yim(kk))
        end do

        do kk = 1, nop
            y1(kk) = yim(kk) + f0(kk)*dtby2
            v1(kk) = vim(kk) + f0v(kk)*dtby2
        end do

        do kk = 1, nop
            kkp1 = kk + 1
            kkn1 = kk - 1
            if (kk == 1) kkn1 = nop
            if (kk == nop) kkp1 = 1
            
            f1(kk) = v1(kk)
            f1v(kk) = km*(y1(kkp1) + y1(kkn1) - 2.0d0*y1(kk))
        end do

        do kk = 1, nop
            y2(kk) = yim(kk) + dtby2*f1(kk)
            v2(kk) = vim(kk) + dtby2*f1v(kk)
        end do

        do kk = 1, nop
            kkp1 = kk + 1
            kkn1 = kk - 1
            if (kk == 1) kkn1 = nop
            if (kk == nop) kkp1 = 1

            f2(kk) = v2(kk)
            f2v(kk) = km*(y2(kkp1) + y2(kkn1) - 2.0d0*y2(kk))
        end do

        do kk = 1, nop
            y3(kk) = yim(kk) + dt*f2(kk)
            v3(kk) = vim(kk) + dt*f2v(kk)
        end do

        do kk = 1, nop
            kkp1 = kk + 1
            kkn1 = kk - 1
            if (kk == 1) kkn1 = nop
            if (kk == nop) kkp1 = 1
            
            f3(kk) = v3(kk)
            f3v(kk) = km*(y3(kkp1) + y3(kkn1) - 2.0d0*y3(kk))
        end do

        yim = yim + dt*(f0 + 2.0d0*f1 + 2.0d0*f2 + f3)/6.0d0
        force = (f0v+ 2.0d0*f1v + 2.0d0*f2v + f3v)/6.0d0
        vim = vim + dt*force
        write(*,*) sum(force), sum(vim)
        ke = 0.0d0
        pe = 0.0d0

        do kk = 1, nop
            kkp1 = kk + 1
            kkn1 = kk - 1
            if (kk == 1) kkn1 = nop
            if (kk == nop) kkp1 = 1

            ke = ke + mass*0.5*vim(kk)*vim(kk)
            pe = pe + km*0.25d0*((yim(kkp1)-yim(kk))**2 + (yim(kkn1)-yim(kk))**2)
        end do
        write(80, fmt = '(9g25.15)') t, ke + pe, ke, pe
        

        if(mod(i, 10) == 0) then
            write(81, *) 50
            write(81, *)
            ! call addnumtostring(charac_a, i)

            do kk = 1, nop
                write(81, '(A,3g25.15)')'N', r*cos(theta(kk)), r*sin(theta(kk)),yim(kk)
            end do

        end if

    end do
    close(81)
    close (80)
end program coupledode

! subroutine addnumtostring(charac_a, i)
!     implicit none
!     character(len = 30) :: charac_a, filename
!     integer :: i, iunit, io

!     write(filename, '(charac_a, io, ".xyz")') i
!     open(newunit = iunit, file = filename)
! end subroutine addnumtostring