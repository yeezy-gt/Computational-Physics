program assignment3
    implicit none
    integer :: i, j, k, L, p, a, b, c, d, e, f,  niter, time, mm, nn, kk, N, T_temp
    real :: r, M, mag, Ei, Ef, dE, u, h
    real*8 :: Energy, avemag, chi, aveen, cv, aveM2, aveE2, binders, aveM4
    real :: T, J_ising = 1.0
    integer, dimension(:,:,:), allocatable :: spin
    real, dimension(:), allocatable :: Magnetisation, Energyarray
    integer :: seed 

    seed = 69420
    print *, 'Enter the value of T'
    read *, T
    print *, 'Enter the number of lattice points in one dimension:'
    read *, L
    print *, 'Enter the number of iterations'
    read *, niter

    allocate(spin(L,L,L))
    Energy = 0.0     !Instantaneous energy of the lattice
    M = 0.0     !Instantaneous magnetisation of the lattice
    N = L*L*L   !Total number of spins in lattice

    call random_seed

    !initialising the lattice
    open(unit = 10, file = "initial_ising.dat")
    p = 0
    do i = 1,L
        do j = 1, L
            do k = 1, L
                call random_number(r)
                spin(k, j, i) = -1
                if (r<0.5) then
                    spin(k, j, i) = 1
                else
                    spin(k, j, i) = -1
                end if
                write(10, *) float(i), float(j), float(k), float(p), float(spin(k,j,i))
            end do
        end do
    end do
    close(10)

    !calculate initial magnetisation and energy
    do i = 1, L
        do j = 1, L
            do k = 1, L
                !Identifying the 6 nearest neighbours
                a = i+1
                b = i-1
                c = j+1
                d = j-1
                e = k+1
                f = k-1

                !Applying Periodic Boundary conditions
                if (i == L) a = 1
                if (i == 1) b = L
                if (j == L) c = 1
                if (j == 1) d = L
                if (k == L) e = 1
                if (k == 1) f = L
                M = M + spin(i,j,k)
                Energy = Energy - J_ising*float(spin(i,j,k)*(spin(a,j,k)+spin(b,j,k)+spin(i,c,k)+spin(i,d,k)+spin(i,j,e) &
                + spin(i,j,f)))
            end do
        end do
    end do

    mag = M/float(N)
    Energy = Energy * 0.5d0
    print *, 'initial energy E, E per spin = ', Energy, Energy/float(N)
    print *, 'initial magnetisation M, M per spin = ', M, mag

    ! INITIALIZATION COMPLETE

    open(unit=11, file = "ising_T2_N40_init_random.dat")
    do time = 1, niter
        do mm = 1, L
            do nn = 1, L
                do kk = 1, L
                    !chosing a lattice site
                    call random_number(r)
                    i = int(r*float(L))+1
                    call random_number(r)
                    j = int(r*float(L))+1
                    call random_number(r)
                    k = int(r*float(L))+1
                    
                    !Identifying the 6 nearest neighbours
                    a = i+1
                    b = i-1
                    c = j+1
                    d = j-1
                    e = k+1
                    f = k-1

                    !Applying Periodic Boundary conditions
                    if (i == L) a = 1
                    if (i == 1) b = L
                    if (j == L) c = 1
                    if (j == 1) d = L
                    if (k == L) e = 1
                    if (k == 1) f = L

                    Ei = - J_ising*float(spin(i,j,k)*(spin(a,j,k)+spin(b,j,k)+spin(i,c,k)+spin(i,d,k)+spin(i,j,e)+ spin(i,j,f)))
                    spin(i, j, k) = -spin(i, j, k)  !trial flip
                    Ef = - J_ising*float(spin(i,j,k)*(spin(a,j,k)+spin(b,j,k)+spin(i,c,k)+spin(i,d,k)+spin(i,j,e)+ spin(i,j,f)))
                    dE = Ef - Ei

                    if (dE<=0.0) then
                        Energy = Energy + dE
                        M = M + (2*float(spin(i,j, k)))
                    else 
                        u = exp(-dE/(T))
                        call random_number(h)
                        if (h<u) then
                            Energy = Energy + dE
                            M = M + (2*float(spin(i,j, k)))
                        else
                            spin(i,j,k) = -spin(i,j,k)
                        end if
                    end if
                end do
            end do
        end do
        write(11,*) real(time), M/float(N), Energy/float(N)
    end do
    close(11) 3.8

!     open(unit = 100, file = 'q7l10.dat')
!     L = 10
!     niter = 1000000
!     allocate(spin(L,L,L))
!     allocate(Magnetisation(niter))
!     allocate(Energyarray(niter))
!     do T_temp = 470, 380, -2

!         T = float(T_temp)/100.0
!         print *, T
!         Energy = 0.0     !Instantaneous energy of the lattice
!         M = 0.0     !Instantaneous magnetisation of the lattice
!         N = L*L*L   !Total number of spins in lattice
!         open(unit = 10, file = "initial_ising.dat")
!         p = 0
!         do i = 1,L
!             do j = 1, L
!                 do k = 1, L
!                     call random_number(r)
!                     spin(k, j, i) = -1
!                     if (r<0.5) then
!                         spin(k, j, i) = 1
!                     else
!                         spin(k, j, i) = -1
!                     end if
!                 end do
!             end do
!         end do
!         close(10)

!         !calculate initial magnetisation and energy
!         do i = 1, L
!             do j = 1, L
!                 do k = 1, L
!                     !Identifying the 6 nearest neighbours
!                     a = i+1
!                     b = i-1
!                     c = j+1
!                     d = j-1
!                     e = k+1
!                     f = k-1

!                     !Applying Periodic Boundary conditions
!                     if (i == L) a = 1
!                     if (i == 1) b = L
!                     if (j == L) c = 1
!                     if (j == 1) d = L
!                     if (k == L) e = 1
!                     if (k == 1) f = L
!                     M = M + spin(i,j,k)
!                     Energy = Energy - J_ising*float(spin(i,j,k)*(spin(a,j,k)+spin(b,j,k)+spin(i,c,k)+spin(i,d,k)+spin(i,j,e) &
!                     + spin(i,j,f)))
!                 end do
!             end do
!         end do

!         mag = M/float(N)
!         Energy = Energy * 0.5d0

!         do time = 1, niter
!             do mm = 1, L
!                 do nn = 1, L
!                     do kk = 1, L
!                         !chosing a lattice site
!                         call random_number(r)
!                         i = int(r*float(L))+1
!                         call random_number(r)
!                         j = int(r*float(L))+1
!                         call random_number(r)
!                         k = int(r*float(L))+1
                        
!                         !Identifying the 6 nearest neighbours
!                         a = i+1
!                         b = i-1
!                         c = j+1
!                         d = j-1
!                         e = k+1
!                         f = k-1
    
!                         !Applying Periodic Boundary conditions
!                         if (i == L) a = 1
!                         if (i == 1) b = L
!                         if (j == L) c = 1
!                         if (j == 1) d = L
!                         if (k == L) e = 1
!                         if (k == 1) f = L
    
!                         Ei = - J_ising*float(spin(i,j,k)*(spin(a,j,k)+spin(b,j,k)+spin(i,c,k)+spin(i,d,k)+spin(i,j,e)+ spin(i,j,f)))
!                         spin(i, j, k) = -spin(i, j, k)  !trial flip
!                         Ef = - J_ising*float(spin(i,j,k)*(spin(a,j,k)+spin(b,j,k)+spin(i,c,k)+spin(i,d,k)+spin(i,j,e)+ spin(i,j,f)))
!                         dE = Ef - Ei
    
!                         if (dE<=0.0) then
!                             Energy = Energy + dE
!                             M = M + (2*float(spin(i,j, k)))
!                         else 
!                             u = exp(-dE/(T))
!                             call random_number(h)
!                             if (h<u) then
!                                 Energy = Energy + dE
!                                 M = M + (2*float(spin(i,j, k)))
!                             else
!                                 spin(i,j,k) = -spin(i,j,k)
!                             end if
!                         end if
!                     end do
!                 end do
!             end do
!             Magnetisation(time) = M/float(N)
!             Energyarray(time) = Energy/float(N)
!         end do
        
!         aveE2 = 0.0
!         aveM2 = 0.0
!         aveen = 0.0
!         avemag = 0.0

!         do i = 10000, niter
!             aveE2 = aveE2 + (Energyarray(i))**2
!         end do
!         aveE2 = aveE2/(niter-10000)

!         do i = 10000, niter
!             aveM2 = aveM2 + (Magnetisation(i))**2
!         end do
!         aveM2 = aveM2/(niter-10000)

!         do i = 10000, niter
!             aveM4 = aveM4 + (Magnetisation(i))**4
!         end do
!         aveM4 = aveM4/(niter-10000)

!         binders = 1 - aveM4/(3*((aveM2)**2))

!         aveen = sum(Energyarray(10000:))/(niter-10000)
!         avemag = sum(abs(Magnetisation(10000:)))/(niter-10000)
!         chi = (aveM2 - ((avemag)**2))/((niter-10000))
!         cv = (aveE2 - ((aveen)**2))/((niter-10000))
!         write(100, *) T, avemag, chi, aveen, cv, binders
!     end do

end program assignment3
