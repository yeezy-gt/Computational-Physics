module sim_para
    implicit none
    integer           :: zzz = 4280145, zzzz=77777
    integer,parameter :: lx=20, ly=20,lz=20,n_part=1200
    integer,parameter :: kb_T=2,niter=100000, mass=1
    real*8            :: pos(3*n_part), vel(3*n_part), force(3*n_part), acc(3*n_part), new_acc(3*n_part)
    real*8, parameter :: rc=2.50d0, rs=4.5d0, sigma=1.0d0 
    real*8, parameter :: sigma6 = sigma**6, eps=4.0d0, sigma12 = sigma**12
    real*8, parameter :: fc=eps*((12.0d0*sigma12/(rc**13)) - (6.0d0*sigma6/(rc**7)))
    real*8, parameter :: ufc = fc*rc + eps*(((sigma/rc)**12) - ((sigma/rc)**6))
    real*8            :: avr_vel_x,avr_vel_y,avr_vel_z
    integer           :: i,j,k,c, p
    integer,parameter :: max_nbr = int(((5*50*3.14/3)*(n_part)/(lx*ly*lz)*(rs)**3))
    integer           :: no_nbr(n_part), nbr_list(n_part, max_nbr)
    real*8            :: x1, y1, z1, x2, y2, z2, x, y , z, dx, dy, dz, r, lj, pot_energy, ke, gy, gz, gx, delta_t = 0.0025
    real*8            :: lj_force
    real*8            :: new_pot_energy, new_force(3*n_part)
    real*8            :: invr, ir2, ir6
    real*8            :: llx, lly, llz, theoryke, scalef
    real*8, parameter :: avg_rho = dfloat(n_part)/(lx*ly*lz), dr = 0.1, dv = 0.05
    integer,parameter :: len_gr = int(dfloat(2*lx)/(2*dr))
    real*8            :: gr(len_gr), vel_arr(n_part), dis(int(lx/(2*dv)))
end module sim_para

program md
    use sim_para
    implicit none

    dis = 0
    llx = dfloat(lx)
    lly = dfloat(ly)
    llz = dfloat(lz)
    gr =  0
    c = 0

    call posinit
    call forceinit
    call velinit

    open(82, file = 'energy1.dat')
    open(69, file = 'momenta1.dat')
    open(70, file = 'gr1.dat')
    open(71, file = 'maxwell1.dat')

    do k = 1, niter
        if (mod(k, 40) == 0) then
            call neigh_list
        end if
        if (mod(k, 100) == 0) then
            call thermostat 
        end if
        if (mod(k, 3000) == 0) then
            write(*,*) k
        end if
        if (mod(k, 100) == 0 .and. k >= 20000) then        
            call pair_corr
            call speed_dis
            c = c+1
        end if
        
        call pos_update
        write(82, '(4f35.15)') dfloat(k), new_pot_energy/(dfloat(n_part)), ke/dfloat(n_part), (ke+new_pot_energy)/(dfloat(n_part))
        write(69, '(3f35.15)') avr_vel_x, avr_vel_y, avr_vel_z
    end do

    do i = 1, len_gr
        gr(i) = gr(i)/(4*3.1415*dr*i*dr*i*dr*avg_rho*c*n_part)
        write(70, *) gr(i)
        write(71, *) dis(i)/(n_part*c)
    end do

end program md

subroutine posinit 
    use sim_para 
    implicit none
    integer :: nx, ny, nz, px, py, pz 
    real*8,parameter :: lattice_spacing=1.2
    
    nx=int(lx/lattice_spacing)
    ny=int(ly/lattice_spacing) 
    nz=int(lz/lattice_spacing)
    px=-2
    py=-1
    pz=0
    do i = 1, nx-1
        do j=1,ny-1
            do k=1, nz-1
                x = dfloat(i)*lattice_spacing
                y = dfloat(j)*lattice_spacing
                z = dfloat(k)*lattice_spacing
                if (x < lx .and. y < ly .and. z < lz) then
                    px = px + 3 
                    py = py + 3
                    pz = pz + 3
                    if (pz <= 3*n_part) then
                        pos(px) = x
                        pos(py) = y
                        pos(pz) = z
                    else
                        return
                    end if 
                else
                    return
                end if
            end do
        end do
    end do
end subroutine posinit

subroutine forceinit
    use sim_para
    implicit none
    
    new_force = 0.0d0
    new_pot_energy = 0.0d0
    do i = 1, n_part - 1
        x1 = pos(3*i - 2)
        y1 = pos(3*i - 1)
        z1 = pos(3*i)

        do j = 1, no_nbr(i)
            p = nbr_list(i, j)
            if (p>i) then
                x2 = pos(3*p - 2)
                y2 = pos(3*p - 1)
                z2 = pos(3*p    )
                
                x = x1 - x2
                y = y1 - y2
                z = z1 - z2

                if (abs(x) >= (dfloat(lx)/2.0d0)) x=(dfloat(lx)-abs(x))*((-1.0d0*x)/abs(x))
                if (abs(y) >= (dfloat(ly)/2.0d0)) y=(dfloat(ly)-abs(y))*((-1.0d0*y)/abs(y))
                if (abs(z) >= (dfloat(lz)/2.0d0)) z=(dfloat(lz)-abs(z))*((-1.0d0*z)/abs(z))

                r = dsqrt(x*x + y*y + z*z)

                if (r <= rc) then
                    lj=eps*((sigma/r)**12 - (sigma/r)**6) - ufc + fc*r
                    new_pot_energy = new_pot_energy+lj
                    lj_force = eps*((12.0d0*((sigma12)/(r**13)))-(6.0d0*((sigma6)/(r**7)))) - fc

                    new_force(3*i - 2) = new_force(3*i - 2) + lj_force*(x/r)
                    new_force(3*i - 1) = new_force(3*i - 1) + lj_force*(y/r)
                    new_force(3*i    ) = new_force(3*i    ) + lj_force*(z/r)
                    new_force(3*p - 2) = new_force(3*p - 2) - lj_force*(x/r)
                    new_force(3*p - 1) = new_force(3*p - 1) - lj_force*(y/r)
                    new_force(3*p    ) = new_force(3*p    ) - lj_force*(z/r)
                end if
            end if
        end do
    end do
    force = new_force
end subroutine forceinit

subroutine velinit
    use sim_para
    implicit none
    real*8 :: vel_const, avg_vx, avg_vy, avg_vz, ran
    integer:: tn_part

    vel_const = dsqrt(12.0d0)*dfloat(kb_T)

    tn_part = 3*n_part

    do i = 1, tn_part
        call random_number(ran)
        vel(i) = vel_const*(ran - 0.5d0)
    end do

    avg_vx = 0.0d0
    avg_vy = 0.0d0
    avg_vz = 0.0d0

    do i = 1, n_part
        avg_vx = vel(3*i - 2) + avg_vx
        avg_vy = vel(3*i - 1) + avg_vy
        avg_vz = vel(3*i) + avg_vz
    end do

    avg_vx = avg_vx/dfloat(n_part)
    avg_vy = avg_vy/dfloat(n_part)
    avg_vz = avg_vz/dfloat(n_part)
    
    do i = 1, n_part
        vel(3*i - 2) = vel(3*i - 2) - avg_vx
        vel(3*i - 1) = vel(3*i - 1) - avg_vy
        vel(3*i) = vel(3*i) - avg_vz
    end do 
    
end subroutine velinit

subroutine pos_update
    use sim_para
    implicit none
    real*8 :: dt2by2

    dt2by2 = 0.50d0*delta_t**2
    
    !Updating positions
    do i = 1, n_part
        pos(3*i - 2) = pos(3*i - 2) + vel(3*i - 2)*delta_t + dt2by2*force(3*i - 2)
        pos(3*i - 1) = pos(3*i - 1) + vel(3*i - 1)*delta_t + dt2by2*force(3*i - 1)
        pos(3*i    ) = pos(3*i    ) + vel(3*i    )*delta_t + dt2by2*force(3*i    )

        pos(3*i - 2) = modulo(pos(3*i - 2), llx)
        pos(3*i - 1) = modulo(pos(3*i - 1), lly)
        pos(3*i    ) = modulo(pos(3*i    ), llz)
    end do

    !Updating Forces
    new_force = 0.0d0
    new_pot_energy = 0.0d0
    
    do i = 1, n_part - 1
        x1 = pos(3*i - 2)
        y1 = pos(3*i - 1)
        z1 = pos(3*i    )

        do j = 1, no_nbr(i)

            p = nbr_list(i, j)
            if(p>i) then
                x2 = pos(3*p - 2)
                y2 = pos(3*p - 1)
                z2 = pos(3*p    )
                
                x = x1 - x2
                y = y1 - y2
                z = z1 - z2

                if (abs(x) >= (dfloat(lx)/2.0d0)) x=(dfloat(lx)-abs(x))*((-1.0d0*x)/abs(x))
                if (abs(y) >= (dfloat(ly)/2.0d0)) y=(dfloat(ly)-abs(y))*((-1.0d0*y)/abs(y))
                if (abs(z) >= (dfloat(lz)/2.0d0)) z=(dfloat(lz)-abs(z))*((-1.0d0*z)/abs(z))

                r = dsqrt(x*x + y*y + z*z)
                
                if (r <= rc) then
                    lj = eps*((sigma/r)**12-(sigma/r)**6)-ufc+fc*r
                    new_pot_energy = new_pot_energy + lj
                    lj_force = eps*((12.0d0*((sigma12)/(r)**13))-(6.0d0*((sigma6)/(r)**7))) - fc

                    new_force(3*i - 2) = new_force(3*i - 2) + lj_force*(x/r)
                    new_force(3*i - 1) = new_force(3*i - 1) + lj_force*(y/r)
                    new_force(3*i    ) = new_force(3*i    ) + lj_force*(z/r)
                    new_force(3*p - 2) = new_force(3*p - 2) - lj_force*(x/r)
                    new_force(3*p - 1) = new_force(3*p - 1) - lj_force*(y/r)
                    new_force(3*p    ) = new_force(3*p    ) - lj_force*(z/r)
                end if
            end if
        end do
    end do

    !Updating velocity
    ke = 0.0d0
    avr_vel_x = 0.0d0
    avr_vel_y = 0.0d0
    avr_vel_z = 0.0d0

    do i = 1, n_part
        vel(3*i-2)= vel(3*i-2) + (delta_t*0.50d0*(force(3*i-2)+ new_force(3*i-2)))
        vel(3*i-1)= vel(3*i-1) + (delta_t*0.50d0*(force(3*i-1)+ new_force(3*i-1)))
        vel(3*i  )= vel(3*i  ) + (delta_t*0.50d0*(force(3*i  )+ new_force(3*i  )))

        avr_vel_x = avr_vel_x + vel(3*i - 2)
        avr_vel_y = avr_vel_y + vel(3*i - 1)
        avr_vel_z = avr_vel_z + vel(3*i    )

        ke = ke + (vel(3*i - 2)*vel(3*i - 2)) + (vel(3*i - 1)*vel(3*i - 1)) + (vel(3*i)*vel(3*i))
    end do
    avr_vel_x = avr_vel_x/n_part
    avr_vel_y = avr_vel_y/n_part
    avr_vel_z = avr_vel_z/n_part
    
    ke = 0.50d0*dfloat(mass)*ke
    force = new_force
end subroutine pos_update

subroutine thermostat
    use sim_para
    implicit none
    
    theoryke = 1.5d0*n_part*kb_T
    scalef = dsqrt(theoryke/ke)
    vel = vel*scalef

end subroutine thermostat

subroutine neigh_list
    use sim_para
    implicit none
    integer :: ij, jk

    no_nbr = 0
    nbr_list = 0
    do ij = 1, n_part - 1
        x1 = pos(3*ij - 2)
        y1 = pos(3*ij - 1)
        z1 = pos(3*ij    )

        do jk = ij+1, n_part
            x2 = pos(3*jk - 2)
            y2 = pos(3*jk - 1)
            z2 = pos(3*jk    )
            
            x = x1 - x2
            y = y1 - y2
            z = z1 - z2

            if (abs(x) >= (dfloat(lx)/2.0d0)) x=(dfloat(lx)-abs(x))*((-1.0d0*x)/abs(x))
            if (abs(y) >= (dfloat(ly)/2.0d0)) y=(dfloat(ly)-abs(y))*((-1.0d0*y)/abs(y))
            if (abs(z) >= (dfloat(lz)/2.0d0)) z=(dfloat(lz)-abs(z))*((-1.0d0*z)/abs(z))

            r = dsqrt(x*x + y*y + z*z)

            if (r < rs) then
                no_nbr(ij) = no_nbr(ij) + 1
                nbr_list(ij, no_nbr(ij)) = jk
            end if 
        end do
    end do
end subroutine neigh_list

subroutine pair_corr
    use sim_para
    implicit none
    integer :: ij, jk, loc

    do ij = 1, n_part
        x1 = pos(3*ij - 2)
        y1 = pos(3*ij - 1)
        z1 = pos(3*ij    )

        do jk = 1, n_part
            if (ij /= jk) then
                x2 = pos(3*jk - 2)
                y2 = pos(3*jk - 1)
                z2 = pos(3*jk    )
                
                x = x1 - x2
                y = y1 - y2
                z = z1 - z2

                if (abs(x) >= (dfloat(lx)/2.0d0)) x=(dfloat(lx)-abs(x))*((-1.0d0*x)/abs(x))
                if (abs(y) >= (dfloat(ly)/2.0d0)) y=(dfloat(ly)-abs(y))*((-1.0d0*y)/abs(y))
                if (abs(z) >= (dfloat(lz)/2.0d0)) z=(dfloat(lz)-abs(z))*((-1.0d0*z)/abs(z))

                r = dsqrt(x*x + y*y + z*z)
                loc = int(r/(dr)) + 1
                gr(loc) = gr(loc) + 1.0d0
            end if
        end do
    end do

end subroutine pair_corr

subroutine speed_dis
    use sim_para
    implicit none
    integer::counter
    real*8::vx,vy,vz, v_r, shell_v
    vel_arr=0
    counter=0
    do i=1,n_part
        counter = counter+1
        vx = vel(3*i-2)
        vy = vel(3*i-1)
        vz = vel(3*i)
        v_r= dsqrt((vx*vx) + (vy*vy) + (vz*vz))
        vel_arr(counter) = v_r
    end do
    
    shell_v=0
    do i=1,size(dis)
        shell_v = shell_v + dv
        do j = 1, size(vel_arr)
            if ((shell_v-dv<vel_arr(j)).and.(vel_arr(j)<shell_v)) then
                dis(i) = dis(i) + 1
            end if
        end do
    end do
endsubroutine speed_dis

! subroutine Diffusion
!     use sim_para
!     implicit none

!                 real*8 :: msd, t, sum_msd, mean_msd, diffusion_const
!                 character(len=50) :: filename
!                 integer :: file_unit
            
!                 npart = sim_para%npart
!                 niter = sim_para%niter
!                 n_calc_av = sim_para%n_calc_av
!                 dt = sim_para%delta_t
            
                
!                 msd = 0.0d0
!                 t = 0.0d0
!                 sum_msd = 0.0d0
            
              
!                 filename = 'msd_vs_time.dat'
!                 open(unit=67, file='diffconst.f90', status='replace')
            
                
!                 do i = 1, niter
!                     do j=1,npart
!                         initial_pos(3*j-2)=pos(3*j)
!                         initial_pos(3*j-1)=pos(3*j)
!                         initial_pos(3*j)=pos(3*j)
                   
!                     call update_pos
            
!                     if (mod(i, n_calc_av) == 0) then
!                         msd = 0.0d0
!                         do j = 1, npart
!                             dx = pos(3*j-2) - initial_pos(3*j-2)
!                             dy = pos(3*j-1) - initial_pos(3*j-1)
!                             dz = pos(3*j) - initial_pos(3*j)
!                             msd = msd + dx**2 + dy**2 + dz**2
!                         end do
!                         msd = msd / dble(npart)
!                         sum_msd = sum_msd + msd
!                         t = t + n_calc_av * dt
            
!                         write(67, '(2E20.10)') t, msd
!                     end if
!                 end do
            
               
!                 close(67)
            
               
!                 mean_msd = sum_msd / dble(niter / n_calc_av)
!                 diffusion_const = mean_msd / (6.0d0 * t)  
            
!                 write(*,*) 'Diffusion constant (D) = ', diffusion_const
            
! end subroutine Diffusion
