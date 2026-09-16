! program ring
!     implicit none
!     integer, parameter :: nop = 50
!     integer :: niter, i, dx1000, kk, kkp1, kkn1
!     real*8 :: init_t, init_y(nop), init_vel(nop)
!     real*8, parameter :: km = 1.0d0, mass = 1.0d0
!     real*8 :: dt, t, dtby2, ke, pe, zpos, twopi, r, theta(nop)
!     real*8 :: f0(nop), f1(nop), f2(nop), f3(nop), f0v(nop), f1v(nop), f2v(nop), f3v(nop)
!     real*8 :: y1(nop), y2(nop), y3(nop), v1(nop), v2(nop), v3(nop), force(nop)
!     character(len = 30) :: charac_a, charac_e, charac_m, charac_i



! end program ring

program nparticle

    implicit none
    
    real(8)::f0(50), f1(50),f2(50),f3(50)
    real(8)::f0v(50),f1v(50),f2v(50),f3v(50)
    real(8)::x0(50),x1(50),x2(50),x3(50)
    real(8)::v0(50),v1(50),v2(50),v3(50)
    real(8)::E(50),V(50),KE(50)

    real(8)::t0,dt
    integer::n_iter,i,j
    t0=0.0d0 
    
    dt=0.02d0
    n_iter=40/dt
    do i=1,50
        v0(i)=0
        if (i== 1 .or. i== 26) then
          x0(i)=0.8d0
        else
          x0(i)=0.0d0
        end if
        print*, x0(i)
    end do
    
    open(1,file='nparticles.dat')
    do i=1,n_iter

      do j=1,50 
        f0(j)=v0(j) 
        if (j==1) then
          f0v(j)=-2*x0(j)+x0(50)+x0(j+1)
        end if
        if (j==50) then
          f0v(j)=-2*x0(j)+x0(j-1)+x0(1)
        else
          f0v(j)= -2*x0(j)+x0(j-1)+x0(j+1) 
        end if
      end do

      do j=1,50
        x1(j)=x0(j)+(dt/2)*f0(j)
        v1(j)=v0(j)+(dt/2)*f0v(j)
      end do
      
      do j=1,50
        f1(j)=v1(j)
        if (j==1) then
            f1v(j)=-2*x1(j)+x1(50)+x1(j+1)
        end if
        if (j==50) then
            f1v(j)=-2*x1(j)+x1(j-1)+x1(1)
        else
            f1v(j)= -2*x1(j)+x1(j-1)+x1(j+1) 
        end if
      end do

      do j=1,50
        x2(j)=x0(j)+(dt/2)*f1(j)
        v2(j)=v0(j)+(dt/2)*f1v(j)
      end do
      
      do j=1,50
        f2(j)=v2(j)
        if (j==1) then
            f2v(j)=-2*x2(j)+x2(50)+x2(j+1)
        end if
        if (j==50) then
            f2v(j)=-2*(x2(j))+x2(j-1)+x2(1)
        else
            f2v(j)= -2*(x2(j))+ x2(j-1) + x2(j+1) 
        end if
      end do

      do j=1,50
        x3(j)=x0(j)+(dt)*f2(j)
        v3(j)=v0(j)+(dt)*f2v(j)
      end do

      do j=1,50
        f3(j)=v3(j)
        if (j==1) then
            f3v=-2*x3(j)+ x3(50) + x3(j+1)
        end if
        if (j==50) then
            f3v=-2*x3(j)+x3(j-1)+x3(1)
        else
            f3v(j)= -2*x3(j)+ x3(j-1)+ x3(j+1) 
        end if
      end do

      do j=1,50
        x0(j)=x0(j)+(dt/6)*((f0(j))+(2*f1(j))+(2*f2(j))+(f3(j)))
        v0(j)=v0(j)+(dt/6)*((f0v(j))+(2*f1v(j))+(2*f2v(j))+(f3v(j)))
        write(1,*) n_iter,j, x0(j),v0(j),t0
      end do

      t0=t0+dt
    
     
    end do
    
end program nparticle