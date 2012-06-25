!program sounding
module sounding

  USE module_model_constants
  USE met_funcs
CONTAINS

!*********************************************************************************
      subroutine sounding_rh2q_zlev(esflag, p_surf, th_surf, qv_surf, rh_surf&
                                   , zk, p, theta, rho   &
                                   , u, v, qv, rh, tk, nl_in )


      implicit none

      integer,parameter    :: nl_max=1000
      integer                 esflag
!f2py intent(in)              esflag

      integer                 nl_in
!f2py intent(out)             nl_in
      real p_surf, th_surf, qv_surf, rh_surf
!f2py intent(out)             p_surf, th_surf, qv_surf, rh_surf 
      real                    zk(nl_max), p(nl_max), theta(nl_max), rho(nl_max), &
                              u(nl_max), v(nl_max), qv(nl_max), rh(nl_max), tk(nl_max)
!f2py intent(out)             zk(nl_max), p(nl_max), theta(nl_max), rho(nl_max)
!f2py intent(out)             u(nl_max), v(nl_max), qv(nl_max), rh(nl_max), tk(nl_max)

      !logical dry

      integer n
      parameter(n=1000)
      integer,parameter    :: iter_max  = 200
      real,parameter       :: drh_thres = 0.001

! input sounding data

      real pi_surf, pi(n)
      real h_input(n), th_input(n), qv_input(n), rh_input(n), u_input(n), v_input(n)

! diagnostics

      real rho_surf, p_input(n), rho_input(n)
      real pm_input(n)  !  this are for full moist sounding



! local data

      real r
      !parameter (r = r_d)
      integer k, it, nl
      integer i_temp
      real qvf, qvf1, dz
      real qvs_surf
      real tk_surf
      real rh_temp, rh_temp2, tk_temp, qv_temp, qvs_temp
      real drh, qv_large, qv_small
      r     = r_d


!  first, read the sounding

      call read_sounding_rh_zlev( p_surf, th_surf, rh_surf, &
                          h_input, th_input, rh_input, u_input, v_input,n, nl)

        nl_in = nl

!  compute diagnostics,
!  first, convert qv(g/kg) to qv(g/g)

      !do k=1,nl
      !  qv_input(k) = 0.001*qv_input(k)
      !enddo


      p_surf = 100.*p_surf  ! convert to pascals

      !-----------------------
      ! calc qv_surf

      tk_surf = th_surf * (p1000mb / p_surf)**(-r/cp)  ! (K)
      qvs_surf = cal_rmixs(p_surf, tk_surf, esflag)
      qv_surf  = qvs_surf * rh_surf

      !------------------------
      

      !qvf = 1. + rvovrd*qv_input(1) 
      qvf = 1. + rvovrd*qv_surf

      rho_surf = 1./((r/p1000mb)*th_surf*qvf*((p_surf/p1000mb)**cvpm))
      pi_surf = (p_surf/p1000mb)**(r/cp)


!  integrate moist sounding hydrostatically, starting from the
!  specified surface pressure
!  -> first, integrate from surface to lowest level

       !---------------------------------
       ! estimate qv_input

          qv_input(1)  = qv_surf * rh_input(1) / rh_surf
          qv_temp      = qv_input(1)
          drh          = 9999.0

          tk_temp      = tk_surf * th_input(1)/th_surf
          qvs_temp     = cal_rmixs(p_surf, tk_temp, esflag)
          qv_surf      = qvs_temp * rh_surf
          qv_large     = qvs_temp * 1.5
          qv_small     = 0.0

          i_temp       = 0
          do while (abs(drh) > drh_thres)
            i_temp     = i_temp + 1

            if (i_temp > iter_max) then
              print *,"too much, drh=",1, drh
              exit
            end if

            qvf = 1. + rvovrd*qv_input(1) 
            qvf1 = 1. + qv_input(1)
            rho_input(1) = rho_surf
            dz = h_input(1)

            do it=1,10
              pm_input(1) = p_surf &
                      - 0.5*dz*(rho_surf+rho_input(1))*g*qvf1
              rho_input(1) = 1./((r/p1000mb)*th_input(1)*qvf*((pm_input(1)/p1000mb)**cvpm))
            enddo

            tk(1)    =  th_input(1) * (p1000mb / pm_input(1))**(-r/cp)  ! (K)
            tk(1)    =  tk_temp
            rh_temp  =  qv_input(1) / cal_rmixs(pm_input(1), tk(1), esflag)
            drh      =  rh_temp - rh_input(1) 
            rh(1)    =  rh_temp

            if (drh .gt. 0.0)then
              qv_large      = qv_input(1)
              qv_input(1)   = 0.5 * (qv_large + qv_small)
            else if (drh .lt. 0.0) then
              qv_small      = qv_input(1)
              qv_input(1)   = 0.5 * (qv_large + qv_small)
            end if
          enddo
       !----------------------------------


! integrate up the column

          do k=2,nl
          !do k=2,2

            qv_input(k)  = qv_input(k-1) * rh_input(k) / rh_input(k-1)
            qv_temp      = qv_input(k)
            drh          = 9999.0
  
            tk_temp      = tk(k-1) * th_input(k)/th_input(k-1)
            qvs_temp     = cal_rmixs(pm_input(k-1), tk_temp, esflag)
            qv_input(k)  = qvs_temp * rh_input(k)
            qv_large     = qvs_temp*1.5
            qv_small     = 0.0
  
            i_temp       = 0
            do while (abs(drh) > drh_thres)
              i_temp     = i_temp + 1
              if (i_temp > iter_max) then
                print *,k, "too much,  drh=", drh
                exit
              end if

              !----------
              rho_input(k) = rho_input(k-1)
              dz = h_input(k)-h_input(k-1)
              qvf1 = 0.5*(2.+(qv_input(k-1)+qv_input(k)))
              qvf = 1. + rvovrd*qv_input(k)   ! qv is in g/kg here

  
              qvf = 1. + rvovrd*qv_input(1) 
              qvf1 = 1. + qv_input(1)
              rho_input(1) = rho_surf
              dz = h_input(1)
  
              do it=1,10
                pm_input(k) = pm_input(k-1) &
                        - 0.5*dz*(rho_input(k)+rho_input(k-1))*g*qvf1
                rho_input(k) = 1./((r/p1000mb)*th_input(k)*qvf*((pm_input(k)/p1000mb)**cvpm))
              enddo
              !-----------


              tk(k)  =  th_input(k) * (p1000mb / pm_input(k))**(-r/cp)  ! (K)
              rh_temp  =  qv_input(k) / cal_rmixs(pm_input(k), tk(k), esflag)
              drh      =  rh_temp - rh_input(k) 
              rh(k)    =  rh_temp

              if (drh .gt. 0.0)then
                qv_large    = qv_input(k)
                qv_input(k) = 0.5 * (qv_large + qv_small)
              else if (drh .lt. 0.0) then
                qv_small    = qv_input(k)
                qv_input(k) = 0.5 * (qv_large + qv_small)
              end if
            enddo
            

            
            rh_temp2        = qv_input(k)/ cal_rmixs(pm_input(k), tk(k), esflag)
            print *, k, pm_input(k), tk(k), rh_temp2,  6.11*10.0**(7.5*(( tk(k)-273.16)/(237.3+tk(k) - 273.16))), cal_es(tk(k), esflag)
             
          enddo
!
!!  we have the moist sounding
!
!!  next, compute the dry sounding using p at the highest level from the
!!  moist sounding and integrating down.
!
!        p_input(nl) = pm_input(nl)
!
!          do k=nl-1,1,-1
!            dz = h_input(k+1)-h_input(k)
!            p_input(k) = p_input(k+1) + 0.5*dz*(rho_input(k)+rho_input(k+1))*g
!          enddo
!
!
        do k=1,nl

          zk(k) = h_input(k)
          p(k) = pm_input(k)
          !p_dry(k) = p_input(k)
          theta(k) = th_input(k)
          rho(k) = rho_input(k)
          u(k) = u_input(k)
          v(k) = v_input(k)
          qv(k) = qv_input(k)

        enddo

      end subroutine sounding_rh2q_zlev

!*********************************************************************************
!---------------------------------------------------------------------------
      subroutine sounding_z2p( zk, p, p_dry, theta, rho, &
                               u, v, qv, nl_in )


      implicit none

      integer,parameter    :: nl_max=1000

      integer                 nl_in
!f2py intent(out)             nl_in
      real                    zk(nl_max), p(nl_max), theta(nl_max), rho(nl_max), &
                              u(nl_max), v(nl_max), qv(nl_max), p_dry(nl_max)
!f2py intent(out)             zk(nl_max), p(nl_max), theta(nl_max), rho(nl_max)
!f2py intent(out)             u(nl_max), v(nl_max), qv(nl_max), p_dry(nl_max)

      !logical dry

      integer n
      parameter(n=1000)

! input sounding data

      real p_surf, th_surf, qv_surf
      real pi_surf, pi(n)
      real h_input(n), th_input(n), qv_input(n), u_input(n), v_input(n)

! diagnostics

      real rho_surf, p_input(n), rho_input(n)
      real pm_input(n)  !  this are for full moist sounding



! local data

      real r
      !parameter (r = r_d)
      integer k, it, nl
      real qvf, qvf1, dz
      r     = r_d


!  first, read the sounding

      call read_sounding_zlev( p_surf, th_surf, qv_surf, &
                          h_input, th_input, qv_input, u_input, v_input,n, nl)

      !if(dry) then
      ! do k=1,nl
      !   qv_input(k) = 0.
      ! enddo
      !endif

      !if(debug) write(6,*) ' number of input levels = ',nl

        nl_in = nl
        !if(nl_in .gt. nl_max ) then
        !  write(6,*) ' too many levels for input arrays ',nl_in,nl_max
        !  call wrf_error_fatal ( ' too many levels for input arrays ' )
        !end if

!  compute diagnostics,
!  first, convert qv(g/kg) to qv(g/g)

      do k=1,nl
        qv_input(k) = 0.001*qv_input(k)
      enddo

      p_surf = 100.*p_surf  ! convert to pascals
      qvf = 1. + rvovrd*qv_input(1) 
      rho_surf = 1./((r/p1000mb)*th_surf*qvf*((p_surf/p1000mb)**cvpm))
      pi_surf = (p_surf/p1000mb)**(r/cp)

      !if(debug) then
      !  write(6,*) ' surface density is ',rho_surf
      !  write(6,*) ' surface pi is      ',pi_surf
      !end if


!  integrate moist sounding hydrostatically, starting from the
!  specified surface pressure
!  -> first, integrate from surface to lowest level

          qvf = 1. + rvovrd*qv_input(1) 
          qvf1 = 1. + qv_input(1)
          rho_input(1) = rho_surf
          dz = h_input(1)
          do it=1,10
            pm_input(1) = p_surf &
                    - 0.5*dz*(rho_surf+rho_input(1))*g*qvf1
            rho_input(1) = 1./((r/p1000mb)*th_input(1)*qvf*((pm_input(1)/p1000mb)**cvpm))
          enddo

! integrate up the column

          do k=2,nl
            rho_input(k) = rho_input(k-1)
            dz = h_input(k)-h_input(k-1)
            qvf1 = 0.5*(2.+(qv_input(k-1)+qv_input(k)))
            qvf = 1. + rvovrd*qv_input(k)   ! qv is in g/kg here
 
            do it=1,10
              pm_input(k) = pm_input(k-1) &
                      - 0.5*dz*(rho_input(k)+rho_input(k-1))*g*qvf1
              rho_input(k) = 1./((r/p1000mb)*th_input(k)*qvf*((pm_input(k)/p1000mb)**cvpm))
            enddo
          enddo

!  we have the moist sounding

!  next, compute the dry sounding using p at the highest level from the
!  moist sounding and integrating down.

        p_input(nl) = pm_input(nl)

          do k=nl-1,1,-1
            dz = h_input(k+1)-h_input(k)
            p_input(k) = p_input(k+1) + 0.5*dz*(rho_input(k)+rho_input(k+1))*g
          enddo


        do k=1,nl

          zk(k) = h_input(k)
          p(k) = pm_input(k)
          p_dry(k) = p_input(k)
          theta(k) = th_input(k)
          rho(k) = rho_input(k)
          u(k) = u_input(k)
          v(k) = v_input(k)
          qv(k) = qv_input(k)

        enddo

     !if(debug) then
     ! write(6,*) ' sounding '
     ! write(6,*) '  k  height(m)  press (Pa) pd(Pa) theta (K) den(kg/m^3)  u(m/s)     v(m/s)    qv(g/g) '
     ! do k=1,nl
     !   write(6,'(1x,i3,8(1x,1pe10.3))') k, zk(k), p(k), p_dry(k), theta(k), rho(k), u(k), v(k), qv(k)
     ! enddo

     !end if

      end subroutine sounding_z2p

!-------------------------------------------------------

      subroutine read_sounding_rh_zlev( ps,ts,rhs,h,th,rh,u,v,n,nl)
      implicit none
      integer          n
!f2py intent(in)       n
      integer          nl
!f2py intent(out)      nl

      real             ps,ts,rhs,h(n),th(n),rh(n),u(n),v(n)
!f2py intent(out)      ps,ts,rhs,h(n),th(n),rh(n),u(n),v(n)
      logical end_of_file

      integer k

      open(unit=10,file='input_sounding_rh.txt',form='formatted',status='old')
      rewind(10)
      read(10,*) ps, ts, rhs
      !if(debug) then
      !  write(6,*) ' input sounding surface parameters '
      !  write(6,*) ' surface pressure (mb) ',ps
      !  write(6,*) ' surface pot. temp (K) ',ts
      !  write(6,*) ' surface mixing ratio (g/kg) ',rhs
      !end if

      end_of_file = .false.
      k = 0

      do while (.not. end_of_file)

        read(10,*,end=100) h(k+1), th(k+1), rh(k+1), u(k+1), v(k+1)
        k = k+1
        !if(debug) write(6,'(1x,i3,5(1x,e10.3))') k, h(k), th(k), rh(k), u(k), v(k)
        go to 110
 100    end_of_file = .true.
 110    continue
      enddo

      nl = k

      close(unit=10,status = 'keep')

      end subroutine read_sounding_rh_zlev

!-----------------------------------------------------

      subroutine read_sounding_zlev( ps,ts,qvs,h,th,qv,u,v,n,nl)
      implicit none
      integer          n
!f2py intent(in)       n
      integer          nl
!f2py intent(out)      nl

      real             ps,ts,qvs,h(n),th(n),qv(n),u(n),v(n)
!f2py intent(out)      ps,ts,qvs,h(n),th(n),qv(n),u(n),v(n)
      logical end_of_file

      integer k

      open(unit=10,file='input_sounding',form='formatted',status='old')
      rewind(10)
      read(10,*) ps, ts, qvs
      !if(debug) then
      !  write(6,*) ' input sounding surface parameters '
      !  write(6,*) ' surface pressure (mb) ',ps
      !  write(6,*) ' surface pot. temp (K) ',ts
      !  write(6,*) ' surface mixing ratio (g/kg) ',qvs
      !end if

      end_of_file = .false.
      k = 0

      do while (.not. end_of_file)

        read(10,*,end=100) h(k+1), th(k+1), qv(k+1), u(k+1), v(k+1)
        k = k+1
        !if(debug) write(6,'(1x,i3,5(1x,e10.3))') k, h(k), th(k), qv(k), u(k), v(k)
        go to 110
 100    end_of_file = .true.
 110    continue
      enddo

      nl = k

      close(unit=10,status = 'keep')

      end subroutine read_sounding_zlev

!-----------------------------------------------------

END MODULE SOUNDING
!end program sounding
