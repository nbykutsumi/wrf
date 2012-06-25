MODULE mk_sounding

   !USE module_domain
   !USE module_io_domain
   !USE module_state_description
   USE module_model_constants
   !USE module_bc
   !USE module_timing
   !USE module_configure
   USE module_init_utilities
   USE met_funcs
!#ifdef DM_PARALLEL
!   USE module_dm
!#endif


CONTAINS

!*****************************************************************
   SUBROUTINE mk_lqv( p_surf_in, th_surf_in, rh_surf_in, zk_in, theta_in, qv_in, rh_target, nl_in, qv_in_top)
   implicit none
   ! input data
   integer                        nl_in
!f2py intent(in)               :: nl_in
   real                           p_surf_in, th_surf_in, rh_surf_in
!f2py intent(in)               :: p_surf_in, th_surf_in, rh_surf_in
   real,dimension(nl_in)       :: zk_in, theta_in, qv_in
!f2py intent(in)               :: zk_in, theta_in, qv_in
   real                           rh_target
!f2py intent(in)                  rh_target
!---- output ---
   real                           qv_in_top
!f2py intent(out)              :: qv_in_top
!
!--- staggered layer -------------
   real,dimension(nl_in-1)     :: ptot_stag
   real,dimension(nl_in-1)     :: p_stag
   real,dimension(nl_in-1)     :: pb_stag       ! base state pressure
   real,dimension(nl_in-1)     :: zk_stag
!--------------------------------
   real,dimension(nl_in)       :: qv_in_temp
   real                           qv_surf_in, t_surf_in
!---- top value -----------------
   real                           qv_in_top_small, qv_in_top_large
   real                           th_in_top, t_in_top, rh_in_top
   real                           p_top_stag
   real                           p_top_full
!---- top value -----------------
   real                           drh
   integer                        esflag
   integer                        i_temp
   integer                        nl_max
   real                           r
   !-------------------------------------------------------
   integer,parameter           :: iter_max = 1000
   real,parameter              :: drh_thres = 1.0e-6
   !-------------------------------------------------------
   nl_max                      = nl_in
   r                           = r_d
   esflag                      = 1
   !-- qv_surf --------------------------------------------
   t_surf_in     = th_surf_in * (p1000mb / p_surf_in) ** (-r/cp) !(K)
   qv_surf_in    = cal_rmixs(p_surf_in, t_surf_in, esflag)*rh_surf_in
   !-------------------------------------------------------
   drh        = 9999.0
   i_temp     = 0
   !------------
   qv_in_top        = cal_rmixs(p_surf_in, t_surf_in, esflag)
   qv_in_top_large  = qv_in_top*1.5
   qv_in_top_small  = 0.0
   !------------

   th_in_top          = theta_in(nl_in)
   do while (abs(drh) > drh_thres)
     i_temp           = i_temp + 1
     !
     if (i_temp > iter_max) then
       print *,"too much, drh=",1, drh
       exit
     end if
     !
     qv_in_temp          = qv_in
     qv_in_temp(nl_in)   = qv_in_top


     CALL init_domain_rk(p_surf_in, th_surf_in, qv_surf_in, zk_in, theta_in, qv_in_temp, nl_max, nl_in, p_stag, pb_stag, zk_stag)

     !--- interpolate p(staggerd) -> full-coordinate  -----
     ptot_stag           = p_stag + pb_stag
     if (nl_in .eq. 2) then
       p_top_full        = ptot_stag(nl_in-1)
     else
       p_top_full = interp_0( ptot_stag, zk_stag, zk_in(nl_in), nl_in-1)
     endif
     !-----------------------------------------------------

     t_in_top           = th_in_top * (p1000mb / p_top_full) ** (-r/cp) !(K)

     rh_in_top          = qv_in_top / cal_rmixs(p_top_full, t_in_top, esflag)

     drh             = rh_in_top - rh_target

     if (drh .gt. 0.0)then
       qv_in_top_large  = qv_in_top
       qv_in_top        = 0.5 * (qv_in_top_large + qv_in_top_small)
     else if (drh .lt. 0.0) then
       qv_in_top_small  = qv_in_top
       qv_in_top        = 0.5 * (qv_in_top_large + qv_in_top_small)
     end if
   end do
     print *,"*******************************************" 
     print *, i_temp, "p_top, rh_top, qv_in_top", p_top_full, rh_in_top, qv_in_top
     print *,"*******************************************" 



   END SUBROUTINE

!*****************************************************************
   SUBROUTINE init_domain_rk(p_surf_in, th_surf_in, qv_surf_in, zk, theta, qv, nl_max, nl_in, p, pb, zk_stag)


   IMPLICIT NONE
!----------------------------------------------------------
   !  Input data.
   integer                        nl_max, nl_in
!f2py intent(in)               :: nl_max, nl_in
   real                           p_surf_in, th_surf_in, qv_surf_in
!f2py intent(in)               :: p_surf_in, th_surf_in, qv_surf_in
   real,dimension(nl_max)      :: zk, theta, qv
!f2py intent(in)               :: zk, theta, qv
!----- staggered layers --------------------
   real,dimension(nl_max-1)    :: p
!f2py intent(out)              :: p
   real,dimension(nl_max-1)    :: pb       ! base state pressure
!f2py intent(out)              :: pb
   real,dimension(nl_max-1)    :: zk_stag
!f2py intent(out)              :: zk_stag

!----------------------------------------------------------
!
   !  Local data
   INTEGER                             ::                       &
                                  ids, ide, jds, jde, kds, kde, &
                                  ims, ime, jms, jme, kms, kme, &
                                  its, ite, jts, jte, kts, kte, &
                                  i, j, k

   ! Local data

   !INTEGER, PARAMETER :: nl_max = 1000
   !REAL, DIMENSION(nl_max) :: zk, p_in, theta, rho, u, v, qv, pd_in
   REAL, DIMENSION(nl_max) :: p_in, pd_in
!f2py intent(out)          :: p_in, pd_in
   REAL, DIMENSION(nl_max) :: rho


   INTEGER :: icm,jcm, ii, im1, jj, jm1, loop, error, fid, nxc, nyc
   REAL    :: u_mean,v_mean, f0, p_surf, p_level, qvf, z_at_v, z_at_u
   REAL    :: z_scale, xrad, yrad, zrad, rad, delt, cof1, cof2
!   REAL, EXTERNAL :: interp_0
   REAL    :: hm
   REAL    :: pi

!  stuff from original initialization that has been dropped from the Registry 
   REAL    :: vnu, xnu, xnus, dinit0, cbh, p0_temp, t0_temp, zd, zt
   REAL    :: qvf1, qvf2, pd_surf
   INTEGER :: it
   real :: thtmp, ptmp, temp(3)

   LOGICAL :: stretch_grid, dry_sounding

  INTEGER :: xs , xe , ys , ye
   REAL :: mtn_ht
!  For LES, add randx
   real :: randx
!********************************************************
   real,dimension(nl_max)  ::  znw, rdnw, dn, rdn, fnp, fnm
   real                        cf1, cf2, cf3, cfn, cfn1
   real                        rdx, dx, rdy, dy
   real                        p_top
   real                        ht       ! terrain height ?
   real,dimension(nl_max)  ::  phb      ! base-state geopotential
   real,dimension(nl_max)  ::  ph0
   real,dimension(nl_max)  ::  ph_1, ph_2
   real                        mub      ! base-state dry air mass in column, p_surf - p_top
   real,dimension(nl_max)  ::  t_init   ! initial potential temperature
   real,dimension(nl_max)  ::  alb      ! inverse base density (m3 kg-1)
   real,dimension(nl_max)  ::  znu      ! eta values on half (mass) levels
   real,dimension(nl_max)  ::  dnw      ! d(eta) values between full (w) levels
   real                        mu0, mu_1, mu_2
   real,dimension(nl_max)  ::  moist    ! d(eta) values between full (w) levels
   real,dimension(nl_max)  ::  alt      ! inverse density (m3 kg-1)
   real,dimension(nl_max)  ::  al       ! inverse density (m3 kg-1)
   real,dimension(nl_max)  ::  t_base, t_1, t_2
   real,dimension(nl_max)  ::  h_diabatic
   real,dimension(nl_max)  ::  qv_base
   real                        tmn, tsk
!--------------------------------------------------------
   logical                     wrf_dm_on_monitor
   real                        ztop
!********************************************************
   wrf_dm_on_monitor        = .true.
!********************************************************
   kts                      = 1
   kte                      = nl_max
   kde                      = nl_max
   ztop                     = zk(nl_max)
!------------------------------------------------------
!  stretch_grid = .true.
!  FOR LES, set stretch to false
   stretch_grid = .false.
   delt = 3.
!   z_scale = .50
   z_scale = .40
   pi = 2.*asin(1.0)
   !write(6,*) ' pi is ',pi
   nxc = (ide-ids)/2
   nyc = (jde-jds)/2
!------------------------------------------------------
!  for LES, include Coriolis force
   !f        = 1.e-4 

! set up the grid

   IF (stretch_grid) THEN ! exponential stretch for eta (nearly constant dz)
     DO k=1, kde
      znw(k) = (exp(-(k-1)/float(kde-1)/z_scale) - exp(-1./z_scale))/ &
                                (1.-exp(-1./z_scale))
     ENDDO
   ELSE
     DO k=1, kde
      znw(k) = 1. - float(k-1)/float(kde-1)
     ENDDO
   ENDIF

   DO k=1, kde-1
    dnw(k) = znw(k+1) - znw(k)
    rdnw(k) = 1./dnw(k)
    znu(k) = 0.5*(znw(k+1)+znw(k))
   ENDDO
   DO k=2, kde-1
    dn(k) = 0.5*(dnw(k)+dnw(k-1))
    rdn(k) = 1./dn(k)
    fnp(k) = .5* dnw(k  )/dn(k)
    fnm(k) = .5* dnw(k-1)/dn(k)
   ENDDO

   cof1 = (2.*dn(2)+dn(3))/(dn(2)+dn(3))*dnw(1)/dn(2) 
   cof2 =     dn(2)        /(dn(2)+dn(3))*dnw(1)/dn(3) 
   cf1  = fnp(2) + cof1
   cf2  = fnm(2) - cof1 - cof2
   cf3  = cof2       

   cfn  = (.5*dnw(kde-1)+dn(kde-1))/dn(kde-1)
   cfn1 = -.5*dnw(kde-1)/dn(kde-1)
   rdx = 1./dx
   rdy = 1./dy


!  get the sounding from the ascii sounding file, first get dry sounding and 
!  calculate base state

  dry_sounding = .true.
  IF ( wrf_dm_on_monitor ) THEN
  !write(6,*) ' getting dry sounding for base state '

  
  !CALL get_sounding( zk, p_in, pd_in, theta, rho, u, v, qv, dry_sounding, nl_max, nl_in )
  CALL get_sounding( p_surf_in, th_surf_in, qv_surf_in, zk, theta, qv, p_in, pd_in, rho, dry_sounding, nl_max, nl_in)
  ENDIF

  !write(6,*) ' returned from reading sounding, nl_in is ',nl_in

!  find ptop for the desired ztop (ztop is input from the namelist),
!  and find surface pressure

  p_top = interp_0( p_in, zk, ztop, nl_in )

  ht = 0.
!!
!!  xs=ide/2 -3
!!  xs=ids   -3
!!  xe=xs + 6
!!  ys=jde/2 -3
!!  ye=ys + 6
!!  mtn_ht = 500
!!#ifdef MTN

     ht = mtn_ht * 0.25 * &
               ( 1. + COS ( 2*pi/(xe-xs) * ( i-xs ) + pi ) ) * &
               ( 1. + COS ( 2*pi/(ye-ys) * ( j-ys ) + pi ) )
!!#endif
!!#ifdef EW_RIDGE
     ht = mtn_ht * 0.50 * &
               ( 1. + COS ( 2*pi/(ye-ys) * ( j-ys ) + pi ) )

!!#endif
!!#ifdef NS_RIDGE

     ht = mtn_ht * 0.50 * &
               ( 1. + COS ( 2*pi/(xe-xs) * ( i-xs ) + pi ) )

!!#endif
    phb(1) = g * ht
    ph0(1) = g * ht

    p_surf = interp_0( p_in, zk, phb(1)/g, nl_in )
    mub    = p_surf-p_top
!!
!  this is dry hydrostatic sounding (base state), so given grid%p (coordinate),
!  interp theta (from interp) and compute 1/rho from eqn. of state
!!

    DO K = 1, kte-1
      p_level = znu(k)*(p_surf - p_top) + p_top

      !*****************************************************
      pb(k) = p_level
      !*****************************************************
      t_init(k) = interp_0( theta, p_in, p_level, nl_in ) - t0
      alb(k) = (r_d/p1000mb)*(t_init(k)+t0)*(pb(k)/p1000mb)**cvpm
    ENDDO
!!
!  calc hydrostatic balance (alternatively we could interp the geopotential from the
!  sounding, but this assures that the base state is in exact hydrostatic balance with
!  respect to the model eqns.

    DO k  = 2,kte
      phb(k) = phb(k-1) - dnw(k-1)*mub*alb(k-1)
    ENDDO

!!
  IF ( wrf_dm_on_monitor ) THEN
    !write(6,*) ' base state mub, p_surf is ',mub,mub+p_top

  ENDIF

!  calculate full state for each column - this includes moisture.

  !write(6,*) ' getting moist sounding for full state '
  dry_sounding = .false.
  !CALL get_sounding( zk, p_in, pd_in, theta, rho, u, v, qv, dry_sounding, nl_max, nl_in )

  CALL get_sounding( p_surf_in, th_surf_in, qv_surf_in, zk, theta, qv, p_in, pd_in, rho, dry_sounding, nl_max, nl_in)

!  At this point grid%p_top is already set. find the DRY mass in the column 
!  by interpolating the DRY pressure.  

   pd_surf = interp_0( pd_in, zk, phb(1)/g, nl_in )

!  compute the perturbation mass and the full mass

    mu_1 = pd_surf-p_top - mub
    mu_2 = mu_1
    mu0  = mu_1 + mub
! given the dry pressure and coordinate system, interp the potential
! temperature and qv

    do k=1,kde-1
      p_level = znu(k)*(pd_surf - p_top) + p_top

      moist(k) = interp_0( qv, pd_in, p_level, nl_in )
      t_1(k)          = interp_0( theta, pd_in, p_level, nl_in ) - t0
      t_2(k)          = t_1(k)

    enddo

!  integrate the hydrostatic equation (from the RHS of the bigstep
!  vertical momentum equation) down from the top to get grid%p.
!  first from the top of the model to the top pressure

    k = kte-1  ! top level

    qvf1 = 0.5*(moist(k)+moist(k))
    qvf2 = 1./(1.+qvf1)
    qvf1 = qvf1*qvf2

    !********************************************
    p(k) = - 0.5*(mu_1+qvf1*mub)/rdnw(k)/qvf2
    !********************************************
    qvf = 1. + rvovrd*moist(k)
    alt(k) = (r_d/p1000mb)*(t_1(k)+t0)*qvf* &
                (((p(k)+pb(k))/p1000mb)**cvpm)
    al(k) = alt(k) - alb(k)


!  down the column

    do k=kte-2,1,-1
      qvf1 = 0.5*(moist(k)+moist(k+1))
      qvf2 = 1./(1.+qvf1)
      qvf1 = qvf1*qvf2
    !********************************************
      p(k) = p(k+1) - (mu_1 + qvf1*mub)/qvf2/rdn(k+1)
    !********************************************
      qvf = 1. + rvovrd*moist(k)
      alt(k) = (r_d/p1000mb)*(t_1(k)+t0)*qvf* &
                  (((p(k)+pb(k))/p1000mb)**cvpm)
      al(k) = alt(k) - alb(k)
    enddo

!!
!  this is the hydrostatic equation used in the model after the
!  small timesteps.  In the model, grid%al (inverse density)
!  is computed from the geopotential.


    ph_1(1) = 0.
    DO k  = 2,kte

      ph_1(k) = ph_1(k-1) - (1./rdnw(k-1))*(       &
                   (mub+mu_1)*al(k-1)+ &
                    mu_1*alb(k-1)  )
                                                   
      ph_2(k) = ph_1(k) 
      ph0(k) = ph_1(k) + phb(k)

    ENDDO

!#if 0
!  set a few more things

  DO K = kts, kte-1
    !grid%h_diabatic(i,k,j) = 0.
    h_diabatic(k) = 0.
  ENDDO
!!
  IF ( wrf_dm_on_monitor ) THEN
  DO k=1,kte-1
    t_base(k) = t_1(k)
    qv_base(k) = moist(k)

  ENDDO
  ENDIF

     thtmp   = t_2(1)+t0
     ptmp    = p(1)+pb(1)
     temp(1) = thtmp * (ptmp/p1000mb)**rcp
     thtmp   = t_2(2)+t0
     ptmp    = p(2)+pb(2)
     temp(2) = thtmp * (ptmp/p1000mb)**rcp
     thtmp   = t_2(3)+t0
     ptmp    = p(3)+pb(3)
     temp(3) = thtmp * (ptmp/p1000mb)**rcp

!    For LES-CBL, add 5 degrees to the surface temperature!
!
!    grid%tsk(I,J)=grid%cf1*temp(1)+grid%cf2*temp(2)+grid%cf3*temp(3)
     tsk=cf1*temp(1)+cf2*temp(2)+cf3*temp(3)+5.
     tmn=tsk-0.5

     !--- make zk_ztag -----------
     zk_stag = ph0 / 9.81

     !----------------------------

 END SUBROUTINE init_domain_rk

   SUBROUTINE init_module_initialize
   END SUBROUTINE init_module_initialize

!---------------------------------------------------------------------

!  test driver for get_sounding
!
!      implicit none
!      integer n
!      parameter(n = 1000)
!      real zk(n),p(n),theta(n),rho(n),u(n),v(n),qv(n),pd(n)
!      logical dry
!      integer nl,k
!
!      dry = .false.
!      dry = .true.
!      call get_sounding( zk, p, pd, theta, rho, u, v, qv, dry, n, nl )
!      write(6,*) ' input levels ',nl
!      write(6,*) ' sounding '
!      write(6,*) '  k  height(m)  press (Pa) pd(Pa) theta (K) den(kg/m^3)  u(m/s)     v(m/s)    qv(g/g) '
!      do k=1,nl
!        write(6,'(1x,i3,8(1x,1pe10.3))') k, zk(k), p(k), pd(k), theta(k), rho(k), u(k), v(k), qv(k)
!      enddo
!      end
!
!---------------------------------------------------------------------------

      subroutine get_sounding( p_surf, th_surf, qv_surf,&
                               h_input, th_input, qv_input,&
                               p, p_dry, rho,&
                               dry, nl_max, nl_in )
      implicit none

      integer nl_max, nl_in
      real zk(nl_max), p(nl_max), theta(nl_max), rho(nl_max), &
           u(nl_max), v(nl_max), qv(nl_max), p_dry(nl_max)
      logical dry

      integer n, nl
      !parameter(n=nl_in)
      logical debug
      parameter( debug = .true.)

! input sounding data

      real p_surf, th_surf, qv_surf
      real pi_surf, pi(nl_in)
      real h_input(nl_in), th_input(nl_in), qv_input(nl_in)

! diagnostics

      real rho_surf, p_input(nl_in), rho_input(nl_in)
      real pm_input(nl_in)  !  this are for full moist sounding

! local data

      real r
      parameter (r = r_d)
      integer k, it
      real qvf, qvf1, dz
!--------------------------------------------------------------
      nl        = nl_in
      n         = nl_in

!--------------------------------------------------------------
!  first, read the sounding
      !call read_sounding( p_surf, th_surf, qv_surf, &
      !                    h_input, th_input, qv_input, u_input, v_input,n, nl, debug )

      if(dry) then
       do k=1,nl
         qv_input(k) = 0.
       enddo
      endif

      !if(debug) write(6,*) ' number of input levels = ',nl

        nl_in = nl
        if(nl_in .gt. nl_max ) then
          write(6,*) ' too many levels for input arrays ',nl_in,nl_max
          !call wrf_error_fatal ( ' too many levels for input arrays ' )
        end if

!  compute diagnostics,
!  first, convert qv(g/kg) to qv(g/g)

      !do k=1,nl
      !  qv_input(k) = 0.001*qv_input(k)
      !enddo

      !p_surf = 100.*p_surf  ! convert to pascals
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
          qv(k) = qv_input(k)

        enddo



     !if(debug) then
     ! write(6,*) ' sounding '
     ! write(6,*) '  k  height(m)  press (Pa) pd(Pa) theta (K) den(kg/m^3)  u(m/s)     v(m/s)    qv(g/g) '
     ! do k=1,nl
     !   write(6,'(1x,i3,8(1x,1pe10.3))') k, zk(k), p(k), p_dry(k), theta(k), rho(k), u(k), v(k), qv(k)
     ! enddo

     !end if

      end subroutine get_sounding

!---------------------------------------------------------------------------



!!      subroutine get_sounding( zk, p, p_dry, theta, rho, &
!!                               u, v, qv, dry, nl_max, nl_in )
!!      implicit none
!!
!!      integer nl_max, nl_in
!!      real zk(nl_max), p(nl_max), theta(nl_max), rho(nl_max), &
!!           u(nl_max), v(nl_max), qv(nl_max), p_dry(nl_max)
!!      logical dry
!!
!!      integer n
!!      parameter(n=1000)
!!      logical debug
!!      parameter( debug = .true.)
!!
!!! input sounding data
!!
!!      real p_surf, th_surf, qv_surf
!!      real pi_surf, pi(n)
!!      real h_input(n), th_input(n), qv_input(n), u_input(n), v_input(n)
!!
!!! diagnostics
!!
!!      real rho_surf, p_input(n), rho_input(n)
!!      real pm_input(n)  !  this are for full moist sounding
!!
!!! local data
!!
!!      real r
!!      parameter (r = r_d)
!!      integer k, it, nl
!!      real qvf, qvf1, dz
!!
!!!  first, read the sounding
!!
!!      call read_sounding( p_surf, th_surf, qv_surf, &
!!                          h_input, th_input, qv_input, u_input, v_input,n, nl, debug )
!!
!!      if(dry) then
!!       do k=1,nl
!!         qv_input(k) = 0.
!!       enddo
!!      endif
!!
!!      if(debug) write(6,*) ' number of input levels = ',nl
!!
!!        nl_in = nl
!!        if(nl_in .gt. nl_max ) then
!!          write(6,*) ' too many levels for input arrays ',nl_in,nl_max
!!          !call wrf_error_fatal ( ' too many levels for input arrays ' )
!!        end if
!!
!!!  compute diagnostics,
!!!  first, convert qv(g/kg) to qv(g/g)
!!
!!      do k=1,nl
!!        qv_input(k) = 0.001*qv_input(k)
!!      enddo
!!
!!      p_surf = 100.*p_surf  ! convert to pascals
!!      qvf = 1. + rvovrd*qv_input(1) 
!!      rho_surf = 1./((r/p1000mb)*th_surf*qvf*((p_surf/p1000mb)**cvpm))
!!      pi_surf = (p_surf/p1000mb)**(r/cp)
!!
!!      if(debug) then
!!        write(6,*) ' surface density is ',rho_surf
!!        write(6,*) ' surface pi is      ',pi_surf
!!      end if
!!
!!
!!!  integrate moist sounding hydrostatically, starting from the
!!!  specified surface pressure
!!!  -> first, integrate from surface to lowest level
!!
!!          qvf = 1. + rvovrd*qv_input(1) 
!!          qvf1 = 1. + qv_input(1)
!!          rho_input(1) = rho_surf
!!          dz = h_input(1)
!!          do it=1,10
!!            pm_input(1) = p_surf &
!!                    - 0.5*dz*(rho_surf+rho_input(1))*g*qvf1
!!            rho_input(1) = 1./((r/p1000mb)*th_input(1)*qvf*((pm_input(1)/p1000mb)**cvpm))
!!          enddo
!!
!!! integrate up the column
!!
!!          do k=2,nl
!!            rho_input(k) = rho_input(k-1)
!!            dz = h_input(k)-h_input(k-1)
!!            qvf1 = 0.5*(2.+(qv_input(k-1)+qv_input(k)))
!!            qvf = 1. + rvovrd*qv_input(k)   ! qv is in g/kg here
!! 
!!            do it=1,10
!!              pm_input(k) = pm_input(k-1) &
!!                      - 0.5*dz*(rho_input(k)+rho_input(k-1))*g*qvf1
!!              rho_input(k) = 1./((r/p1000mb)*th_input(k)*qvf*((pm_input(k)/p1000mb)**cvpm))
!!            enddo
!!          enddo
!!
!!!  we have the moist sounding
!!
!!!  next, compute the dry sounding using p at the highest level from the
!!!  moist sounding and integrating down.
!!
!!        p_input(nl) = pm_input(nl)
!!
!!          do k=nl-1,1,-1
!!            dz = h_input(k+1)-h_input(k)
!!            p_input(k) = p_input(k+1) + 0.5*dz*(rho_input(k)+rho_input(k+1))*g
!!          enddo
!!
!!
!!        do k=1,nl
!!
!!          zk(k) = h_input(k)
!!          p(k) = pm_input(k)
!!          p_dry(k) = p_input(k)
!!          theta(k) = th_input(k)
!!          rho(k) = rho_input(k)
!!          u(k) = u_input(k)
!!          v(k) = v_input(k)
!!          qv(k) = qv_input(k)
!!
!!        enddo
!!
!!     !if(debug) then
!!     ! write(6,*) ' sounding '
!!     ! write(6,*) '  k  height(m)  press (Pa) pd(Pa) theta (K) den(kg/m^3)  u(m/s)     v(m/s)    qv(g/g) '
!!     ! do k=1,nl
!!     !   write(6,'(1x,i3,8(1x,1pe10.3))') k, zk(k), p(k), p_dry(k), theta(k), rho(k), u(k), v(k), qv(k)
!!     ! enddo
!!
!!     !end if
!!
!!      end subroutine get_sounding

!-------------------------------------------------------

      subroutine read_sounding( ps,ts,qvs,h,th,qv,u,v,n,nl,debug )
      implicit none
      integer n,nl
      real ps,ts,qvs,h(n),th(n),qv(n),u(n),v(n)
      logical end_of_file
      logical debug

      integer k

      open(unit=10,file='input_sounding',form='formatted',status='old')
      rewind(10)
      read(10,*) ps, ts, qvs
      if(debug) then
        write(6,*) ' input sounding surface parameters '
        write(6,*) ' surface pressure (mb) ',ps
        write(6,*) ' surface pot. temp (K) ',ts
        write(6,*) ' surface mixing ratio (g/kg) ',qvs
      end if

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

      end subroutine read_sounding

END MODULE mk_sounding
