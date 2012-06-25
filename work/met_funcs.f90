MODULE met_funcs

CONTAINS

!***************************************************
FUNCTION cal_rmixs(P, T, flag)
  implicit none
!-----------------------
! calc mixing ratio (kg kg-1)
! P: (Pa)
! T: (K)
!-- in -----------------
  integer             flag      ! 1:water,   -1:ice,   0:mix
!f2py intent(in)      flag
  real                P, T
!f2py intent(in)                  P, T
!-- out ----------------
  real                cal_rmixs
!f2py intent(out)                 cal_rmixs
!-- parameter ----------
  real,parameter   :: epsi = 0.62185e0
!-- calc ---------------
  real                es
!-----------------------
es        = cal_es(T, flag)
cal_rmixs = epsi * es / (P - es)
RETURN
END FUNCTION cal_rmixs
!************************************************************
FUNCTION cal_es(rT, flag)
  implicit none
!-- calc es (Pa), tetens equation ------------
  integer                      flag   ! 1:water ,  -1:ice,   0:mix
!f2py intent(in)               flag
  real                         rT
!f2py intent(in)                           rT
  real                         cal_es     ![Pa]
!f2py intent(out)                          cal_es
!
  real                          Tdeg
  real                          es_liq, es_ice
  real,parameter             :: a_liq = 7.5
  real,parameter             :: b_liq = 237.3
  real,parameter             :: a_ice = 9.5
  real,parameter             :: b_ice = 265.3
!
  real,parameter             :: rTliq = 273.15  !   0 deg.C
  real,parameter             :: rTice = 250.15  ! -23 deg.C
!
Tdeg     =  rT - 273.16

if (flag == 1)then
  cal_es = 6.1078 * 10.0**(a_liq * Tdeg/ (b_liq + Tdeg))
else if (flag == -1)then
  cal_es = 6.1078 * 10.0**(a_ice * Tdeg/ (b_ice + Tdeg))
else if (flag == 0)then
  if ( rT .ge. rTliq) then
    cal_es = 6.1078 * 10.0**(a_liq * Tdeg/ (b_liq + Tdeg))
  else if ( rT .le. rTice ) then
    cal_es = 6.1078 * 10.0**(a_ice * Tdeg/ (b_ice + Tdeg))
  else
    es_liq = 6.1078 * 10.0**(a_liq * Tdeg/ (b_liq + Tdeg))
    es_ice = 6.1078 * 10.0**(a_ice * Tdeg/ (b_ice + Tdeg))
    cal_es = ((rT - rTice)*es_liq + (rTliq - rT)*es_ice)/(rTliq - rTice)
  end if
end if

cal_es     = cal_es * 100.0     ! hPa --> Pa   
RETURN
END FUNCTION cal_es
!************************************************************
FUNCTION cal_latentheat(rT, flag)
  implicit none

  integer               flag       ! flag=1 : water, flag=-1: ice, flag= 0, mix
!f2py intent(in)        flag
  real                  rT
!f2py intent(in)                    rT
  real,parameter     :: Lv = 2.5e6  ! for vaporization
  real,parameter     :: Ld = 2.834e6 ! for sublimation
  real,parameter     :: rTliq = 273.15  !   0 deg.C
  real,parameter     :: rTice = 250.15   ! -23 deg.C
  real               cal_latentheat
!f2py intent(out)                cal_latentheat
!
if (flag == 1) then
  cal_latentheat = Lv
else if (flag == -1) then
  cal_latentheat = Ld
else if (flag == 0) then
  if ( rT .ge. rTliq) then
    cal_latentheat = Lv
  else if ( rT .le. rTice ) then
    cal_latentheat = Ld
  else
    cal_latentheat = ((rT - rTice)*Lv + (rTliq - rT)*Ld)/(rTliq - rTice)
  end if
end if
RETURN
END FUNCTION cal_latentheat
!***************************************************



END MODULE met_funcs
