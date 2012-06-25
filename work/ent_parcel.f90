module ent_parcel

CONTAINS

!!***************************************************
SUBROUTINE mixair(P, T1, T2, q1, q2, m1, m2, T_o, md_o, mm_o)

implicit none
!-- memo -----------
! q : specific humidity [kg kg-1]
!
!-- input ----------
double precision                P, T1, T2, q1, q2, m1, m2
!f2py intent(in)                P, T1, T2, q1, q2, m1, m2
!-- output ---------
double precision                T_o, md_o, mm_o
!f2py intent(out)               T_o, md_o, mm_o
!-------------------
double precision,parameter   :: c_d = 1004.0  ! specific heat of dry air [J K-1 kg-1]
double precision,parameter   :: c_m = 1952.0  ! specific heat of moisture [J K-1 kg-1]
!- for calc --------
double precision                md1, md2, mm1, mm2
double precision                A, B, C
double precision                qs_o
!-------------------
md1   = m1 * (1.0 - q1)
md2   = m2 * (1.0 - q2)
mm1   = m1 * q1
mm2   = m2 * q2
!
A     = c_d * (md1 + md2) + c_m * (mm1 + mm2)
B     = md1*c_d + mm1*c_m
C     = md2*c_d + mm2*c_m
!
T_o   = (T1*B + T2*C)/A
md_o  = md1 + md2
mm_o  = mm1 + mm2


END SUBROUTINE mixair
!***************************************************
FUNCTION lcl(rPsfc, rTsfc, rqsfc)
!###########################################################
! original code was obtained from
! http://www1.doshisha.ac.jp/~jmizushi/program/Fortran90/4.2.txt
! modified by: N.Utsumi
! f(x)=x**3+6*x**2+21*x+32
!###########################################################
implicit none
double precision                  rPsfc, rTsfc, rqsfc   ! rPsfc:[Pa]
!f2py intent(in)      rPsfc, rTsfc, rqsfc
double precision                  lcl
!f2py intent(out)     lcl
double precision      dPsfc_hPa, dTsfc, dq
double precision      x, xk, fx
double precision      delta
integer               k
INTEGER,PARAMETER :: KMAX=200
!-------------
!Psfc = 1000   !(hPa)
!Tsfc = 293.15 !(K)
!q    = 0.0087268029 !(kg/kg)
!-------------
dPsfc_hPa = dble(rPsfc)*0.01d0  ! Pa -> hPa
dTsfc = dble(rTsfc)
dq    = dble(rqsfc)
!-------------
x=1000.d0
delta=1.d-10
!-------------
fx=func(x, dPsfc_hPa, dTsfc, dq)
k=0
!WRITE(*,"('x(',i2,')=',1PE15.8,', f(',i2,')=',1PE15.8)") k,x,k,fx
!WRITE(*,*)

DO k=1,KMAX

xk=fnewton(x, dPsfc_hPa, dTsfc, dq)
fx=func(xk, dPsfc_hPa, dTsfc, dq)
!WRITE(*,"('x(',i2,')=',1PE15.8,', f(',i2,')=',1PE15.8)") k,xk,k,fx

    IF(abs(fx)<delta)GOTO 100

x=xk    ! LCL [hPa]

END DO

WRITE(*,*) 'could not solve.'
print *, "Psfc=",dPsfc_hPa
print *, "Tsfc=",dTsfc
print *, "q=",dq
print *, "fx=",fx
if (.not.isnan(x)) then
  STOP
endif

100 CONTINUE

!
if (isnan(x) ) then
  lcl = x    ! lcl = nan
else
  lcl = dble(x) *100.0d0  ! [hPa] -> [Pa]
endif
!-----------------
! for the case: lcl is lower than the surface (RH > 100%)
!-----------------
if (-lcl .lt. -rPsfc) then
  lcl = rPsfc
endif
!-----------------
return
END FUNCTION lcl

!**************************************************************
FUNCTION func(P, Psfc, Tsfc, q)
  implicit none
  double precision      P, Psfc, Tsfc, q
  double precision      f1, f2, func
  double precision      L
!
  double precision :: T0    = 273.16d0  !(K)
  double precision :: e0    = 6.1173d0  !(hPa)
  double precision :: Rv    = 461.7d0   !(J kg^-1 K^-1)
  !double precision :: Lv    = 2.500d6 !(J kg^-1)
  double precision :: epsi  = 0.62185d0 !(-)
  double precision :: Rd    = 287.04d0  !(J kg^-1 K^-1)
  double precision :: Cpd   = 1004.67d0 !(J kg^-1 K^-1)
!
L = dble(cal_latentheat( dble(Tsfc) ))
f1 = (1d0/T0 - Rv/L *log( q * P /( e0*(epsi + q) ) ) )**-1d0
f2 = Tsfc * ( P / Psfc )**(Rd/Cpd)
func = f1 - f2
RETURN
END FUNCTION func

!**************************************************************
FUNCTION fnewton(P, Psfc, Tsfc, q)
  implicit none
  double precision       P, Psfc, Tsfc, q
  double precision       f1, f2, func
  double precision       df1_P, df2_P, df_P
  double precision       fnewton

!
  double precision    L
  double precision :: T0    = 273.16d0  !(K)
  double precision :: e0    = 6.1173d0  !(hPa)
  double precision :: Rv    = 461.7d0   !(J kg^-1 K^-1)
  !double precision :: Lv    = 2.500d6 !(J kg^-1)
  double precision :: epsi  = 0.62185d0 !(-)
  double precision :: Rd    = 287.04d0  !(J kg^-1 K^-1)
  double precision :: Cpd   = 1004.67d0 !(J kg^-1 K^-1)
!
L = dble(cal_latentheat( dble(Tsfc) ))
f1 = (1d0/T0 - Rv/L *log( q * P /( e0*(epsi + q) ) ) )**-1d0
f2 = Tsfc * ( P / Psfc )**(Rd/Cpd)
func = f1 - f2
!
df1_P = 1d0/P * Rv/L *(1/T0 - Rv/L*log( q*P /(e0*(epsi + q)) ) )**-2d0
df2_P = Tsfc* (1d0/Psfc)**(Rd/Cpd) * Rd/Cpd * (P **(Rd/Cpd -1d0))
df_P  = df1_P - df2_P
!
fnewton = P - func / df_P
RETURN
END FUNCTION fnewton

!###########################################################


FUNCTION cal_dqdz(P, T)
  implicit none
!-- input ---------------------------------------
  double precision                 P, T
!-- output --------------------------------------
  double precision                 cal_dqdz
!-- parameter -----------------------------------
  double precision,parameter            ::  g    = 9.80665d0  ! [m/s^2]
  double precision,parameter            ::  epsi = 0.62185d0
  double precision,parameter            ::  cp   = 1004.67d0
  double precision,parameter            ::  Rv   = 461.7d0   !(J kg^-1 K^-1)
  double precision,parameter            ::  Rd   = 287.04d0
  double precision,parameter            ::  Lv   = 2.5d6  ! for vaporization
  double precision,parameter            ::  Ld   = 2.834d6 ! for sublimation
!-- calc ----------------------------------------
  double precision                          AA, BB
  double precision                          qs, dtdz
!------------------------------------------------
qs           = cal_qs(T, P)
dtdz         = dt_dz_moist(P, T)
AA           = Lv * qs *dtdz / (Rv * T * T)
BB           = qs * g / (Rd * T)
cal_dqdz     = AA + BB

RETURN
END FUNCTION cal_dqdz
!***************************************************
FUNCTION dt_dz_moist(P, T)
  implicit none
!-- input ---------------------------------------
  double precision                          P, T 
!f2py intent(in)                            P, T
!-- output --------------------------------------
  double precision                          dt_dz_moist
!f2py intent(out)                           dt_dz_moist
!-- parameter ----------------------------------- 
  double precision,parameter            ::  g    = 9.80665d0  ! [m/s^2]
  double precision,parameter            ::  epsi = 0.62185d0
  double precision,parameter            ::  cp   = 1004.67d0
  double precision,parameter            ::  Rd   = 287.04d0
  double precision,parameter            ::  Lv = 2.5d6  ! for vaporization
  double precision,parameter            ::  Ld = 2.834d6 ! for sublimation

!-- calc ----------------------------------------
  double precision                          qs
  double precision                          AA, BB
!------------------------------------------------
qs           = cal_qs(T, P) 
AA           = 1.0 + Lv * qs / (Rd * T)
BB           = 1.0 + epsi * Lv * Lv * qs / (cp * Rd * T * T)
dT_dZ_moist  = -g/cp * AA /BB
RETURN

END FUNCTION dt_dz_moist

!***************************************************
FUNCTION dt_dp_moist(rP, rT)
  implicit none
  double precision                        rP, rT
!f2py                         rP, rT
  double precision                        res, rqs        ! rP:[Pa], not [hPa]
  double precision                        dT_dP_moist     ! [K/Pa], not [K/hPa]
!f2py                         dT_dP_moist
!** parameters ******
  double precision                          L, a, b, c
  double precision,parameter            ::  epsi = 0.62185d0
  double precision,parameter            ::  cp   = 1004.67d0
  double precision,parameter            ::  Rd   = 287.04d0
  !double precision,parameter            ::  a0 = 0.28571d0
  !double precision,parameter            ::  b0 = 1.347e7d0
  !double precision,parameter            ::  c0 = 2488.4d0
  double precision                          rtemp
!********************
L = cal_latentheat(rT)
a = Rd / cp
b = epsi *(L**2d0)/(cp*Rd)
c = L/cp
rqs = cal_qs(rT, rP)
dT_dP_moist = (a * rT + c *rqs)/( rP *(1d0 + b*rqs/(rT**2d0) ) )
!
RETURN
END FUNCTION dT_dP_moist


!***************************************************
FUNCTION moistadiabat(rP1,rT1, rP2, dP)
  implicit none
  double precision                       rP1, rP2, rT1, dP
!f2py intent(in)             rP1, rP2, rT1, dP
  double precision                       rP, rT
  double precision                       rsign
  double precision                       rTnext, rT2, dT_dP
  double precision                       moistadiabat
!f2py intent(out)            moistadiabat
  integer                    ip, np
!
  double precision                       rtemp
!
if (rP1 .ge. rP2) then
  rsign = 1.0d0
else
  rsign = -1.0d0
end if
np = int( (rP1 - rP2)/dP )
rP = rP1
rT = rT1
do ip = 1,abs(np)
  dT_dP = dT_dP_moist(rP, rT)
  rT = rT - rsign *dT_dP * dP
  rP = rP - rsign *dP
end do
rT2 = rT - rsign * dT_dP_moist( rP, rT ) * abs((rP1 - np*dP) - rP2)
moistadiabat = rT2
RETURN
END FUNCTION moistadiabat
!***************************************************
FUNCTION t1tot2dry(rT1, rP1, rP2)
  implicit none
  double precision                 rT1, rP1, rP2
!f2py intent(in)       rT1, rP1, rP2
  double precision                 t1tot2dry, rT2
!f2py intent(out)      T1toT2dry
  double precision              :: Rd    = 287.04d0  !(J kg^-1 K^-1)
  double precision              :: Cpd   = 1004.67d0 !(J kg^-1 K^-1)
!
rT2 = rT1 * (rP2/rP1)**(Rd/Cpd)
t1tot2dry = rT2
END FUNCTION t1tot2dry
!***************************************************
!***************************************************

FUNCTION cal_rdqdP(rP, rT, dP)
  implicit none
  double precision                   rP, rT, dP
!f2py intent(in)         rP, rT, dP         ! rP : [Pa], not in [hPa]
!--------
  double precision                   rP1, rP2, rT1, rT2, rqs1, rqs2
  double precision                   cal_rdqdP          ! [(g/g)/Pa], not in [(g/g)/hPa]
!f2py intent(out)        cal_rdqdP
!-------------------
rP1 = rP
rP2 = rP - dP
rT1 = rT
rT2 = moistadiabat(rP1, rT1, rP2, dP)
rqs1 = cal_qs(rT1, rP1)
rqs2 = cal_qs(rT2, rP2)
cal_rdqdP = (rqs2 -rqs1)/(rP2 - rP1)
!!
RETURN
END FUNCTION cal_rdqdP
!***************************************************
FUNCTION cal_es(rT)
  implicit none
  double precision                         rT
!f2py intent(in)                           rT 
  double precision                         cal_es     ![Pa]
!f2py intent(out)                          cal_es  
!
  double precision                          L
  double precision,parameter            ::  rT0 = 273.16
  double precision,parameter            ::  res0= 611.73 ![Pa]
  !double precision,parameter            ::  Lv  = 2.5e6  ![J kg-1]
  double precision,parameter            ::  Rv  = 461.7 ![J K-1 kg -1]
!
L = cal_latentheat(rT)
cal_es = res0 * exp( L/Rv *(1.0/rT0 - 1.0/rT))
RETURN
END FUNCTION cal_es
!***************************************************
FUNCTION cal_latentheat(rT)
  implicit none

  double precision                  rT
!f2py intent(in)                    rT
  double precision,parameter     :: Lv = 2.5e6  ! for vaporization
  double precision,parameter     :: Ld = 2.834e6 ! for sublimation
  double precision,parameter     :: rTliq = 273.15  !   0 deg.C
  double precision,parameter     :: rTice = 250.15   ! -23 deg.C
  double precision               cal_latentheat
!f2py intent(out)                cal_latentheat
!
if ( rT .ge. rTliq) then
  cal_latentheat = Lv
else if ( rT .le. rTice ) then
  cal_latentheat = Ld
else
  cal_latentheat = ((rT - rTice)*Lv + (rTliq - rT)*Ld)/(rTliq - rTice)
end if
RETURN
END FUNCTION cal_latentheat
!***************************************************
FUNCTION cal_qs(rT, rP)
  implicit none
!---------------------------------------
! calculate specific humidity [kg kg-1]
!---------------------------------------
  double precision                 rT, rP
!f2py intent(in)       rT, rP
  double precision                 res
  double precision                 cal_qs
!f2py intent(out)      cal_qs
  double precision,parameter    :: repsi = 0.62185d0
!
res = cal_es(rT)
cal_qs = repsi * res / (rP - (1.0-repsi)* res)
RETURN
END FUNCTION cal_qs
!***************************************************
FUNCTION cal_rmixs(P, T)
  implicit none
!-- in -----------------
  double precision                P, T
!f2py intent(in)                  P, T
!-- out ----------------
  double precision                cal_rmixs
!f2py intent(out)                 cal_rmixs
!-- parameter ----------
  double precision,parameter   :: epsi = 0.62185d0
!-- calc ---------------
  double precision                es
!-----------------------
es        = cal_es(T)
cal_rmixs = epsi * es / (P - es)
RETURN
END FUNCTION cal_rmixs
!***************************************************
SUBROUTINE evap_cool(P, T, mliq, mm, md, T_new, mevap)
!###########################################################
! original code was obtained from
! http://www1.doshisha.ac.jp/~jmizushi/program/Fortran90/4.2.txt
! modified by: N.Utsumi
!###########################################################
implicit none
double precision      P, T, mliq, mm, md
!f2py intent(in)      P, T, mliq, mm, md
double precision      T_new, mevap
!f2py intent(out)     T_new, mevap
double precision      x, xk, fx
double precision      delta
!--- calc ------------
double precision      Lv, q, qs, RH
double precision      T_temp, mm_temp, md_temp, q_temp, qs_temp
integer               k
!--- para ------------
INTEGER,PARAMETER :: KMAX=10000000
double precision :: c_d   = 1004.0  ! specific heat of dry air [J K-1 kg-1]
double precision :: c_m   = 1952.0  ! specific heat of moisture [J K-1 kg-1]
!-------------
x=T-0.3
delta=1.d-7
!-------------
Lv = cal_latentheat( T )
!-------------
fx = func_evapcool(x, P, T, mm, md)
k=0
!WRITE(*,"('x(',i2,')=',1PE15.8,', f(',i2,')=',1PE15.8)") k,x,k,fx
!WRITE(*,*)
!-----------------------------------------
q   = mm / (mm + md)
qs  = cal_qs(T, P)
RH  = q / qs
if (RH <100.0) then
  T_temp   = T - mliq * Lv / (mm*c_m + md* c_d)
  q_temp   = (mm + mliq)/(md + mm + mliq)
  qs_temp  = cal_qs(T_temp, P)
  if (q_temp <= qs_temp) then
    T_new  = T_temp
    mevap  = mliq
    !print *, "q_temp, qs_temp", q_temp, qs_temp
    return
  end if
end if
!-----------------------------------------
DO k=1,KMAX

xk =xnew_evapcool(x, P, T, mm, md)
fx =func_evapcool(xk, P, T, mm, md)

!print *,k, xk,fx

!WRITE(*,"('x(',i2,')=',1PE15.8,', f(',i2,')=',1PE15.8)") k,xk,k,fx

    IF(abs(fx)<delta)GOTO 100

x=xk    ! T_new [K]

END DO

100 CONTINUE

T_new = x
mevap = (mm * c_m + md * c_d)* (T - T_new)/Lv

return
END SUBROUTINE evap_cool

!**************************************************************
FUNCTION func_evapcool(T2, P, T1, mm, md)
  implicit none
  double precision      T2, P, T1, mm, md
  double precision      f1, f2, func_evapcool, es2
  double precision      Lv
!
  double precision :: T0    = 273.16d0  !(K)
  double precision :: e0    = 6.1173d0  !(hPa)
  double precision :: Rv    = 461.7d0   !(J kg^-1 K^-1)
  !double precision :: Lv    = 2.500d6 !(J kg^-1)
  double precision :: epsi  = 0.62185d0 !(-)
  double precision :: Rd    = 287.04d0  !(J kg^-1 K^-1)
  double precision :: Cpd   = 1004.67d0 !(J kg^-1 K^-1)
  double precision :: c_d   = 1004.0  ! specific heat of dry air [J K-1 kg-1]
  double precision :: c_m   = 1952.0  ! specific heat of moisture [J K-1 kg-1]


!
Lv = cal_latentheat( T2 )
es2    = cal_es(T2)
f1 = md * epsi * es2 / (P -es2)
f2 = mm + (T1 - T2)* (mm*c_m + md*c_d) / Lv
func_evapcool = f1 - f2
RETURN
END FUNCTION func_evapcool

!**************************************************************
FUNCTION xnew_evapcool(T2, P, T1, mm, md)
  implicit none
  double precision       T2, P, T1, mm, md
  double precision       xnew_evapcool
!--- calc ---------------
  double precision       es2, f1, f2, func, des_dt2
  double precision       df1_T2, df2_T2, df_T2
  double precision       Lv
!
  double precision :: Rv    = 461.7d0   !(J kg^-1 K^-1)
  !double precision :: Lv    = 2.500d6 !(J kg^-1)
  double precision :: epsi  = 0.62185d0 !(-)
  double precision :: Rd    = 287.04d0  !(J kg^-1 K^-1)
  !double precision :: Cpd   = 1004.67d0 !(J kg^-1 K^-1)
  double precision :: c_d   = 1004.0  ! specific heat of dry air [J K-1 kg-1]
  double precision :: c_m   = 1952.0  ! specific heat of moisture [J K-1 kg-1]
!-------------------
Lv     = dble(cal_latentheat( dble(T1) ))
es2    = cal_es(T2)
!
func   = func_evapcool(T2, P, T1, mm, md)
!
des_dt2 = Lv * es2  / (Rv * T2 * T2) 
df1_T2 = md * epsi * P * des_dt2 / (P - es2)
df2_T2 = -(mm*c_m + md*c_d)/Lv

df_T2   = df1_T2 - df2_T2
!
xnew_evapcool = T2 - func / df_T2
RETURN
END FUNCTION xnew_evapcool

!***************************************************


!***************************************************
!***************************************************
!***************************************************
!***************************************************


END MODULE ent_parcel
