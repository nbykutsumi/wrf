module temp

USE met_funcs
CONTAINS

subroutine subf(a)
  implicit none
  integer           a
!f2py intent(out)   a
!  real              r_d


  a = 1
  print *,cal_es(273.16)

end subroutine subf

end module temp
