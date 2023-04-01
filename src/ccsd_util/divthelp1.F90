!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************

subroutine divthelp1(t1,dima,dimi,dp)
! this routine does:
! t1(a,i) = t1(a,i)/(dp(i)-dp(a))
!
! t1   - T1 matrix (I/O
! dima - v dimension of T1 (I)
! dimi - o dimension of T1 (I)
! dp   - diagonal part of Fok (I)
!
! N.B. Since for T1 i and a are of the same spin, there is no reason
! to specify spin of dp. It must be automatically the same spin as i and a.

integer dima, dimi
real*8 t1(1:dima,1:dimi)
real*8 dp(*)
! help variables
integer a, i
real*8 dpi, den

do i=1,dimi
  dpi = dp(i)
  do a=1,dima
    ! t1(a,i)=t1(a,i)/(dpi-dp(dimi+a))

    den = dpi-dp(dimi+a)
    if (abs(den) < 1.0d-7) then
      if (abs(t1(a,i)) > 1.0d-10) then
        t1(a,i) = t1(a,i)/den
      end if
    else
      t1(a,i) = t1(a,i)/den
    end if

  end do
end do

return

end subroutine divthelp1
