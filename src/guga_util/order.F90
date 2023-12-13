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

subroutine ORDER(C,D,N)

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: N
real(kind=wp), intent(inout) :: C(N,N), D(N)
integer(kind=iwp) :: I, IMIN, J, K
real(kind=wp) :: DMIN, TMP

do I=1,N-1
  IMIN = I
  DMIN = D(I)
  do J=I+1,N
    if (D(J) < DMIN) then
      DMIN = D(J)
      IMIN = J
    end if
  end do
  if (I == IMIN) cycle
  D(IMIN) = D(I)
  D(I) = DMIN
  do K=1,N
    TMP = C(K,I)
    C(K,I) = C(K,IMIN)
    C(K,IMIN) = TMP
  end do
end do

return

end subroutine ORDER
