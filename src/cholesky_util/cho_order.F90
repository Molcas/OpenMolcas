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

subroutine CHO_ORDER(VEC,LVEC,IJOB)
!
! Purpose: sort elements of VEC in
!          IJOB = -1: descending order
!          IJOB =  1: ascending  order
!          (all other IJOB values are ignored)

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: LVEC, IJOB
real(kind=wp), intent(inout) :: VEC(LVEC)
integer(kind=iwp) :: I, IMAX, IMIN, J
real(kind=wp) :: VMAX, VMIN

if (IJOB == -1) then

  do I=1,LVEC-1
    VMAX = VEC(I)
    IMAX = I
    do J=I+1,LVEC
      if (VEC(J) > VMAX) then
        VMAX = VEC(J)
        IMAX = J
      end if
    end do
    if (IMAX /= I) then
      VEC(IMAX) = VEC(I)
      VEC(I) = VMAX
    end if
  end do

else if (IJOB == 1) then

  do I=1,LVEC-1
    VMIN = VEC(I)
    IMIN = I
    do J=I+1,LVEC
      if (VEC(J) < VMIN) then
        VMIN = VEC(J)
        IMIN = J
      end if
    end do
    if (IMIN /= I) then
      VEC(IMIN) = VEC(I)
      VEC(I) = VMIN
    end if
  end do

end if

end subroutine CHO_ORDER
