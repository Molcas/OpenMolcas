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

subroutine PKWLK(N,IPWLK,NWALK,IWALK,ICASE)
! PURPOSE: PACK THE GUGA STEP NUMBERS INTO THE ARRAY IWALK.
! EACH OF THE NWALK WALKS HAS N STEP NUMBERS, 2 BITS EACH,
! AT MOST 15 TO AN INTEGER ELEMENT OF IWALK, EACH NEW WALK
! IS ALIGNED ON INTEGERS.
! NOTE: Can be used for upper, lower, or total walks, so the
! number of integers used for each walk (IPWLK) is given as
! call parameter.

use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: N, IPWLK, NWALK, ICASE(N,NWALK)
integer(kind=iwp), intent(out) :: IWALK(*)
integer(kind=iwp) :: I, IPOS, IWORD, J, L, LEND, LSTA

IPOS = 0
do I=1,NWALK
  LEND = 0
  do J=1,IPWLK
    LSTA = LEND+1
    LEND = min(LSTA+14,N)
    IPOS = IPOS+1
    IWORD = 0
    do L=LEND,LSTA,-1
      IWORD = 4*IWORD+ICASE(L,I)
    end do
    IWALK(IPOS) = IWORD
  end do
end do

end subroutine PKWLK
