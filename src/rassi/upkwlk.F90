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

subroutine UPKWLK(N,IPWLK,NWALK,IWALK,ICASE)

use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: N, IPWLK, NWALK, IWALK(*)
integer(kind=iwp), intent(out) :: ICASE(N,NWALK)
integer(kind=iwp) :: I, IPOS, IWORD, J, L, LEND, LSTA, NEXT

! See companion subroutine PKWLK.
IPOS = 0
do I=1,NWALK
  LEND = 0
  do J=1,IPWLK
    LSTA = LEND+1
    LEND = min(LSTA+14,N)
    IPOS = IPOS+1
    IWORD = IWALK(IPOS)
    do L=LSTA,LEND
      NEXT = IWORD/4
      ICASE(L,I) = IWORD-4*NEXT
      IWORD = NEXT
    end do
  end do
end do

end subroutine UPKWLK
