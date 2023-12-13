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

subroutine BestMatch(nConstr,nOrb,Occ,Match,ldim)

use Constants, only: Two
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nConstr, nOrb, ldim
real(kind=wp), intent(inout) :: Occ(nOrb)
integer(kind=iwp), intent(out) :: Match(2,ldim)
integer(kind=iwp) :: i, j, jConstr, k, l
real(kind=wp) :: delta, delta0, xOcc

jConstr = 1

do

  delta0 = Two
  do i=1,nOrb
    do j=1,i-1
      xOcc = Occ(i)+Occ(j)
      delta = abs(Two-xOcc)
      if (delta < delta0) then
        delta0 = delta
        if (Occ(i) > Occ(j)) then
          Match(1,jConstr) = i
          Match(2,jConstr) = j
        else
          Match(1,jConstr) = j
          Match(2,jConstr) = i
        end if
      end if
    end do
  end do

  if (jConstr >= nConstr) exit

  ! note: Occ array destroyed here
  k = Match(1,jConstr)
  Occ(k) = -42.0_wp
  l = Match(2,jConstr)
  Occ(l) = -42.0_wp
  jConstr = jConstr+1

end do

return

end subroutine BestMatch
