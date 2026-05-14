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

subroutine multy0(wrk,wrksize,mvec,ix,y,key)
! This routine realizes multiplying according mvec
! for Y=A*B
! N.B. if key=0, Y file is not vanished (i.e. can be used for
! adding to some existing file)

use ccsd_global, only: Map_Type
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: wrksize, mvec(4096,7), ix, key
real(kind=wp), intent(inout) :: wrk(wrksize)
type(Map_Type), intent(in) :: y
integer(kind=iwp) :: iix, iy, nhelp1, nhelp2, nhelp3, nhelp4, nhelp5

!1 set C=0

if (key == 1) then

  ! Y vector must be vanished

  do iy=1,y%d(0,5)
    nhelp1 = y%d(iy,1)
    nhelp2 = nhelp1+y%d(iy,2)-1
    wrk(nhelp1:nhelp2) = Zero
  end do

end if

!2 Y=Y+A*B

if (ix == 0) return

do iix=1,ix

  ! skip this summation if yes/no=0
  if (mvec(iix,1) == 0) cycle

  ! realize individual summation

  ! def positions of A,B,Y
  nhelp1 = mvec(iix,2)
  nhelp2 = mvec(iix,3)
  nhelp3 = mvec(iix,4)

  ! def rowA(rowY), colA(sum)
  nhelp4 = mvec(iix,5)
  nhelp5 = mvec(iix,6)

  call mv0v1a3u(nhelp4,nhelp5,nhelp5,nhelp4,nhelp4,nhelp5,1,1,wrk(nhelp1),wrk(nhelp2),wrk(nhelp3))

end do

return

end subroutine multy0
