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

subroutine getw3(wrk,wrksize,lunw3xxxx,nxxxx)
! This routine reconstructs W3(m,e,a,j)xxxx from lunw3xxxx file
! to V1(m,e,a,j) and defines corresponding v1%d and v1%i
! This routine also closes lunw3xxxx file
!
! lunw3xxxx - lun of opened w3xxxx file
! nxxxx     - xxxx identifier
!             1 - aaaa
!             2 - bbbb
!             3 - aabb
!             4 - abba
!             5 - baab
!             6 - bbaa
!
! the structure of lunw3xxxx is:
!
! do syma=1,nsym
!   map of H _a(m,e,j)bbbb (W3(m,e,a,j))
!   skip cycle over a if length of all files is zero
!   do a=1,nvx(syma) [ x is a or b ]
!     if (h1length > 0) write H(m,e,j)
!   end do
! end do

use ccsd_global, only: h1, nsym, nva, nvb, v1
use Constants, only: Zero, One
use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
integer(kind=iwp), intent(in) :: wrksize, nxxxx
real(kind=wp), intent(inout) :: wrk(wrksize)
integer(kind=iwp), intent(_IN_) :: lunw3xxxx
integer(kind=iwp) :: a, aalphayes, aup, h1length, iiv1, post, rc, syma

!0 def aalphayes

if ((nxxxx == 1) .or. (nxxxx == 5) .or. (nxxxx == 6)) then
  aalphayes = 1
else
  aalphayes = 0
end if

!1.1 def maps of V1(m,e,a,j)

if (nxxxx == 1) then
  call grc0(4,0,1,3,3,1,1,post,v1)
else if (nxxxx == 2) then
  call grc0(4,0,2,4,4,2,1,post,v1)
else if (nxxxx == 3) then
  call grc0(4,0,1,3,4,2,1,post,v1)
else if (nxxxx == 4) then
  call grc0(4,0,1,4,4,1,1,post,v1)
else if (nxxxx == 5) then
  call grc0(4,0,2,3,3,2,1,post,v1)
else if (nxxxx == 6) then
  call grc0(4,0,2,4,3,1,1,post,v1)
end if

!1.2 vanish V1
iiv1 = v1%d(0,5)
wrk(v1%pos0:v1%d(iiv1,1)+v1%d(iiv1,2)-1) = Zero

!2 rewind tape lunw3xxxx
call filemanager(2,lunw3xxxx,rc)

!3 loop over symA

do syma=1,nsym

  !3.1 get map of H _a(m,e,j) to %d,%i H1
  call getmap(lunw3xxxx,h1length,h1,rc)

  !3.2 skip cycle over a if length of H1 is 0
  if (h1length == 0) cycle

  !3.3 loop over all a in this symmetry

  if (aalphayes == 1) then
    aup = nva(syma)
  else
    aup = nvb(syma)
  end if

  do a=1,aup

    if (h1length > 0) then

      !3.3.1 read H1 if any
      call rea(lunw3xxxx,h1length,wrk(h1%pos0))

      !3.3.2 insert H1 into V1 for given a and syma
      call add(wrk,wrksize,3,4,1,3,a,0,syma,1,One,h1,syma,v1,1,rc)

    end if

  end do

end do

!4 close lunw3xxxx file
call filemanager(3,lunw3xxxx,rc)

return

end subroutine getw3
