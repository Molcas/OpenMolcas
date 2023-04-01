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
! to V1(m,e,a,j) and defines corresponding mapdv1 and mapiv1
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

#include "ccsd1.fh"
#include "ccsd2.fh"
#include "wrk.fh"
integer lunw3xxxx, nxxxx
! help variables
integer rc, syma, a, h1length, posst, aup, aalphayes, iiv1, v1length

!0 def aalphayes

if ((nxxxx == 1) .or. (nxxxx == 5) .or. (nxxxx == 6)) then
  aalphayes = 1
else
  aalphayes = 0
end if

!1.1 def maps of V1(m,e,a,j)

if (nxxxx == 1) then
  call grc0(4,0,1,3,3,1,1,possv10,posst,mapdv1,mapiv1)
else if (nxxxx == 2) then
  call grc0(4,0,2,4,4,2,1,possv10,posst,mapdv1,mapiv1)
else if (nxxxx == 3) then
  call grc0(4,0,1,3,4,2,1,possv10,posst,mapdv1,mapiv1)
else if (nxxxx == 4) then
  call grc0(4,0,1,4,4,1,1,possv10,posst,mapdv1,mapiv1)
else if (nxxxx == 5) then
  call grc0(4,0,2,3,3,2,1,possv10,posst,mapdv1,mapiv1)
else if (nxxxx == 6) then
  call grc0(4,0,2,4,3,1,1,possv10,posst,mapdv1,mapiv1)
end if

!1.2 vanish V1
iiv1 = mapdv1(0,5)
v1length = mapdv1(iiv1,1)+mapdv1(iiv1,2)-possv10
call mv0zero(v1length,v1length,wrk(possv10))

!2 rewind tape lunw3xxxx
call filemanager(2,lunw3xxxx,rc)

!3 loop over symA

do syma=1,nsym

  !3.1 get map of H _a(m,e,j) to mapd,i H1
  call getmap(lunw3xxxx,possh10,h1length,mapdh1,mapih1,rc)

  !3.2 skip cycle over a if length of H1 is 0
  if (h1length == 0) goto 3000

  !3.3 loop over all a in this symmetry

  if (aalphayes == 1) then
    aup = nva(syma)
  else
    aup = nvb(syma)
  end if

  do a=1,aup

    if (h1length > 0) then

      !3.3.1 read H1 if any
      call rea(lunw3xxxx,h1length,wrk(possh10))

      !3.3.2 insert H1 into V1 for given a and syma
      call add(wrk,wrksize,3,4,1,3,a,0,syma,1,1.0d0,mapdh1,syma,mapdv1,mapiv1,1,rc)

    end if

  end do

  3000 continue
end do

!4 close lunw3xxxx file
call filemanager(3,lunw3xxxx,rc)

return

end subroutine getw3
