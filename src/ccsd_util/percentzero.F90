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

subroutine percentzero(wrk,wrksize,mapd,pz)
! this routine tests % of small elements in mediate, decribed by mpd
!
! mapd - direct map of required mediate (I)

#include "wrk.fh"
integer mapd(0:512,1:6)
real*8 pz
! help variables
integer poss, length
integer nhelp, nzero
real*8 zerolim

! def length, poss, zerolim

poss = mapd(1,1)
nhelp = mapd(0,5)
length = mapd(nhelp,1)+mapd(nhelp,2)-mapd(1,1)
zerolim = 1.0d-6

if (length > 0) then
  nzero = 0
  do nhelp=poss,poss+length-1
    if (abs(wrk(nhelp)) < zerolim) nzero = nzero+1
  end do
  pz = dble(100*nzero)/dble(length)
else
  pz = 1.0d0
end if

return

end subroutine percentzero
