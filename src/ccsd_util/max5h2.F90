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

subroutine max5h2(wrk,wrksize,nind,mapd,mapi,rmax,imax,text)
! this routine writes:
! a) note
! b) 5 maximal elements with their indices in given vector V
! c) euclidian norm
!
! nind - number of indices in V (I)
! mapd - direct map of V (I)
! mapi - inverse map of V (I)
! rmax - store of maximal values (I)
! imax - store of corr. indices (I)
! text - notice (I)

#include "ccsd1.fh"
#include "wrk.fh"
integer nind
integer mapd(0:512,1:6)
integer mapi(1:8,1:8,1:8)
integer imax(1:8,1:5)
real*8 rmax(1:5)
character*8 text
! help variables
integer nhelp1, nhelp2, rc
real*8 scalar

!1 write notice

write(6,101) text
101 format(' Five largest amplitudes of :',a8)

!2 write 5 maximal amplitudes

write(6,102)
102 format('  SYMA   SYMB   SYMI   SYMJ     A      B      I      J     VALUE')
do nhelp1=1,5
  write(6,103) (imax(nhelp2,nhelp1),nhelp2=1,8),rmax(nhelp1)
103 format(8(2x,i3,2x),f15.10)
end do

!3 write euclidian norm

!3.1 calc euclidian norm
call multdot(wrk,wrksize,nind,mapd,mapi,1,mapd,mapi,1,scalar,rc)
scalar = sqrt(scalar)

write(6,104) scalar
104 format(' Euclidian norm is :',f17.10)

write(6,*)

return

end subroutine max5h2
