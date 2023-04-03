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

use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp) :: wrksize, nind, mapd(0:512,6), mapi(8,8,8), imax(8,5)
real(kind=wp) :: wrk(wrksize), rmax(5)
character(len=8) :: text
integer(kind=iwp) :: nhelp1, nhelp2, rc
real(kind=wp) :: scalar

!1 write notice

write(u6,101) text

!2 write 5 maximal amplitudes

write(u6,102)
do nhelp1=1,5
  write(u6,103) (imax(nhelp2,nhelp1),nhelp2=1,8),rmax(nhelp1)
end do

!3 write euclidian norm

!3.1 calc euclidian norm
call multdot(wrk,wrksize,nind,mapd,mapi,1,mapd,mapi,1,scalar,rc)
scalar = sqrt(scalar)

write(u6,104) scalar

write(u6,*)

return

101 format(' Five largest amplitudes of :',a8)
102 format('  SYMA   SYMB   SYMI   SYMJ     A      B      I      J     VALUE')
103 format(8(2x,i3,2x),f15.10)
104 format(' Euclidian norm is :',f17.10)

end subroutine max5h2
