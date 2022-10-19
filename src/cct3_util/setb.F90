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

subroutine setb(wrk,wrksize,mapda,mapdb,factor)
! this routine does
! B = factor . A
!
! mapda  - direct map of A (I)
! mapdb  - direct map of B (I)
! factor - numerical factor (I)
!
! mediate B must have defined maps, and they must be
! of identical type as those for A. If they are not defined,
! use grc0 before setb
!
! N.B. this routine should be done using matrix operations

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: wrksize, mapda(0:512,6), mapdb(0:512,6)
real(kind=wp) :: wrk(wrksize), factor
integer(kind=iwp) :: length, nhelp, posa0, posb0

!1 def the length of the mediate
nhelp = mapda(0,5)
length = mapda(nhelp,1)+mapda(nhelp,2)-mapda(1,1)
if (length == 0) return

!2 def initial positions
posa0 = mapda(1,1)
posb0 = mapdb(1,1)

!3 set B=f.A
do nhelp=0,length-1
  wrk(posb0+nhelp) = factor*wrk(posa0+nhelp)
end do

return

end subroutine setb
