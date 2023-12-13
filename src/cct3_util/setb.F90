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

subroutine setb(wrk,wrksize,a,b,factor)
! this routine does
! B = factor . A
!
! a      - A (I)
! b      - B (I)
! factor - numerical factor (I)
!
! mediate B must have defined maps, and they must be
! of identical type as those for A. If they are not defined,
! use grc0 before setb
!
! N.B. this routine should be done using matrix operations

use CCT3_global, only: Map_Type
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: wrksize
real(kind=wp), intent(inout) :: wrk(wrksize)
type(Map_Type), intent(in) :: a, b
real(kind=wp), intent(in) :: factor
integer(kind=iwp) :: length, nhelp, posa0, posb0

!1 def the length of the mediate
nhelp = a%d(0,5)
length = a%d(nhelp,1)+a%d(nhelp,2)-a%d(1,1)
if (length == 0) return

!2 def initial positions
posa0 = a%d(1,1)
posb0 = b%d(1,1)

!3 set B=f.A
wrk(posb0:posb0+length-1) = factor*wrk(posa0:posa0+length-1)

return

end subroutine setb
