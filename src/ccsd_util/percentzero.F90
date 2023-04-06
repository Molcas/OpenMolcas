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

subroutine percentzero(wrk,wrksize,map,pz)
! this routine tests % of small elements in mediate, decribed by map
!
! map - map type of required mediate (I)

use ccsd_global, only: Map_Type
use Constants, only: One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: wrksize
real(kind=wp), intent(in) :: wrk(wrksize)
type(Map_Type), intent(in) :: map
real(kind=wp), intent(out) :: pz
integer(kind=iwp) :: length, nhelp, nzero, pos
real(kind=wp) :: zerolim

! def length, pos, zerolim

pos = map%d(1,1)
nhelp = map%d(0,5)
length = map%d(nhelp,1)+map%d(nhelp,2)-map%d(1,1)
zerolim = 1.0e-6_wp

if (length > 0) then
  nzero = 0
  do nhelp=pos,pos+length-1
    if (abs(wrk(nhelp)) < zerolim) nzero = nzero+1
  end do
  pz = real(100*nzero,kind=wp)/real(length,kind=wp)
else
  pz = One
end if

return

end subroutine percentzero
