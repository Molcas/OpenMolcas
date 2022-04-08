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

subroutine getLIMIT(l1,l2,l3,l4,Lanf,Lend)
!bs get the minimum and maximum L-values
!bs of the the coulomb-potential to interact
!bs with l1-l4

use Definitions, only: iwp, u6

implicit none
integer(kind=iwp), intent(in) :: l1, l2, l3, l4
integer(kind=iwp), intent(out) :: Lanf, Lend
integer(kind=iwp) :: lower1, lower2, lsum, lupper1, lupper2

lower1 = abs(l1-l3)
lower2 = abs(l2-l4)
lupper1 = l1+l3
lupper2 = l2+l4
Lanf = max(lower1,lower2)
Lend = min(lupper1,lupper2)
!bs check for parity
lsum = Lanf+l1+l3
if (mod(lsum,2) == 1) Lanf = Lanf+1
lsum = Lend+l1+l3
if (mod(lsum,2) == 1) Lend = Lend-1
!bs check the other parity
lsum = Lanf+l2+l4
if (mod(lsum,2) == 1) then
  write(u6,*) ' error in getLIMIT: '
  write(u6,*) ' parity inconsistency for '
  write(u6,*) 'l1,l2,l3,l4= ',l1,l2,l3,l4
  call Abend()
end if

return

end subroutine getLIMIT
