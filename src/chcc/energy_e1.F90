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

subroutine Energy_E1(T1n,Fvo,no,nv,E1)
! this routine does:
! E1 = sum(a,i) T1n(a,i) . Fvo(a,i)
!
! calculate E1 component of energy

use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: no, nv
real(kind=wp), intent(in) :: T1n(nv*no), Fvo(nv*no)
real(kind=wp), intent(out) :: e1
integer(kind=iwp) :: dim_1

e1 = Zero
dim_1 = nv*no
call mr0u3wt(dim_1,dim_1,dim_1,1,1,T1n,Fvo,e1)

return

end subroutine Energy_E1
