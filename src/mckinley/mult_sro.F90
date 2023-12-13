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

subroutine mult_sro(A,nA,C,nC,B,nB,Fact,Res,Tmp)

use Constants, only: Zero, One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nA, nC, nB
real(kind=wp), intent(in) :: A(nA,nC), C(nC,nC), B(nC,nB), Fact
real(kind=wp), intent(inout) :: Res(nA,nB)
real(kind=wp), intent(out) :: Tmp(nA,nC)

call DGEMM_('N','N',nA,nC,nC,One,A,nA,C,nC,Zero,Tmp,nA)
call DGEMM_('N','N',nA,nB,nC,Fact,Tmp,nA,B,nC,One,Res,nA)

return

end subroutine mult_sro
