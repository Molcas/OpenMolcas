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

! This subroutine should be in a module, to avoid explicit interfaces
#ifndef _IN_MODULE_
#error "This file must be compiled inside a module"
#endif

subroutine Cho_P_GetLQ(QVec,l_QVec,LstQSP,nQSP)

use Cholesky, only: Cho_Real_Par
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: l_QVec, nQSP, LstQSP(nQSP)
real(kind=wp), target, intent(out) :: QVec(l_Qvec)
character(len=*), parameter :: SecNam = 'Cho_P_GetLQ'

! In parallel:
! This code only works if MxShpr is set to 1
!      otherwise each node computes a slice
!      of QVec and thus a more sophisticated
!      "synchronization" would be needed

if (Cho_Real_Par) then
  if (nQSP > 1) call Cho_Quit('Oops! Bug detected in '//SecNam,103)
  QVec(:) = Zero
  call Cho_p_QualSwp()
  call Cho_GetLQ(QVec,l_QVec,LstQSP,nQSP)
  call Cho_p_QualSwp()
  call Cho_GAdGOp(QVec,l_QVec,'+') ! sync. array
else
  call Cho_GetLQ(QVec,l_QVec,LstQSP,nQSP)
end if

end subroutine Cho_P_GetLQ
