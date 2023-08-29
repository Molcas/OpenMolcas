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

subroutine Cho_P_GetMQ(MQ,l_MQ,LstQSP,nQSP)

use Cholesky, only: Cho_Real_Par
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: l_MQ, nQSP, LstQSP(nQSP)
real(kind=wp), intent(out) :: MQ(l_MQ)
character(len=*), parameter :: SecNam = 'Cho_P_GetMQ'

! In parallel:
! This code only works if MxShpr is set to 1
!      otherwise each node reads a slice
!      of MQ and thus a more sophisticated
!      "synchronization" would be needed
! ------------------------------------------

if (Cho_Real_Par) then
  if (nQSP > 1) call Cho_Quit('Oops! Bug detected in '//SecNam,103)
  MQ(:) = Zero
  call Cho_p_QualSwp()
  call Cho_GetMQ(MQ,l_MQ,LstQSP,nQSP)
  call Cho_p_QualSwp()
  call Cho_GAdGop(MQ,l_MQ,'+')
else
  call Cho_GetMQ(MQ,l_MQ,LstQSP,nQSP)
end if

end subroutine Cho_P_GetMQ
