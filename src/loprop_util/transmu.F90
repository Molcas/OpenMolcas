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

subroutine TransMu(SqMu,nDim,TTot,Temp)
! Transform with TTOT the multipole moment integral matrix

use Constants, only: Zero, One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nDim
real(kind=wp), intent(inout) :: SqMu(nDim*nDim)
real(kind=wp), intent(in) :: TTot(nDim*nDim)
real(kind=wp), intent(out) :: Temp(nDim*nDim)

!                                                                      *
!***********************************************************************
!                                                                      *
call DGEMM_('N','N',nDim,nDim,nDim,One,SqMu,nDim,TTot,nDim,Zero,Temp,nDim)
call DGEMM_('T','N',nDim,nDim,nDim,One,Ttot,nDim,Temp,nDim,Zero,SqMu,nDim)
!call RecPrt('Overout',' ',SqMu,nDim,nDim)
!                                                                      *
!***********************************************************************
!                                                                      *
return

end subroutine TransMu
