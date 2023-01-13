!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 2022, Jie J. Bao                                       *
!***********************************************************************
! ****************************************************************
! history:                                                       *
! Jie J. Bao, on Jul 01, 2022, created this file.                *
! ****************************************************************

subroutine GetDiagScr(nScr,Mat,EigVal,nDim)

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(out) :: nScr
integer(kind=iwp), intent(in) :: nDim
real(kind=wp), intent(inout) :: Mat(nDim**2)
real(kind=wp), intent(out) :: EigVal(nDim)
integer(kind=iwp) :: INFO
real(kind=wp) :: Scr(2)

call DSYEV_('V','U',nDim,Mat,nDim,EigVal,Scr,-1,INFO)
NScr = int(Scr(1))

return

end subroutine GetDiagScr

