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
! Copyright (C) 2020, Roland Lindh                                     *
!***********************************************************************

subroutine Trans_K(X,Y,nInter,nIter)

use kriging_mod, only: layer_U
use Constants, only: Zero, One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nInter, nIter
real(kind=wp), intent(in) :: X(nInter,nIter)
real(kind=wp), intent(out) :: Y(nInter,nIter)

!call RecPrt('layer_U',' ',layer_U,nInter,nInter)
!call RecPrt('X',' ',X,nInter,nIter)
call DGEMM_('T','N',nInter,nIter,nInter,One,layer_U,nInter,X,nInter,Zero,Y,nInter)
!call RecPrt('Y',' ',Y,nInter,nIter)

return

end subroutine Trans_K
