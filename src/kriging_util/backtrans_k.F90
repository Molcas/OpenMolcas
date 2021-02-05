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

subroutine BackTrans_K(X,Y,nInter,nIter)

use kriging_mod, only: layer_U

implicit none
integer nInter, nIter
real*8 X(nInter,nIter), Y(nInter,nIter)

!call RecPrt('layer_U',' ',layer_U,nInter,nInter)
!call RecPrt('X',' ',X,nInter,nIter)
call DGEMM_('N','N',nInter,nIter,nInter,1.0d0,layer_U,nInter,X,nInter,0.0d0,Y,nInter)
!call RecPrt('Y',' ',Y,nInter,nIter)

end subroutine BackTrans_K

!-------------------------------------------------------------------------

subroutine BackTrans_Kt(X,Y,nInter,nIter)

use kriging_mod, only: layer_U

implicit none
integer nInter, nIter
real*8 X(nInter,nIter), Y(nInter,nIter)

!call RecPrt('layer_U',' ',layer_U,nInter,nInter)
!call RecPrt('X',' ',X,nInter,nIter)
call DGEMM_('N','T',nInter,nIter,nInter,1.0d0,X,nInter,layer_U,nInter,0.0d0,Y,nInter)
!call RecPrt('Y',' ',Y,nInter,nIter)

end subroutine BackTrans_Kt
