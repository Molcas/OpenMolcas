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

subroutine Put_Grad(Grad,nGrad)

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: nGrad
real(kind=wp) :: Grad(nGrad)
integer(kind=iwp) :: iGO
character(len=24) :: Label

Label = 'GRAD'
call Put_dArray(Label,Grad,nGrad)

call Get_iScalar('Grad ready',iGO)
iGO = ior(iGO,2**0)
call Put_iScalar('Grad ready',iGO)

return

end subroutine Put_Grad
