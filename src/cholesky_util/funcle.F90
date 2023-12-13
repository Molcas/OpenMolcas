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
! Copyright (C) 2007, Ten-no Research Group                            *
!               2012, Thomas Bondo Pedersen                            *
!***********************************************************************

function FuncLe(X,Delta)
!----------------------------------------------------------------------
! Function : Define function after CoV for Gauss-Legendre Method
!            H. Koch, A. Sanchez de Meras (J. Chem. Phys. 113, 508)
!-----------------------------------------------------------------------

use Constants, only: One
use Definitions, only: wp

implicit none
real(kind=wp) :: FuncLe
real(kind=wp), intent(in) :: X, Delta
real(kind=wp) :: DerV, Varl

Varl = X/(One-X)
Derv = One/((One-X)*(One-X))
FuncLe = Derv*exp(-Delta*Varl)

return

end function FuncLe
