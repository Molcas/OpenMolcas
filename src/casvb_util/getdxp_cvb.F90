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
! Copyright (C) 1996-2006, Thorstein Thorsteinsson                     *
!               1996-2006, David L. Cooper                             *
!***********************************************************************

subroutine getdxp_cvb(dxp,gradp,heigval,nnegeig,npr,alfa)

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nnegeig, npr
real(kind=wp), intent(out) :: dxp(npr)
real(kind=wp), intent(in) :: gradp(npr), heigval(npr), alfa

dxp(1:nnegeig) = -gradp(1:nnegeig)/(heigval(1:nnegeig)-alfa)
dxp(nnegeig+1:) = -gradp(nnegeig+1:)/(heigval(nnegeig+1:)+alfa)

return

end subroutine getdxp_cvb
