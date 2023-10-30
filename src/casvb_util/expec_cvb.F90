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

subroutine expec_cvb(dxp,gradp,heigval,nnegeig,npr,expc,exp1,exp2)

use Constants, only: Zero, Half
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nnegeig, npr
real(kind=wp), intent(in) :: dxp(npr), gradp(npr), heigval(npr)
real(kind=wp), intent(out) :: expc, exp1, exp2
integer(kind=iwp) :: i
real(kind=wp) :: exp1l, exp2l

exp1l = Zero
do i=1,nnegeig
  exp1l = exp1l+dxp(i)*(gradp(i)+Half*dxp(i)*heigval(i))
end do
exp2l = Zero
do i=nnegeig+1,npr
  exp2l = exp2l+dxp(i)*(gradp(i)+Half*dxp(i)*heigval(i))
end do
expc = exp1l+exp2l
exp1 = exp1l
exp2 = exp2l

return

end subroutine expec_cvb
