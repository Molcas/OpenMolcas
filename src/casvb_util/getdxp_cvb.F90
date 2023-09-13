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

implicit real*8(a-h,o-z)
dimension dxp(npr), gradp(npr), heigval(npr)

do i=1,nnegeig
  dxp(i) = -gradp(i)/(heigval(i)-alfa)
end do
do i=nnegeig+1,npr
  dxp(i) = -gradp(i)/(heigval(i)+alfa)
end do

return

end subroutine getdxp_cvb
