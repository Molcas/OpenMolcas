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
! Copyright (C) 2003, Per-Olof Widmark                                 *
!***********************************************************************

real*8 function optim_E(C,G,H,n,nDim)

use Constants, only: Zero, Half

implicit none
integer n, nDim
real*8 C(nDim), G(nDim), H(nDim,nDim)
integer k, m
real*8 Tmp

Optim_E = Zero
do k=1,n
  Tmp = Zero
  do m=1,n
    Tmp = Tmp+Half*(C(k)*C(m)*H(k,m))
  end do
  Optim_E = Optim_E+C(k)*G(k)+Tmp
end do

end function optim_E
