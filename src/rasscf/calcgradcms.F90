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
!*****************************************************************
! history:                                                       *
! Jie J. Bao, on Apr. 12, 2022, created this file.               *
!*****************************************************************

subroutine CalcGradCMS(Grad,DDg,nDDg,lRoots,nSPair)

use Constants, only: Two

implicit none
integer nDDg, lRoots, nSPair
real*8 Grad(nSPair), DDg(nDDg)
integer K, L, iKL, lRoots2, lRoots3, iLoc1, iLoc2

lRoots2 = lRoots**2
lRoots3 = lRoots*lRoots2

do K=2,lRoots
  do L=1,K-1
    iLoc1 = K+(K-1)*lRoots+(K-1)*lRoots2+(L-1)*lRoots3
    iLoc2 = L+(L-1)*lRoots+(K-1)*lRoots2+(L-1)*lRoots3
    iKL = (K-1)*(K-2)/2+L
    Grad(iKL) = DDg(iLoc1)-DDg(iLoc2)
  end do
end do
call DSCal_(nSPair,Two,Grad,1)

return

end subroutine CalcGradCMS
