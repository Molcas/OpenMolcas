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
! Copyright (C) 1992, Roland Lindh                                     *
!***********************************************************************

subroutine Fold2(nSym,nBas,A,B)

implicit none
integer nSym, nBas(*)
real*8 A(*), B(*)
integer iOff1, iOff2, iSym, mBas, iBas, jBas

iOff1 = 0
iOff2 = 0
do iSym=1,nSym
  mBas = nBas(iSym)
  do iBas=1,mBas
    do jBas=1,iBas-1
      B(iOff2+jBas) = A(iOff1+jBas)
    end do
    B(iOff2+iBas) = A(iOff1+iBas)
    iOff1 = iOff1+mBas
    iOff2 = iOff2+iBas
  end do
end do

return

end subroutine Fold2
