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
! Copyright (C) 2021, Jie J. Bao                                       *
!***********************************************************************
! ****************************************************************
! history:                                                       *
! Jie J. Bao, on Aug. 06, 2020, created this file.               *
! ****************************************************************

subroutine CalcDacc(Dacc,GDMat,M,nnA,nRoots,zx)

use Constants, only: Zero

implicit none
integer nnA, nRoots, M
real*8, dimension((nRoots+1)*nRoots/2,nnA,nnA) :: GDMat
real*8, dimension(nnA**2) :: Dacc
real*8, dimension((nRoots-1)*nRoots/2) :: zx
integer it, iu, K, IKM, IKM2, iLoc1, iLoc2
real*8 Fact

Dacc(:) = Zero

do K=1,nRoots
  if (K == M) cycle
  if (M > K) then
    IKM = (M-1)*M/2+K
    IKM2 = (M-1)*(M-2)/2+K
  else
    IKM = (K-1)*K/2+M
    iKM2 = (K-1)*(K-2)/2+M
  end if
  Fact = 4.0d0*zx(IKM2)
  if (K > M) Fact = -Fact
  do it=1,nnA
    iLoc1 = (it-1)*nnA
    do iu=1,nnA
      iLoc2 = iLoc1+iu
      Dacc(iLoc2) = Dacc(iLoc2)+GDMat(IKM,it,iu)*Fact
    end do
  end do
end do

end subroutine CalcDacc
