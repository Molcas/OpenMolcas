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

use Index_Functions, only: iTri, nTri_Elem
use Constants, only: Zero, Four
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: M, nnA, nRoots
real(kind=wp), intent(out) :: Dacc(nnA,nnA)
real(kind=wp), intent(in) :: GDMat(nTri_Elem(nRoots),nnA,nnA), zx(nTri_Elem(nRoots-1))
integer(kind=iwp) :: IKM, IKM2, it, K
real(kind=wp) :: Fact

Dacc(:,:) = Zero

do K=1,nRoots
  if (K == M) cycle
  IKM = iTri(K,M)
  IKM2 = nTri_Elem(max(K,M)-2)+min(K,M)
  Fact = Four*zx(IKM2)
  if (K > M) Fact = -Fact
  do it=1,nnA
    Dacc(:,it) = Dacc(:,it)+Fact*GDMat(IKM,it,:)
  end do
end do

end subroutine CalcDacc
