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

subroutine CalcbX(bX,LOK,R,H)

use Index_Functions, only: nTri_Elem
use MCLR_Data, only: IRLXROOT
use input_mclr, only: nRoots
use Constants, only: Zero, Two
use Definitions, only: wp, iwp

implicit none
real(kind=wp), intent(out) :: bX(nTri_Elem(nRoots-1))
real(kind=wp), intent(in) :: LOK(nRoots,nRoots), R(nRoots,nRoots), H(nRoots,nRoots)
integer(kind=iwp) :: I, IKL, K, L, M, N
real(kind=wp) :: TempD

bX(:) = Zero
I = irlxroot
do K=2,nRoots
  do L=1,K-1
    IKL = nTri_Elem(K-2)+L
    bX(IKL) = R(K,I)**2*LOK(L,K)-R(L,I)**2*LOK(K,L)
    do M=2,nRoots
      do N=1,M-1
        TempD = Zero
        if (M == K) TempD = TempD+H(N,L)
        if (N == K) TempD = TempD+H(L,M)
        if (M == L) TempD = TempD-H(N,K)
        if (N == L) TempD = TempD-H(K,M)
        bX(IKL) = bX(IKL)+TempD*R(M,I)*R(N,I)
      end do
    end do
  end do
end do
bX(:) = Two*bX(:)

end subroutine CalcbX
