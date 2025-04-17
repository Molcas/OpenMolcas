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
! Copyright (C) 2021, Paul B Calio                                     *
!***********************************************************************
! ****************************************************************
! history:                                                       *
! Based on cmsbxbp.f from Jie J. Bao                             *
! ****************************************************************

subroutine CalcbX_CMSNAC(bX,LOK,R,H,E_Final)

use Index_Functions, only: nTri_Elem
use MCLR_Data, only: ISMECIMSPD, NACSTATES
use input_mclr, only: nRoots
use Constants, only: Zero, Two
use Definitions, only: wp, iwp

implicit none
real(kind=wp), intent(out) :: bX(nTri_Elem(nRoots-1))
real(kind=wp), intent(in) :: LOK(nRoots,nRoots), R(nRoots,nRoots), H(nRoots,nRoots), E_Final(nRoots)
integer(kind=iwp) :: I, IKL, J, K, L, M, N
real(kind=wp) :: dE_IJ, TempD

bX(:) = Zero
I = NACstates(1)
J = NACstates(2)
dE_IJ = E_Final(I)-E_Final(J)

! R_JK*H_KL*R_IL
do K=2,nRoots
  do L=1,K-1
    IKL = nTri_Elem(K-2)+L
    ! Diagonal elements of R_JK * H_KL * R_IK
    bX(IKL) = Two*(R(K,J)*R(K,I)*LOK(L,K)-R(L,J)*R(L,I)*LOK(K,L))
    ! Additional NAC term (Requires only one line)
    ! R_JK * <K|L> * R_IL =
    ! R_JL * R_IK - R_JK * R_IL
    if (.not. isMECIMSPD) bX(IKL) = bX(IKL)+dE_IJ*(R(L,J)*R(K,I)-R(K,J)*R(L,I))

    ! Off-Diagonal elements of R_JK * H_KL * R_IK
    do M=1,nRoots
      do N=1,nRoots
        if (M == N) cycle
        TempD = Zero
        if (M == K) TempD = TempD+H(N,L)
        if (N == K) TempD = TempD+H(L,M)
        if (M == L) TempD = TempD-H(N,K)
        if (N == L) TempD = TempD-H(K,M)
        bX(IKL) = bX(IKL)+TempD*R(M,J)*R(N,I)
      end do
    end do
  end do
end do

end subroutine CalcbX_CMSNAC
