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

implicit none
! Output
real*8, dimension(nTri_Elem(nRoots-1)) :: bX
! Input
real*8, dimension(nRoots**2) :: R, H
real*8, dimension(nRoots) :: E_Final
real*8, dimension(nRoots**2) :: LOK
!Auxiliaries
integer I, J, K, L, M, N, IKL, IIN, IJM, IKOL, IIK, IJK, IIL, IJL, ILOK
real*8 TempD, dE_IJ

bX(:) = Zero
I = NACstates(1)
J = NACstates(2)
dE_IJ = E_Final(I)-E_Final(J)

! R_JK*HKL*RIL
do K=2,nRoots
  IIK = (I-1)*nRoots+K
  IJK = (J-1)*nRoots+K
  do L=1,K-1
    IIL = IIK-K+L
    IJL = IJK-K+L
    IKL = nTri_Elem(K-2)+L
    IKOL = (L-1)*nRoots+K
    ILOK = (K-1)*nRoots+L
    ! Diagonal elements of R_JK * H_KL * R_IK
    bX(IKL) = Two*(R(IJK)*R(IIK)*LOK(ILOK)-R(IJL)*R(IIL)*LOK(IKOL))
    ! Additional NAC term (Requires only one line)
    ! R_JK * <K|L> * R_IL =
    ! R_JL * R_IK - R_JK * R_IL
    if (.not. isMECIMSPD) bX(IKL) = bX(IKL)+dE_IJ*(R(IJL)*R(IIK)-R(IJK)*R(IIL))

    ! Off-Diagonal elements of R_JK * H_KL * R_IK
    do M=1,nRoots
      IJM = IJK-K+M
      do N=1,nRoots
        if (M == N) cycle
        TempD = Zero
        IIN = IIK-K+N
        if (M == K) TempD = TempD+H((L-1)*nRoots+N)
        if (N == K) TempD = TempD+H((M-1)*nRoots+L)
        if (M == L) TempD = TempD-H((K-1)*nRoots+N)
        if (N == L) TempD = TempD-H((M-1)*nRoots+K)
        bX(IKL) = bX(IKL)+TempD*R(IJM)*R(IIN)
      end do
    end do
  end do
end do

end subroutine CalcbX_CMSNAC
