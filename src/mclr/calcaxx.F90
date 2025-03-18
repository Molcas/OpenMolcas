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

subroutine CalcAXX(AXX,W)

use Constants, only: Zero, Two, Four
use input_mclr, only: nRoots

implicit none
! Input
real*8, dimension((nRoots+1)*nRoots/2,(nRoots+1)*nRoots/2) :: W
! Output
real*8, dimension(((nRoots-1)*nRoots/2)**2) :: AXX
! Auxiliary Quantities
integer K, L, M, N, IKL, IMN, IKL2, IMN2, IKK, ILL, IMM, INN, IC, nRTri
real*8 VKLMN, VLKNM, VKLNM, VLKMN

nRTri = (nRoots-1)*nRoots/2
do K=1,nRoots
  do L=1,K-1
    IKL = (K-1)*K/2+L
    IKK = (K+1)*K/2
    ILL = (L+1)*L/2
    IKL2 = (K-2)*(K-1)/2+L
    do M=1,nRoots
      do N=1,M-1
        IMN = (M-1)*M/2+N
        IMM = (M+1)*M/2
        INN = (N+1)*N/2
        IMN2 = (M-2)*(M-1)/2+N
        VKLMN = Zero
        VLKNM = Zero
        VLKMN = Zero
        VKLNM = Zero
        if (L == M) then
          if (N < K) then
            IC = (K-1)*K/2+N
          else
            IC = (N-1)*N/2+K
          end if
          VKLMN = W(IC,IKK)+W(IC,INN)-Two*W(IC,ILL)-Four*W(IKL,IMN)
        end if
        if (K == N) then
          if (M < L) then
            IC = (L-1)*L/2+M
          else
            IC = (M-1)*M/2+L
          end if
          VLKNM = W(IC,ILL)+W(IC,IMM)-Two*W(IC,IKK)-Four*W(IKL,IMN)
        end if
        if (K == M) then
          if (N < L) then
            IC = (L-1)*L/2+N
          else
            IC = (N-1)*N/2+L
          end if
          VLKMN = W(IC,ILL)+W(IC,INN)-Two*W(IC,IKK)-Four*W(IKL,IMN)
        end if
        if (L == N) then
          if (M < K) then
            IC = (K-1)*K/2+M
          else
            IC = (M-1)*M/2+K
          end if
          VKLNM = W(IC,IKK)+W(IC,IMM)-Two*W(IC,ILL)-Four*W(IKL,IMN)
        end if
        AXX((IKL2-1)*nRTri+IMN2) = VKLMN+VLKNM-VKLNM-VLKMN
      end do
    end do
  end do
end do

end subroutine CalcAXX
