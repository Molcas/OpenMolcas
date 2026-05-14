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

use Index_Functions, only: iTri, nTri_Elem
use input_mclr, only: nRoots
use Constants, only: Zero, Two, Four
use Definitions, only: wp, iwp

implicit none
real(kind=wp), intent(out) :: AXX(nTri_Elem(nRoots-1),nTri_Elem(nRoots-1))
real(kind=wp), intent(in) :: W(nTri_Elem(nRoots),nTri_Elem(nRoots))
integer(kind=iwp) :: IC, IKK, IKL, IKL2, ILL, IMM, IMN, IMN2, INN, K, L, M, N
real(kind=wp) :: VKLMN, VKLNM, VLKMN, VLKNM

do K=1,nRoots
  do L=1,K-1
    IKL = iTri(K,L)
    IKK = nTri_Elem(K)
    ILL = nTri_Elem(L)
    IKL2 = nTri_Elem(K-2)+L
    do M=1,nRoots
      do N=1,M-1
        IMN = iTri(M,N)
        IMM = nTri_Elem(M)
        INN = nTri_Elem(N)
        IMN2 = nTri_Elem(M-2)+N
        VKLMN = Zero
        VLKNM = Zero
        VLKMN = Zero
        VKLNM = Zero
        if (L == M) then
          IC = iTri(K,N)
          VKLMN = W(IC,IKK)+W(IC,INN)-Two*W(IC,ILL)-Four*W(IKL,IMN)
        end if
        if (K == N) then
          IC = iTri(L,M)
          VLKNM = W(IC,ILL)+W(IC,IMM)-Two*W(IC,IKK)-Four*W(IKL,IMN)
        end if
        if (K == M) then
          IC = iTri(L,N)
          VLKMN = W(IC,ILL)+W(IC,INN)-Two*W(IC,IKK)-Four*W(IKL,IMN)
        end if
        if (L == N) then
          IC = iTri(K,M)
          VKLNM = W(IC,IKK)+W(IC,IMM)-Two*W(IC,ILL)-Four*W(IKL,IMN)
        end if
        AXX(IMN2,IKL2) = VKLMN+VLKNM-VKLNM-VLKMN
      end do
    end do
  end do
end do

end subroutine CalcAXX
