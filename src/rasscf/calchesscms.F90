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

subroutine CalcHessCMS(Hess,DDg,nDDg,lRoots,nSPair)

use Constants, only: Zero, Two, Four
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: nDDg, lRoots, nSPair
real(kind=wp) :: Hess(nSPair**2), DDg(nDDg)
integer(kind=iwp) :: K, L, M, N, iKL, iMN, iLoc1, iLoc2, iLoc3, iLoc4, iLoc5, lRoots2, lRoots3, lRoots23
real(kind=wp) :: Vklmn, Vklnm, Vlkmn, Vlknm

lRoots2 = lRoots**2
lRoots3 = lRoots2*lRoots
lRoots23 = lRoots2+lRoots3

do K=2,lRoots
  do L=1,K-1
    iKL = (K-2)*(K-1)/2+L
    do M=2,lRoots
      do N=1,M-1
        iMN = (M-2)*(M-1)/2+N
        Vklmn = Zero
        Vlknm = Zero
        Vlkmn = Zero
        Vklnm = Zero
        iLoc5 = K+(L-1)*lRoots+(M-1)*lRoots2+(N-1)*lRoots3
        if (L == M) then
          iLoc4 = K+(N-1)*lRoots
          iLoc1 = iLoc4+(K-1)*lRoots23
          iLoc2 = iLoc4+(N-1)*lRoots23
          iLoc3 = iLoc4+(L-1)*lRoots23
          Vklmn = DDg(iLoc1)+DDg(iLoc2)-Two*DDg(iLoc3)-Four*DDg(iLoc5)
          !Vklmn = DDg(K,N,K,K)+DDg(K,N,N,N)-Two*DDg(K,N,L,L)-Four*DDg(K,L,M,N)
        end if
        if (K == N) then
          iLoc4 = L+(M-1)*lRoots
          iLoc1 = iLoc4+(L-1)*lRoots23
          iLoc2 = iLoc4+(M-1)*lRoots23
          iLoc3 = iLoc4+(K-1)*lRoots23
          Vlknm = DDg(iLoc1)+DDg(iLoc2)-Two*DDg(iLoc3)-Four*DDg(iLoc5)
          !Vlknm = DDg(L,M,L,L)+DDg(L,M,M,M)-Two*DDg(L,M,K,K)-Four*DDg(K,L,M,N)
        end if
        if (K == M) then
          iLoc4 = L+(N-1)*lRoots
          iLoc1 = iLoc4+(L-1)*lRoots23
          iLoc2 = iLoc4+(N-1)*lRoots23
          iLoc3 = iLoc4+(K-1)*lRoots23
          Vlkmn = DDg(iLoc1)+DDg(iLoc2)-Two*DDg(iLoc3)-Four*DDg(iLoc5)
          !Vlkmn = DDg(L,N,L,L)+DDg(L,N,N,N)-Two*DDg(L,N,K,K)-Four*DDg(K,L,M,N)
        end if
        if (L == N) then
          iLoc4 = K+(M-1)*lRoots
          iLoc1 = iLoc4+(K-1)*lRoots23
          iLoc2 = iLoc4+(M-1)*lRoots23
          iLoc3 = iLoc4+(L-1)*lRoots23
          Vklnm = DDg(iLoc1)+DDg(iLoc2)-Two*DDg(iLoc3)-Four*DDg(iLoc5)
          !Vklnm = DDg(K,M,K,K)+DDg(K,M,M,M)-Two*DDg(K,M,L,L)-Four*DDg(K,L,M,N)
        end if
        Hess((iKL-1)*nSPair+iMN) = Vklmn+Vlknm-Vklnm-Vlkmn
      end do
    end do
  end do
end do

return

end subroutine CalcHessCMS
