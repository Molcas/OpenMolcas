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

subroutine CalcHessCMS(Hess,DDg,lRoots,nSPair)

use Index_Functions, only: iTri
use Constants, only: Zero, Two, Four
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: lRoots, nSPair
real(kind=wp), intent(out) :: Hess(nSPair,nSPair)
real(kind=wp), intent(in) :: DDg(lRoots,lRoots,lRoots,lRoots)
integer(kind=iwp) :: K, L, M, N, iKL, iMN

Hess(:,:) = Zero
do K=2,lRoots
  do L=1,K-1
    iKL = iTri(K-1,L)
    do M=2,lRoots
      do N=1,M-1
        iMN = iTri(M-1,N)
        if (L == M) Hess(iMN,iKL) = Hess(iMN,iKL)+DDg(K,N,K,K)+DDg(K,N,N,N)-Two*DDg(K,N,L,L)-Four*DDg(K,L,M,N)
        if (K == N) Hess(iMN,iKL) = Hess(iMN,iKL)+DDg(L,M,L,L)+DDg(L,M,M,M)-Two*DDg(L,M,K,K)-Four*DDg(K,L,M,N)
        if (K == M) Hess(iMN,iKL) = Hess(iMN,iKL)-DDg(L,N,L,L)-DDg(L,N,N,N)+Two*DDg(L,N,K,K)+Four*DDg(K,L,M,N)
        if (L == N) Hess(iMN,iKL) = Hess(iMN,iKL)-DDg(K,M,K,K)-DDg(K,M,M,M)+Two*DDg(K,M,L,L)+Four*DDg(K,L,M,N)
      end do
    end do
  end do
end do

return

end subroutine CalcHessCMS
