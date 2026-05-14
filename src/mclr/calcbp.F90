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

subroutine CalcbP(bP,CSFOK,LOK,R)

use ipPage, only: W
use MCLR_Data, only: ipCI, IRLXROOT, nConf1
use input_mclr, only: nRoots
use Constants, only: Two
use Definitions, only: wp, iwp

implicit none
real(kind=wp), intent(out) :: bP(nConf1,nRoots)
real(kind=wp), intent(in) :: CSFOK(nConf1,nRoots), LOK(nRoots,nRoots), R(nRoots,nRoots)
integer(kind=iwp) :: I, K, L

bP(:,:) = CSFOK(:,:)
I = irlxroot
do K=1,nRoots
  do L=1,nRoots
    bP(:,K) = bP(:,K)-LOK(L,K)*W(ipci)%A((L-1)*nConf1+1:L*nConf1)
  end do

  bP(:,K) = Two*R(K,I)**2*bP(:,K)

end do

end subroutine CalcbP
