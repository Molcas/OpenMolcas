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
use MCLR_Data, only: nConf1, ipCI
use MCLR_Data, only: IRLXROOT
use input_mclr, only: nRoots
use Constants, only: Two

implicit none
! Output
real*8, dimension(nConf1*nRoots) :: bP
! Input
real*8, dimension(nRoots*nConf1) :: CSFOK
real*8, dimension(nRoots**2) :: LOK
real*8, dimension(nRoots**2) :: R
! Kind quantities that help
integer I, L, K, iLoc1, iLoc2
real*8 tempd

I = irlxroot
do K=1,nRoots
  iLoc1 = (K-1)*nConf1+1
  call DCopy_(nConf1,CSFOK(iLoc1),1,bP(iLoc1),1)
  do L=K,K
    tempd = -LOK((K-1)*nRoots+L)
    iLoc2 = (L-1)*nConf1+1
    call dAXpY_(nConf1,tempd,W(ipci)%A(iLoc2),1,bP(iLoc1),1)
  end do

  do L=1,nRoots
    if (L == K) cycle
    tempd = -LOK((K-1)*nRoots+L)
    iLoc2 = (L-1)*nConf1+1
    call dAXpY_(nConf1,tempd,W(ipci)%A(iLoc2),1,bP(iLoc1),1)
  end do
  call DScal_(nConf1,Two*R((I-1)*nRoots+K)**2,bP(iLoc1),1)
end do

end subroutine CalcbP
