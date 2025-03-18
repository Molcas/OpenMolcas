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
      subroutine CalcbP_CMSNAC(bP,CSFOK,LOK,R)
      use ipPage, only: W
      use MCLR_Data, only: nConf1, ipCI
      use MCLR_Data, only: NACSTATES
      use input_mclr, only: nRoots
      implicit none
!**** Output
      Real*8,DIMENSION(nConf1*nRoots)::bP
!**** Input
      Real*8,DIMENSION(nRoots*nConf1)::CSFOK
      Real*8,DIMENSION(nRoots**2)::LOK
      Real*8,DIMENSION(nRoots**2)::R
!**** Kind quantities that help
      INTEGER I,J,L,K,iLoc1,iLoc2
      Real*8 tempd

      I=NACstates(1)
      J=NACstates(2)
      DO K=1,nRoots
       iLoc1=(K-1)*nConf1+1
       CALL DCopy_(nConf1,CSFOK(iLoc1),1,bP(iLoc1),1)
       Do L=1, nRoots
        tempd=-LOK((K-1)*nRoots+L)
        iLoc2=(L-1)*nConf1+1
        CALL dAXpY_(nConf1,tempd,W(ipci)%Vec(iLoc2),1,bP(iLoc1),1)
       End Do

       CALL DScal_(nConf1,2*R((J-1)*nRoots+K)*R((I-1)*nRoots+K),        &
     & bP(iLoc1),1)

      END DO
      End Subroutine CalcbP_CMSNAC
