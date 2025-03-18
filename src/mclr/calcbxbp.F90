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
      subroutine CalcbXbP(bX,bP,FMO1t,FMO2t,R,H,nTri)
      use stdalloc, only : mma_allocate, mma_deallocate
      use MCLR_Data, only: nConf1, nAcPr2
      use input_mclr, only: nRoots
      Implicit None


!***** Output
       Real*8,DIMENSION((nRoots-1)*nRoots/2)::bX
       Real*8,DIMENSION(nConf1*nRoots)::bP
!***** Input
       INTEGER nTri
       Real*8,DIMENSION(nRoots*nTri)::FMO1t
       Real*8,DIMENSION(nRoots*nacpr2)::FMO2t
       Real*8,DIMENSION(nRoots**2)::R,H
!***** Auxiliaries
       Real*8,DIMENSION(:),Allocatable::LOK,CSFOK

       CALL mma_allocate(CSFOK,nRoots*nConf1)
       CALL mma_allocate(LOK,nRoots**2)
       CALL CalcOMat(CSFOK,LOK,FMO1t,FMO2t,nTri)
       CALL CalcbP(bP,CSFOK,LOK,R)
       CALL CalcbX(bX,LOK,R,H)
       CALL mma_deallocate(CSFOK)
       CALL mma_deallocate(LOK)
       end subroutine CalcbXbP
