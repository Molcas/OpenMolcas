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
      Subroutine SolveforRHS(Fock,CICSF,AXkzx,AXPzx,bk,bP)
      use MCLR_Data, only: nDens2, nConf1
      use input_mclr, only: nRoots
      Implicit None
!***** Output
      Real*8,DIMENSION(nDens2+6)::Fock
      Real*8,DIMENSION(nconf1*nroots)::CICSF
!***** Input
      Real*8,DIMENSION(nDens2)::AXkzx
      Real*8,DIMENSION(NConf1*nRoots)::AXPzx
      Real*8,DIMENSION(nDens2)::bk
      Real*8,DIMENSION(nConf1*nRoots)::bP
!***** Assistants
      INTEGER nRow

!****  Orbital Rotation Part
      nRow=nDens2
      CALL FZero(Fock,nDens2)
      CALL DCopy_(nRow,Axkzx,1,Fock,1)
      CALL DAXPY_(nRow,1.0d0,bk,1,Fock,1)

!***** State-CSF Rotation Part
      nRow=nRoots*nConf1
      CALL FZero(CICSF,nRow)
      CALL DCopy_(nRow,AXPzx,1,CICSF,1)
      CALL DAXPY_(nRow,-1.0d0,bP,1,CICSF,1)

      END SUBROUTINE SolveforRHS
