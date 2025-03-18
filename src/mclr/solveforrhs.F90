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

subroutine SolveforRHS(Fock,CICSF,AXkzx,AXPzx,bk,bP)

use MCLR_Data, only: nDens2, nConf1
use input_mclr, only: nRoots

implicit none
! Output
real*8, dimension(nDens2+6) :: Fock
real*8, dimension(nconf1*nroots) :: CICSF
! Input
real*8, dimension(nDens2) :: AXkzx
real*8, dimension(NConf1*nRoots) :: AXPzx
real*8, dimension(nDens2) :: bk
real*8, dimension(nConf1*nRoots) :: bP
! Assistants
integer nRow

! Orbital Rotation Part
nRow = nDens2
call FZero(Fock,nDens2)
call DCopy_(nRow,Axkzx,1,Fock,1)
call DAXPY_(nRow,1.0d0,bk,1,Fock,1)

! State-CSF Rotation Part
nRow = nRoots*nConf1
call FZero(CICSF,nRow)
call DCopy_(nRow,AXPzx,1,CICSF,1)
call DAXPY_(nRow,-1.0d0,bP,1,CICSF,1)

end subroutine SolveforRHS
