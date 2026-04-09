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

use MCLR_Data, only: nConf1, nDens
use input_mclr, only: nRoots
use Definitions, only: wp

implicit none
real(kind=wp), intent(out) :: Fock(nDens), CICSF(nConf1*nRoots)
real(kind=wp), intent(in) :: AXkzx(nDens), AXPzx(nConf1*nRoots), bk(nDens), bP(nConf1*nRoots)

! Orbital Rotation Part
Fock(:) = Axkzx(:)+bk(:)

! State-CSF Rotation Part
CICSF(:) = AXPzx(:)-bP(:)

end subroutine SolveforRHS
