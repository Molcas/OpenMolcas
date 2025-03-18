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

subroutine CalcbXbP_CMSNAC(bX,bP,FMO1t,FMO2t,R,H,E_Final,nTri)

use stdalloc, only: mma_allocate, mma_deallocate
use MCLR_Data, only: nCOnf1, nAcPr2
use input_mclr, only: nRoots

implicit none
! Output
real*8, dimension((nRoots-1)*nRoots/2) :: bX
real*8, dimension(nConf1*nRoots) :: bP
! Input
integer nTri
real*8, dimension(nRoots*nTri) :: FMO1t
real*8, dimension(nRoots*nacpr2) :: FMO2t
real*8, dimension(nRoots**2) :: R, H
real*8, dimension(nRoots) :: E_Final
! Auxiliaries
real*8, dimension(:), allocatable :: LOK, CSFOK

call mma_allocate(CSFOK,nRoots*nConf1)
call mma_allocate(LOK,nRoots**2)
! Using CalcOMat in original CalcbXbP
call CalcOMat(CSFOK,LOK,FMO1t,FMO2t,nTri)
call CalcbP_CMSNAC(bP,CSFOK,LOK,R)
call CalcbX_CMSNAC(bX,LOK,R,H,E_Final)
call mma_deallocate(CSFOK)
call mma_deallocate(LOK)

end subroutine CalcbXbP_CMSNAC
