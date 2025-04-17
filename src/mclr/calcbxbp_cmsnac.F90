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

use Index_Functions, only: nTri_Elem
use MCLR_Data, only: nAcPr2, nConf1
use input_mclr, only: nRoots
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp

implicit none
real(kind=wp), intent(out) :: bX(nTri_Elem(nRoots-1)), bP(nConf1,nRoots)
integer(kind=iwp), intent(in) :: nTri
real(kind=wp), intent(in) :: FMO1t(nTri,nRoots), FMO2t(nAcPr2,nRoots), R(nRoots,nRoots), H(nRoots,nRoots), E_Final(nRoots)
real(kind=wp), allocatable :: CSFOK(:,:), LOK(:,:)

call mma_allocate(CSFOK,nConf1,nRoots)
call mma_allocate(LOK,nRoots,nRoots)
! Using CalcOMat in original CalcbXbP
call CalcOMat(CSFOK,LOK,FMO1t,FMO2t,nTri)
call CalcbP_CMSNAC(bP,CSFOK,LOK,R)
call CalcbX_CMSNAC(bX,LOK,R,H,E_Final)
call mma_deallocate(CSFOK)
call mma_deallocate(LOK)

end subroutine CalcbXbP_CMSNAC
