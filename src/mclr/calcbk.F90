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

subroutine Calcbk(bk,R,nTri,GDMat,zX)

use Index_Functions, only: nTri_Elem
use MCLR_Data, only: nDens, nNA
use input_mclr, only: nRoots, ntAsh
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
real(kind=wp), intent(out) :: bk(nDens)
real(kind=wp), intent(in) :: R(nRoots,nRoots), GDMat(nTri_Elem(nRoots),nnA,nnA), zX(nTri_Elem(nRoots-1))
integer(kind=iwp), intent(in) :: nTri
integer(kind=iwp) :: nG1, nP2
real(kind=wp), allocatable :: FOccMO(:), P2MOt(:)

ng1 = nTri_Elem(ntash)
nP2 = nTri_Elem(ng1)
call mma_allocate(FOccMO,nDens)
call mma_allocate(P2MOt,nP2)
bk(:) = Zero
call GetWFFock(FOccMO,bk,R,nTri,P2MOt,nP2)
call GetQaaFock(FOccMO,P2MOt,GDMat,zX,nP2)
call GetPDFTFock(bk)
call PutCMSFockOcc(FOccMO,nTri)

call mma_deallocate(FOccMO)
call mma_deallocate(P2MOt)

end subroutine Calcbk
