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
! Based on cmsbk.f from Jie J. Bao &&                            *
! Additional work from  rhs_nac.f                                *
! ****************************************************************

subroutine Calcbk_CMSNAC(bk,R,nTri,GDMat,zX)

use stdalloc, only: mma_allocate, mma_deallocate
use MCLR_Data, only: nDens2, nNA
use input_mclr, only: nRoots, ntAsh

implicit none
! Output
real*8, dimension(nDens2) :: bk
! Input
real*8, dimension(nRoots**2) :: R
integer nTri
real*8, dimension((nRoots-1)*nRoots/2) :: zX
real*8, dimension(nRoots*(nRoots+1)/2,nnA,nnA) :: GDMat
! Auxiliaries
real*8, dimension(:), allocatable :: FOccMO, P2MOt
integer nP2, nG1
integer i, j, iTri
! Statement function
itri(i,j) = max(i,j)*(max(i,j)-1)/2+min(i,j)

ng1 = itri(ntash,ntash)
nP2 = itri(ng1,ng1)
call mma_allocate(FOccMO,nDens2)
call mma_allocate(P2MOt,nP2)
call FZero(bk,nDens2)
call GetWFFock_NAC(FOccMO,bk,R,nTri,P2MOt,nP2)
call GetQaaFock(FOccMO,P2MOt,GDMat,zX,nP2)
call GetPDFTFock_NAC(bk)
call PutCMSFockOcc(FOccMO,nTri)

call mma_deallocate(FOccMO)
call mma_deallocate(P2MOt)

end subroutine Calcbk_CMSNAC
