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
! Copyright (C) 2022, Jie J. Bao                                       *
!***********************************************************************

subroutine CMSOpt(TUVX)
! ****************************************************************
! history:                                                       *
! Jie J. Bao, on Apr. 07, 2022, created this file.               *
! ****************************************************************

use CMS, only: CMSNotConverged, RGD
use rasscf_global, only: NACPR2, CMSStartMat, lRoots, NAC
use stdalloc, only: mma_allocate, mma_deallocate

implicit none
real*8, dimension(NACPR2) :: TUVX
real*8, dimension(:), allocatable :: Gtuvx, R, GDstate, GDorbit, Dgstate, Dgorbit
real*8, dimension(:,:), allocatable :: RotMat
integer nTUVX, nGD, lRoots2, NAC2
character(len=16) :: VecName
#include "warnings.h"

!*****************************************************************
!     some notes on the arrays:                                  *
!     Gtuvx : two-electron integral, g_tuvx                      *
!     GD    : "generalized 1-e density matrix"                   *
!              GD^KL: transition density matrix from L to K      *
!              GD^KK: density matrix for state K                 *
!     Dg    :  sum_{vx}{GD^KL_vx * g_tuvx}                       *
!     In GDorbit and Dgorbit, the leading index is orbital index;*
!     In GDstate and Dgstate, the leading index is state index.  *
!                                                                *
!     DDg   :  sum_{tuvx}{GD^KL_tu * GD^MN_vx * g_tuvx}          *
!      namely, sum_{tu}{GD^KL_tu * Dg^MN_tu}                     *
!*****************************************************************

NAC2 = NAC**2
nTUVX = NAC2**2
lRoots2 = lRoots**2
nGD = lRoots2*NAC2

CMSNotConverged = .true.

! Memory Allocation
call mma_allocate(R,lRoots2)
call mma_allocate(GDstate,nGD)
call mma_allocate(Dgstate,nGD)
call mma_allocate(GDorbit,nGD)
call mma_allocate(Dgorbit,nGD)
call mma_allocate(Gtuvx,nTUVX)
call mma_allocate(RGD,lRoots2)
call mma_allocate(RotMat,lRoots,lRoots)

! Calculate generalized density mtrix
call UnzipTUVX(TUVX,Gtuvx,nTUVX)

!write(6,*) 'Gtuvx matrix'
!call RecPrt(' ',' ',Gtuvx,NAC2,NAC2)

call CalcGD(GDorbit,nGD)
!write(6,*) 'GD matrix orbital-leading'
!call RecPrt(' ',' ',GDorbit,NAC2,lRoots2)
call CalcDg(Dgorbit,GDorbit,Gtuvx,nGD,nTUVX,NAC,lRoots)
!write(6,*) 'Dg matrix orbital-leading'
!call RecPrt(' ',' ',Dgorbit,NAC2,lRoots2)

call mma_deallocate(Gtuvx)

call TransposeMat(Dgstate,Dgorbit,nGD,NAC2,lRoots2)
call TransposeMat(GDstate,GDorbit,nGD,NAC2,lRoots2)

! Load initial rotation matrix
call InitRotMat(RotMat,lRoots,trim(CMSStartMat),len_trim(CMSStartMat))

call OneDFoil(R,RotMat,lRoots,lRoots)

! Print header of CMS iterations
call CMSHeader(trim(CMSStartMat),len_trim(CMSStartMat))

! Start CMS Optimization
CMSNotConverged = .true.
call CMSNewton(R,GDorbit,GDstate,Dgorbit,Dgstate,nGD)

! Print end of CMS intermediate-state optimization
call CMSTail()

! Save rotation matrix
call AntiOneDFoil(RotMat,R,lRoots,lRoots)
VecName = 'CMS-PDFT'
call PrintMat('ROT_VEC',VecName,RotMat,lroots,lroots,7,16,'T')

! releasing memory
call mma_deallocate(R)
call mma_deallocate(GDstate)
call mma_deallocate(Dgstate)
call mma_deallocate(GDorbit)
call mma_deallocate(Dgorbit)
call mma_deallocate(RGD)
call mma_deallocate(RotMat)

! check convergence
if (CMSNotConverged) then
  call WarningMessage(2,'CMS Intermediate States Not Converged')
  call Quit(_RC_NOT_CONVERGED_)
end if

end subroutine CMSOpt
