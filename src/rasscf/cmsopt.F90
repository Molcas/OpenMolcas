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
use rasscf_global, only: CMSStartMat, lRoots, NAC, NACPR2
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp

implicit none
real(kind=wp), intent(in) :: TUVX(NACPR2)
integer(kind=iwp) :: lRoots2, NAC2, nGD, nTUVX
real(kind=wp), allocatable :: Dgorbit(:), Dgstate(:), GDorbit(:), GDstate(:), Gtuvx(:), RotMat(:,:)
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
call mma_allocate(GDstate,nGD)
call mma_allocate(Dgstate,nGD)
call mma_allocate(GDorbit,nGD)
call mma_allocate(Dgorbit,nGD)
call mma_allocate(Gtuvx,nTUVX)
call mma_allocate(RGD,lRoots2)
call mma_allocate(RotMat,lRoots,lRoots)

! Calculate generalized density mtrix
call LoadGtuvx(TUVX,Gtuvx)

!write(u6,*) 'Gtuvx matrix'
!call RecPrt(' ',' ',Gtuvx,NAC2,NAC2)

call CalcGD(GDorbit,nGD)
!write(u6,*) 'GD matrix orbital-leading'
!call RecPrt(' ',' ',GDorbit,NAC2,lRoots2)
call DGEMM_('T','N',NAC**2,lRoots**2,NAC**2,One,Gtuvx,NAC**2,GDorbit,NAC**2,Zero,Dgorbit,NAC**2)
!write(u6,*) 'Dg matrix orbital-leading'
!call RecPrt(' ',' ',Dgorbit,NAC2,lRoots2)

call mma_deallocate(Gtuvx)

call Trnsps(NAC2,lRoots2,Dgorbit,Dgstate)
call Trnsps(NAC2,lRoots2,GDorbit,GDstate)

! Load initial rotation matrix
call InitRotMat(RotMat,lRoots,trim(CMSStartMat),len_trim(CMSStartMat))

! Print header of CMS iterations
call CMSHeader(trim(CMSStartMat))

! Start CMS Optimization
CMSNotConverged = .true.
call CMSNewton(RotMat,GDorbit,GDstate,Dgorbit,Dgstate,nGD)

! Print end of CMS intermediate-state optimization
call CMSTail()

! Save rotation matrix
call PrintMat('ROT_VEC','CMS-PDFT',RotMat,lroots,lroots,7,8,'T')

! releasing memory
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
