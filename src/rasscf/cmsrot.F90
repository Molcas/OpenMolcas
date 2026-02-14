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
! Copyright (C) 2020, Jie J. Bao                                       *
!***********************************************************************

subroutine CMSRot(TUVX)
! ****************************************************************
! history:                                                       *
! Jie J. Bao, on Aug. 06, 2020, created this file.               *
! ****************************************************************

use CMS, only: CMSNotConverged
use rasscf_global, only: CMSStartMat, lRoots, NAC, NACPR2
use PrintLevel, only: USUAL
use output_ras, only: IPRLOC
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp, u6

implicit none
real(kind=wp) :: TUVX(NACPR2)
integer(kind=iwp) :: iPrLev
character(len=16) :: VecName
real(kind=wp), allocatable :: DDG(:,:,:,:), GDMat(:,:,:), Gtuvx(:,:,:,:), RotMat(:,:)
#include "warnings.h"

! Allocating Memory
call mma_allocate(GDMat,LRoots*(LRoots+1)/2,NAC,NAC)
call mma_allocate(RotMat,lRoots,lRoots)
call mma_allocate(Gtuvx,NAC,NAC,NAC,NAC)
call mma_allocate(DDG,lRoots,lRoots,lRoots,lRoots)

IPRLEV = IPRLOC(6)

! printing header
if (IPRLEV >= USUAL) then
  write(u6,*)
  write(u6,*)
  write(u6,*) '    CMS INTERMEDIATE-STATE OPTIMIZATION'
end if
if (trim(CMSStartMat) == 'XMS') then
  call ReadMat('ROT_VEC',VecName,RotMat,lroots,lroots,7,16,'N')
else
  call ReadMat(trim(CMSStartMat),VecName,RotMat,lroots,lroots,len_trim(CMSStartMat),16,'N')
end if
if (IPRLEV >= USUAL) call CMSHeader(trim(CMSStartMat),len_trim(CMSStartMat))

call LoadGtuvx(TUVX,Gtuvx)

CMSNotConverged = .false.
call GetGDMat(GDMat)
if (lRoots < NAC) then
  !write(u6,*) 'Optimization Approach 1'
  call GetDDgMat(DDg,GDMat,Gtuvx)
  call NStateOpt(RotMat,DDg)
else
  !write(u6,*) 'Optimization Approach 2'
  call NStateOpt2(RotMat,GDMat,Gtuvx)
end if
call PrintMat('ROT_VEC','CMS-PDFT',RotMat,lroots,lroots,7,8,'N')

! Deallocating Memory
call mma_deallocate(GDMat)
call mma_deallocate(RotMat)
call mma_deallocate(Gtuvx)
call mma_deallocate(DDg)
if (CMSNotConverged) then
  call WarningMessage(2,'CMS Intermediate States Not Converged')
  call Quit(_RC_NOT_CONVERGED_)
end if

end subroutine CMSRot
