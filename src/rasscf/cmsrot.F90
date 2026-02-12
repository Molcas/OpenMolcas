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
use rasscf_global, only: NACPR2, CMSStartMat, lRoots, NAC
use PrintLevel, only: USUAL
use output_ras, only: LF, IPRLOC
use stdalloc, only: mma_allocate, mma_deallocate

implicit none
real*8, dimension(NACPR2) :: TUVX
character(len=16) :: VecName
real*8, dimension(:,:,:,:), allocatable :: Gtuvx
real*8, dimension(:,:,:,:), allocatable :: DDG
real*8, dimension(:,:,:), allocatable :: GDMat
real*8, dimension(:,:), allocatable :: RotMat
integer iPrLev
#include "warnings.h"

! Allocating Memory
call mma_allocate(GDMat,LRoots*(LRoots+1)/2,NAC,NAC)
call mma_allocate(RotMat,lRoots,lRoots)
call mma_allocate(Gtuvx,NAC,NAC,NAC,NAC)
call mma_allocate(DDG,lRoots,lRoots,lRoots,lRoots)

IPRLEV = IPRLOC(6)

! printing header
if (IPRLEV >= USUAL) then
  write(LF,*)
  write(LF,*)
  write(LF,*) '    CMS INTERMEDIATE-STATE OPTIMIZATION'
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
  !write(6,*) 'Optimization Approach 1'
  call GetDDgMat(DDg,GDMat,Gtuvx)
  call NStateOpt(RotMat,DDg)
else
  !write(6,*) 'Optimization Approach 2'
  call NStateOpt2(RotMat,GDMat,Gtuvx)
end if
VecName = 'CMS-PDFT'
call PrintMat('ROT_VEC',VecName,RotMat,lroots,lroots,7,16,'N')

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
