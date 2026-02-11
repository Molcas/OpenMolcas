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
      Subroutine CMSRot(TUVX)
! ****************************************************************
! history:                                                       *
! Jie J. Bao, on Aug. 06, 2020, created this file.               *
! ****************************************************************
      use stdalloc, only : mma_allocate, mma_deallocate
      use CMS, only: CMSNotConverged
      use rasscf_global, only: NACPR2, CMSStartMat, lRoots, NAC
      use PrintLevel, only: USUAL
      use output_ras, only: LF,IPRLOC
      Implicit None

#include "warnings.h"

      Real*8,DIMENSION(NACPR2)::TUVX
      CHARACTER(len=16)::VecName
      Real*8,DIMENSION(:,:,:,:),Allocatable::Gtuvx
      Real*8,DIMENSION(:,:,:,:),Allocatable::DDG
      Real*8,DIMENSION(:,:,:),Allocatable::GDMat
      Real*8,DIMENSION(:,:),Allocatable::RotMat
      Integer iPrLev

!     Allocating Memory
      CALL mma_allocate(GDMat,LRoots*(LRoots+1)/2,NAC,NAC)
      CALL mma_allocate(RotMat,lRoots,lRoots)
      CALL mma_allocate(Gtuvx,NAC,NAC,NAC,NAC)
      CALL mma_allocate(DDG,lRoots,lRoots,lRoots,lRoots)

      IPRLEV=IPRLOC(6)

!     printing header
      IF(IPRLEV.ge.USUAL) THEN
      write(LF,*)
      write(LF,*)
      write(LF,*) '    CMS INTERMEDIATE-STATE OPTIMIZATION'
      END IF
      IF(trim(CMSStartMat).eq.'XMS') THEN
       CALL ReadMat('ROT_VEC',VecName,RotMat,lroots,lroots,7,16,'N')
      ELSE
       CALL ReadMat(trim(CMSStartMat),VecName,RotMat,lroots,lroots,     &
     &              len_trim(CMSStartMat),16,'N')
      END IF
      IF(IPRLEV.ge.USUAL)                                               &
     &  CALL CMSHeader(trim(CMSStartMat),len_trim(CMSStartMat))


      CALL LoadGtuvx(TUVX,Gtuvx)

      CMSNotConverged=.false.
      CALL GetGDMat(GDMat)
      IF(lRoots.lt.NAC) THEN
!       write(6,*)"Optimization Approach 1"
       CALL GetDDgMat(DDg,GDMat,Gtuvx)
       CALL NStateOpt(RotMat,DDg)
      ELSE
!       write(6,*)"Optimization Approach 2"
       CALL NStateOpt2(RotMat,GDMat,Gtuvx)
      END IF
      VecName='CMS-PDFT'
      CALL PrintMat('ROT_VEC',VecName,RotMat,lroots,lroots,7,16,'N')

!     Deallocating Memory
      CALL mma_deallocate(GDMat)
      CALL mma_deallocate(RotMat)
      CALL mma_deallocate(Gtuvx)
      CALL mma_deallocate(DDg)
      IF(CMSNotConverged) THEN
       Call WarningMessage(2,'CMS Intermediate States Not Converged')
       Call Quit(_RC_NOT_CONVERGED_)
      END IF
      End Subroutine CMSRot
