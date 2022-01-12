!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************
!
! Some incomplete documentation
!
! Vlct       : logical flag to indicate that the velocity integrals will
!              be computed
! lRel       : logical flag to indicate that the 1-electron Darwin and
!              mass-velocity integrals should be computed
! UnNorm     : logical flag to indicate that the primitive basis
!              functions should not be normalized
! lSchw      : logical flag to indicate that the Cauchy-Schwartz
!              inequality should be used in the prescreening
! lAMFI      : logical flag for Atomic Mean Field Integrals
! GIAO       : integrals for first derivative with respect to B
! Do_FckInt  : logical flag for FckInt
! Do_GuessOrb: logical flag for GuessOrb
! EMFR       : electromagnetic field radiation
! FNMC       : finite nuclear mass correction

module Logical_Info

use Definitions, only: iwp

implicit none
private

integer(kind=iwp), parameter :: lLen = 19 ! number of elements
logical(kind=iwp) :: Align_Only = .false., &
                     Do_Align = .true., &
                     Do_FckInt = .true., &
                     Do_GuessOrb = .true., &
                     DoFMM = .false., &
                     EMFR = .false., &
                     FNMC = .false., &
                     GIAO = .false., &
                     lAMFI = .false., &
                     lDOWNONLY = .false., &
                     lMXTC = .false., &
                     lRel = .false., &
                     lRP = .false., &
                     lRP_Post = .false., &
                     lSchw = .true., &
                     lUPONLY = .false., &
                     NEMO = .false., &
                     UnNorm = .false., &
                     Vlct = .true.

public :: Align_Only, DoFMM, Do_Align, Do_FckInt, Do_GuessOrb, EMFR, FNMC, GIAO, Logical_Info_Dmp, Logical_Info_Get, NEMO, UnNorm, &
          Vlct, lAMFI, lDOWNONLY, lMXTC, lRP, lRP_Post, lRel, lSchw, lUPONLY

!interface
!  subroutine Put_iArray(Label,data,nData)
!    character*(*) Label
!    integer nData
!    integer data(nData)
!  end subroutine Put_iArray
!  subroutine Get_iArray(Label,data,nData)
!    character*(*) Label
!    integer nData
!    integer data(nData)
!  end subroutine Get_iArray
!end interface

contains

subroutine Logical_Info_Dmp()

  use stdalloc, only: mma_allocate, mma_deallocate

  integer(kind=iwp), allocatable :: iDmp(:)

  call mma_allocate(iDmp,lLen,Label='iDmp:Logical')

  iDmp(01) = merge(1,0,Vlct)
  iDmp(02) = merge(1,0,lRel)
  iDmp(03) = merge(1,0,UnNorm)
  iDmp(04) = merge(1,0,lSchw)
  iDmp(05) = merge(1,0,lAMFI)
  iDmp(06) = merge(1,0,NEMO)
  iDmp(07) = merge(1,0,Do_GuessOrb)
  iDmp(08) = merge(1,0,Do_FckInt)
  iDmp(09) = merge(1,0,Align_Only)
  iDmp(10) = merge(1,0,DoFMM)
  iDmp(11) = merge(1,0,lRP)
  iDmp(12) = merge(1,0,lRP_Post)
  iDmp(13) = merge(1,0,Do_Align)
  iDmp(14) = merge(1,0,EMFR)
  iDmp(15) = merge(1,0,GIAO)
  iDmp(16) = merge(1,0,lUPONLY)
  iDmp(17) = merge(1,0,lDOWNONLY)
  iDmp(18) = merge(1,0,FNMC)
  iDmp(19) = merge(1,0,lMXTC)

  call Put_iArray('Logical_Info',iDmp,lLen)
  call mma_deallocate(iDmp)

end subroutine Logical_Info_Dmp

subroutine Logical_Info_Get()

  use stdalloc, only: mma_allocate, mma_deallocate

  integer(kind=iwp), allocatable :: iDmp(:)

  call mma_allocate(iDmp,lLen,Label='iDmp:Logical')
  call Get_iArray('Logical_Info',iDmp,lLen)

  Vlct = iDmp(1) > 0
  lRel = iDmp(2) > 0
  UnNorm = iDmp(3) > 0
  lSchw = iDmp(4) > 0
  lAMFI = iDmp(5) > 0
  NEMO = iDmp(6) > 0
  Do_GuessOrb = iDmp(7) > 0
  Do_FckInt = iDmp(8) > 0
  Align_Only = iDmp(9) > 0
  DoFMM = iDmp(10) > 0
  lRP = iDmp(11) > 0
  lRP_Post = iDmp(12) > 0
  Do_Align = iDmp(13) > 0
  EMFR = iDmp(14) > 0
  GIAO = iDmp(15) > 0
  lUPONLY = iDmp(16) > 0
  lDOWNONLY = iDmp(17) > 0
  FNMC = iDmp(18) > 0
  lMXTC = iDmp(19) > 0

  call mma_deallocate(iDmp)

end subroutine Logical_Info_Get

end module Logical_Info
