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

private
public :: Vlct, lRel, UnNorm, lSchw, lAMFI, NEMO, Do_GuessOrb, Do_FckInt, Align_Only, DoFMM, lRP, lRP_Post, EMFR, Do_Align, GIAO, &
          lUPONLY, lDOWNONLY, FNMC, lMXTC, Logical_Info_Dmp, Logical_Info_Get

#include "stdalloc.fh"
integer i
logical :: Vlct = .true.
logical :: lRel = .false.
logical :: UnNorm = .false.
logical :: lSchw = .true.
logical :: lAMFI = .false.
logical :: NEMO = .false.
logical :: Do_GuessOrb = .true.
logical :: Do_FckInt = .true.
logical :: Align_Only = .false.
logical :: DoFMM = .false.
logical :: lRP = .false.
logical :: lRP_Post = .false.
logical :: Do_Align = .true.
logical :: EMFR = .false.
logical :: GIAO = .false.
logical :: lUPONLY = .false.
logical :: lDOWNONLY = .false.
logical :: FNMC = .false.
logical :: lMXTC = .false.

interface
  subroutine Put_iArray(Label,data,nData)
    character*(*) Label
    integer nData
    integer data(nData)
  end subroutine Put_iArray
  subroutine Get_iArray(Label,data,nData)
    character*(*) Label
    integer nData
    integer data(nData)
  end subroutine Get_iArray
end interface

contains

subroutine Logical_Info_Dmp()

  integer, allocatable :: iDmp(:)
  integer :: Len = 19

  call mma_allocate(iDmp,Len,Label='iDmp:Logical')

  i = 0
  if (Vlct) i = 1
  iDmp(01) = i
  i = 0
  if (lRel) i = 1
  iDmp(02) = i
  i = 0
  if (UnNorm) i = 1
  iDmp(03) = i
  i = 0
  if (lSchw) i = 1
  iDmp(04) = i
  i = 0
  if (lAMFI) i = 1
  iDmp(05) = i
  i = 0
  if (NEMO) i = 1
  iDmp(06) = i
  i = 0
  if (Do_GuessOrb) i = 1
  iDmp(07) = i
  i = 0
  if (Do_FckInt) i = 1
  iDmp(08) = i
  i = 0
  if (Align_Only) i = 1
  iDmp(09) = i
  i = 0
  if (DoFMM) i = 1
  iDmp(10) = i
  i = 0
  if (lRP) i = 1
  iDmp(11) = i
  i = 0
  if (lRP_Post) i = 1
  iDmp(12) = i
  i = 0
  if (Do_Align) i = 1
  iDmp(13) = i
  i = 0
  if (EMFR) i = 1
  iDmp(14) = i
  i = 0
  if (GIAO) i = 1
  iDmp(15) = i
  i = 0
  if (lUPONLY) i = 1
  iDmp(16) = i
  i = 0
  if (lDOWNONLY) i = 1
  iDmp(17) = i
  i = 0
  if (FNMC) i = 1
  iDmp(18) = i
  i = 0
  if (lMXTC) i = 1
  iDmp(19) = i

  call Put_iArray('Logical_Info',iDmp,Len)
  call mma_deallocate(iDmp)

end subroutine Logical_Info_Dmp

subroutine Logical_Info_Get()

  integer, allocatable :: iDmp(:)
  integer :: Len = 19

  call mma_allocate(iDmp,Len,Label='iDmp:Logical')
  call Get_iArray('Logical_Info',iDmp,Len)

  Vlct = iDmp(01) == 1
  lRel = iDmp(02) == 1
  UnNorm = iDmp(03) == 1
  lSchw = iDmp(04) == 1
  lAMFI = iDmp(05) == 1
  NEMO = iDmp(06) == 1
  Do_GuessOrb = iDmp(07) == 1
  Do_FckInt = iDmp(08) == 1
  Align_Only = iDmp(09) == 1
  DoFMM = iDmp(10) == 1
  lRP = iDmp(11) == 1
  lRP_Post = iDmp(12) == 1
  Do_Align = iDmp(13) == 1
  EMFR = iDmp(14) == 1
  GIAO = iDmp(15) == 1
  lUPONLY = iDmp(16) == 1
  lDOWNONLY = iDmp(17) == 1
  FNMC = iDmp(18) == 1
  lMXTC = iDmp(19) == 1

  call mma_deallocate(iDmp)

end subroutine Logical_Info_Get

end module Logical_Info
