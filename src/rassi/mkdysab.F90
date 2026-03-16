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
! Copyright (C) 1989, Per Ake Malmqvist                                *
!               2018, Jesper Norell                                    *
!               2020, Bruno Tenorio                                    *
!***********************************************************************

!*********************************************************************
! Modified from MKTAB to MKDYSAB by Jesper Norell, 2018
!  SUBROUTINE MKDYSAB
!  PURPOSE: CALCULATE DYSON ORBITAL COEFFICIENTS FOR CI EXPANSIONS IN
!  BIORTHONORMAL ORBITAL BASE A,
!  IN ANALOGUE TO MKTDAB FOR TRANSITION DENSITY MATRIX.
!*********************************************************************
!  MODIFIED BY BRUNO TENORIO TO ADDRESS SYMMETRY
!  SEPTEMBER 2020
!*********************************************************************

subroutine MKDYSAB(DYSCOF,DYSAB)

use Constants, only: Zero
use stdalloc, only: mma_allocate, mma_deallocate
use Symmetry_Info, only: nSym => nIrrep
use rassi_data, only: NASHT, NASH, NISH, NOSH

implicit none
real*8 DYSCOF(*), DYSAB(*)
integer :: IOFFA(8)
real*8 GAA, GBB, OVLP
integer IORB, ISORB, I, IOFFTD, ISY, II, IPOS, ICOFF, ISY1, NO1, NA1, NI1
real*8, allocatable :: DYSCOF2(:)

!+++BRN Create a scalar spin summed Dyson coefficients DYSCOF2
!Alpha and beta contributions are added up here
call mma_allocate(DYSCOF2,NASHT,Label='DYSCOF2')
do IORB=1,NASHT
  ISORB = 2*IORB-1
  GAA = DYSCOF(ISORB)
  GBB = DYSCOF(ISORB+1)
  OVLP = GBB+GAA
  !normally GAA gives just zeros...
  DYSCOF2(IORB) = OVLP
end do
! IOFFA=NR OF ACTIVE ORBITALS IN PREVIOUS SYMMETRY BLOCKS.
IOFFA(1) = 0
do I=1,NSYM-1
  IOFFA(I+1) = IOFFA(I)+NASH(I)
end do

! CONTRIBUTION FROM INACTIVE ORBITALS:
! (By definition 0 for Dyson orbitals,
! but we need to fill out the full vector for easier
! transformation.)
IOFFTD = 0
do ISY=1,NSYM
  if (NISH(ISY) /= 0) then
    II = 0
    do I=1,NISH(ISY)
      II = II+1
      IPOS = IOFFTD+II
      DYSAB(IPOS) = Zero
    end do
    IOFFTD = IOFFTD+NOSH(ISY)
  end if
end do
! THEN ADD CONTRIBUTION FROM ACTIVE SPACE.
IOFFTD = 0
ICOFF = 1
do ISY1=1,NSYM
  NO1 = NOSH(ISY1)
  if (NO1 == 0) cycle
  NA1 = NASH(ISY1)
  if (NA1 /= 0) then
    NI1 = NISH(ISY1)
    do I=1,NA1
      II = NI1+I
      IPOS = IOFFTD+II
      DYSAB(IPOS) = DYSCOF2(ICOFF)
      ICOFF = ICOFF+1
    end do
  end if
  IOFFTD = IOFFTD+NO1
end do
call mma_deallocate(DYSCOF2)

end subroutine MKDYSAB
