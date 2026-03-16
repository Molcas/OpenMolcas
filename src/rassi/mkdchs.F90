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
! Copyright (C) 2021, Bruno Tenorio                                    *
!***********************************************************************

subroutine MKDCHS(IFSBTAB1,IFSBTAB2,ISSTAB,MAPORB,DET1,DET2,IF20,IF02,NDCHSM,DCHSM,OrbTab)
! Given two CI expansions, using a biorthonormal set of SD's,
! calculate the matrix elements relevant to DCH state intensities
! in the biorthonormal active orbital basis.
!
!  'I,J,|< N-2 | anni_right anni_right | N >|**2'

use Constants, only: Zero
use stdalloc, only: mma_allocate, mma_deallocate

implicit none
integer IFSBTAB1(*), IFSBTAB2(*)
integer ISSTAB(*), MAPORB(*), NDCHSM
real*8 DET1(*), DET2(*)
real*8 DCHSM(NDCHSM)
integer OrbTab(*)
integer NASHT, NASORB
real*8 GVAL, GAB, GBA
integer IAJB, IBJA
integer JORB, IORB
integer JORBA, JORBB, IORBA, IORBB
integer ITABS, JTABS, IJTABS
integer NSDCHSM
logical IF20, IF02
real*8, allocatable :: SDCHSM(:)

! Pick out nr of active orbitals from orbital table:
NASORB = ORBTAB(4)
NASHT = NASORB/2
NSDCHSM = NASORB*(NASORB-1)/2
call mma_allocate(SDCHSM,nSDCHSM,Label='SDCHSM')
SDCHSM(:) = Zero

call SDCHS(ORBTAB,ISSTAB,IFSBTAB1,IFSBTAB2,DET1,DET2,IF20,IF02,SDCHSM)

! Mapping from active spin-orbital to active orbital in external order.
! Note that these differ, not just because of the existence of two
! spin-orbitals for each orbital, but also because the active orbitals
! (external order) are grouped by symmetry and then RAS space, but the
! spin orbitals are grouped by subpartition.

IAJB = 0 ! dummy initialize
IBJA = 0 ! dummy initialize

do IORB=1,NASHT
  IORBA = 2*IORB-1
  IORBB = 2*IORB
  ITABS = MAPORB(IORBA)
  do JORB=1,NASHT
    JORBA = 2*JORB-1
    JORBB = 2*JORB
    JTABS = MAPORB(JORBA)
    GVAL = Zero
    if (IORB > JORB) then
      IAJB = ((IORBA-1)*(IORBA-2)/2)+JORBB
      IBJA = ((IORBB-1)*(IORBB-2)/2)+JORBA
    else if (JORB == IORB) then
      IAJB = ((IORBA-1)*(IORBA-2)/2)+JORBB
      IBJA = ((IORBB-1)*(IORBB-2)/2)+JORBA
      GAB = SDCHSM(IAJB)
      GBA = SDCHSM(IBJA)
      GVAL = GAB+GBA
    else if (IORB < JORB) then
      IBJA = ((JORBA-1)*(JORBA-2)/2)+IORBB
      IAJB = ((JORBB-1)*(JORBB-2)/2)+IORBA
    end if
    IJTABS = JTABS+NASHT*(ITABS-1)
    DCHSM(IJTABS) = GVAL**2
  end do
end do

call mma_deallocate(SDCHSM)

end subroutine MKDCHS
