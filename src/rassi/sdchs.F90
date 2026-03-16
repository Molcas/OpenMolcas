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

subroutine SDCHS(IORBTAB,ISSTAB,IFSBTAB1,IFSBTAB2,PSI1,PSI2,IF20,IF02,SDCHSM)
! Calculates the DCH  matrix elements between two states with
! N and N-1 electrons, defined as:
! IF02 D = < N-2 | anni_right anni_right | N >, or
! IF20 D = < N | anni_left anni_left  | N-2 >
! reduced 2-electron tdm in the space of active spin-orbitals

use stdalloc, only: mma_allocate, mma_deallocate
use rassi_global_arrays, only: FSBANN1, FSBANN2
use Constants, only: Zero, One
use Definitions, only: u6

implicit none
integer IORBTAB(*), ISSTAB(*)
integer IFSBTAB1(*), IFSBTAB2(*)
real*8 PSI1(*), PSI2(*), SDCHSM(*)
logical IF20, IF02
real*8 COEFF, OVLP
integer NASORB
integer IMODE
integer ISORB, JSORB, IJ
integer ND1, ND2
real*8, allocatable :: ANN1(:), ANN2(:)
real*8, external :: OVERLAP_RASSI

NASORB = IORBTAB(4)

! IF02 = eliminte one electron to the right: < N-2 | anni_right
! anni_right | N >
if (IF02) then

  do ISORB=1,NASORB
    ! Symmetry properties:
    !ISMLAB = IORBTAB(KOINFO+1+8*(ISORB-1))
    !ISPLAB = IORBTAB(KOINFO+3+8*(ISORB-1))
    ! Annihilate a single spin orbital from PSI2, the spin orbital ISORB:
    IMODE = -1
    call FSBOP(IMODE,ISORB,IORBTAB,ISSTAB,IFSBTAB2,1)
    ND1 = FSBANN1(5)
    COEFF = One
    call mma_allocate(ANN1,ND1,Label='ANN1')
    ANN1(:) = Zero
    call PRIMSGM(IMODE,ISORB,IORBTAB,ISSTAB,FSBANN1,IFSBTAB2,COEFF,ANN1,PSI2)

    do JSORB=1,ISORB-1
      ! Symmetry properties:
      !JSMLAB = IORBTAB(KOINFO+1+8*(JSORB-1))
      !JSPLAB = IORBTAB(KOINFO+3+8*(JSORB-1))
      ! Pair index J,L:

      OVLP = Zero
      ! Annihilate another spin orbital from PSI2, LSORB:
      IMODE = -1
      call FSBOP(IMODE,JSORB,IORBTAB,ISSTAB,FSBANN1,2)
      ND2 = FSBANN2(5)
      COEFF = One
      call mma_allocate(ANN2,ND2,Label='ANN2')
      ANN2(:) = Zero
      call PRIMSGM(IMODE,JSORB,IORBTAB,ISSTAB,FSBANN2,FSBANN1,COEFF,ANN2,ANN1)

      ! Compute the spin transition density matrix element:
      OVLP = OVERLAP_RASSI(IFSBTAB1,FSBANN2,PSI1,ANN2)

      IJ = ((ISORB-1)*(ISORB-2))/2+JSORB
      SDCHSM(IJ) = SDCHSM(IJ)+OVLP

      call mma_deallocate(ANN2)
      call mma_deallocate(FSBANN2)
    end do
    call mma_deallocate(ANN1)
    call mma_deallocate(FSBANN1)
  end do

else if (IF20) then
  ! ################################################################################
  ! IF02 = Eliminate to the right (state 2)
  write(u6,*) 'Invalid state combination. Please, give PSI1=(N-2) and PSI2=(N)'
else
  write(u6,*) 'Invalid state combination in DCH states'
end if

end subroutine SDCHS
