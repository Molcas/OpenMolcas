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
! Copyright (C) 2018, Jesper Norell                                    *
!***********************************************************************

subroutine MKDYSORB(IORBTAB,ISSTAB,IFSBTAB1,IFSBTAB2,PSI1,PSI2,IF10,IF01,DYSAMP,DYSCOF)

use Constants, only: One, Zero
use stdalloc, only: mma_allocate, mma_deallocate
use rassi_global_arrays, only: FSBANN1, FSBANN2
use Definitions, only: u6

implicit none
integer IORBTAB(*)
integer ISSTAB(*)
integer IFSBTAB1(*), IFSBTAB2(*)
real*8 PSI1(*), PSI2(*)
logical IF10, IF01
real*8 DYSAMP, DYSCOF(*)
real*8 COEFF, OVLP
real*8, external :: OVERLAP_RASSI
integer NASORB
integer IMODE, ISORB
integer NDETS1, NDETS2
integer JSORB
real*8, allocatable :: ANN1(:), ANN2(:)

! +++ J. Norell 12/7 - 2018
! Calculates the Dyson orbital between two states with
! N and N-1 electrons, defined as:
! D = < N-1 | anni_right | N >, or
! D = < N | anni_left | N-1 >

! Nr of active spin-orbitals
NASORB = IORBTAB(4)
DYSAMP = Zero
do ISORB=1,NASORB
  DYSCOF(ISORB) = Zero
end do

if (IF10) then
  ! IF10 = Eliminate to the left (state 1)

  ! Loop over all spin orbitals ISORB:
  do ISORB=1,NASORB
    OVLP = Zero

    ! Annihilate a single orbital:
    COEFF = One
    IMODE = -1
    call FSBOP(IMODE,ISORB,IORBTAB,ISSTAB,IFSBTAB1,1)
    NDETS1 = FSBANN1(5)
    call mma_allocate(ANN1,NDETS1,Label='ANN1')
    ANN1(:) = Zero
    call PRIMSGM(IMODE,ISORB,IORBTAB,ISSTAB,FSBANN1,IFSBTAB1,COEFF,ANN1,PSI1)

    ! Compute the coefficient as the overlap between the N-1 electron w.f.s
    OVLP = OVERLAP_RASSI(FSBANN1,IFSBTAB2,ANN1,PSI2)
    call mma_deallocate(ANN1)
    call mma_deallocate(FSBANN1)
    DYSCOF(ISORB) = OVLP

    ! Collect the squared norm of the Dyson orbital
    DYSAMP = DYSAMP+OVLP*OVLP

  end do ! ISORB LOOP

else if (IF01) then
  ! IF01 = Eliminate to the right (state 2)

  ! Loop over all spin orbitals JSORB:
  do JSORB=1,NASORB
    OVLP = Zero

    ! Annihilate a single orbital:
    COEFF = One
    IMODE = -1
    call FSBOP(IMODE,JSORB,IORBTAB,ISSTAB,IFSBTAB2,2)
    NDETS2 = FSBANN2(5)
    ! BRN
    call mma_allocate(ANN2,NDETS2,Label='ANN2')
    ANN2(:) = Zero
    call PRIMSGM(IMODE,JSORB,IORBTAB,ISSTAB,FSBANN2,IFSBTAB2,COEFF,ANN2,PSI2)

    ! Compute the coefficient as the overlap between the N-1 electron w.f.s
    OVLP = OVERLAP_RASSI(IFSBTAB1,FSBANN2,PSI1,ANN2)
    call mma_deallocate(ANN2)
    call mma_deallocate(FSBANN2)
    DYSCOF(JSORB) = OVLP

    ! Collect the squared norm of the Dyson orbital
    DYSAMP = DYSAMP+OVLP*OVLP

  end do ! JSORB LOOP

else
  write(u6,*) 'Invalid state combination in MKDYSORB'
  write(u6,*) '(No such Dyson orbital can exist!)'

end if ! IF10 or IF01

! The eventual PES amplitude is given by the squared norm,
! but for transformation of the D_ij elements we need to remove the
! square for now
DYSAMP = sqrt(DYSAMP)

end subroutine MKDYSORB
