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
! Copyright (C) 2020, Bruno Tenorio                                    *
!***********************************************************************

subroutine SRTDM2(IORBTAB,ISSTAB,IFSBTAB1,IFSBTAB2,PSI1,PSI2,IF21,IF12,SRT2M)
! Calculates the 2-electron Dyson matrix between two states with
! N and N-1 electrons, defined as:
! IF12 D = < N-1 | anni_left anni_right anni_right | N >, or
! IF21 D = < N | anni_left anni_left anni_right | N-1 >
! reduced 2-electron tdm in the space of active spin-orbitals

use rassi_global_arrays, only: FSBANN1, FSBANN2, FSBANN3
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp, u6

#include "intent.fh"

implicit none
integer(kind=iwp), intent(in) :: IORBTAB(*), ISSTAB(*), IFSBTAB1(*), IFSBTAB2(*)
real(kind=wp), intent(in) :: PSI1(*), PSI2(*)
logical(kind=iwp), intent(in) :: IF21, IF12
real(kind=wp), intent(_OUT_) :: SRT2M(*)
integer(kind=iwp) :: IJL, IMODE, ISORB, JLSORB, JSORB, LSORB, NASORB, ND1, ND2, ND3
real(kind=wp) :: COEFF, OVLP
real(kind=wp), allocatable :: ANN1(:), ANN2(:), ANN3(:)
real(kind=wp), external :: OVERLAP_RASSI

NASORB = IORBTAB(4)

! IF12 = eliminte one electron to the left: < N-1 | anni_left (PSI1)
! and then eliminate two to the left (PSI2) anni_right anni_right | N >
if (IF12) then
  do ISORB=1,NASORB
    ! Annihilate a single spin orbital from PSI1 (N-1), ISORB:
    IMODE = -1
    call FSBOP(IMODE,ISORB,IORBTAB,ISSTAB,IFSBTAB1,1)
    ND1 = FSBANN1(5)
    COEFF = One
    call mma_allocate(ANN1,ND1,Label='ANN1')
    ANN1(:) = Zero
    call PRIMSGM(IMODE,ISORB,IORBTAB,ISSTAB,FSBANN1,IFSBTAB1,COEFF,ANN1,PSI1)
    !TEST write(u6,*) ' The ANN1 wave function, with ISORB=',ISORB
    !TEST PRTHR = 0.01_wp
    !TEST call PRWVF(IORBTAB,ISSTAB,FSBANN1,PRTHR,ANN1)
    do JSORB=1,NASORB
      ! Annihilate a single spin orbital from PSI2, the spin orbital JSORB:
      IMODE = -1
      call FSBOP(IMODE,JSORB,IORBTAB,ISSTAB,IFSBTAB2,2)
      ND2 = FSBANN2(5)
      COEFF = One
      call mma_allocate(ANN2,ND2,Label='ANN2')
      ANN2(:) = Zero
      call PRIMSGM(IMODE,JSORB,IORBTAB,ISSTAB,FSBANN2,IFSBTAB2,COEFF,ANN2,PSI2)
      !TEST write(u6,*) ' The ANN2 wave function, with JSORB=',JSORB
      !TEST PRTHR = 0.01_wp
      !TEST call PRWVF(IORBTAB,ISSTAB,FSBANN2,PRTHR,ANN2)
      do LSORB=1,NASORB
        ! Pair index J,L:
        JLSORB = (NASORB*(JSORB-1))+LSORB-1
        OVLP = Zero
        ! Annihilate another spin orbital from PSI2, LSORB:
        IMODE = -1
        call FSBOP(IMODE,LSORB,IORBTAB,ISSTAB,FSBANN2,3)
        ND3 = FSBANN3(5)
        COEFF = One
        call mma_allocate(ANN3,ND3,Label='ANN3')
        ANN3(:) = Zero
        if (JSORB /= LSORB) then
          call PRIMSGM(IMODE,LSORB,IORBTAB,ISSTAB,FSBANN3,FSBANN2,COEFF,ANN3,ANN2)
          ! Compute the spin transition density matrix element:
          OVLP = OVERLAP_RASSI(FSBANN1,FSBANN3,ANN1,ANN3)
        else
          OVLP = Zero
        end if
        IJL = ISORB+(NASORB*JLSORB)
        SRT2M(IJL) = OVLP
        call mma_deallocate(ANN3)
        call mma_deallocate(FSBANN3)
      end do
      call mma_deallocate(ANN2)
      call mma_deallocate(FSBANN2)
    end do
    call mma_deallocate(ANN1)
    call mma_deallocate(FSBANN1)
  end do
else if (IF21) then
  ! ################################################################################
  ! IF12 = Eliminate to the right (state 2)
  write(u6,*) 'Invalid state combination. Please, give PSI1=(N-1) and PSI2=(N)'
else
  write(u6,*) 'Invalid state combination in 2particle DYSON'
end if ! IF10 or IF01
! ################################################################################
end subroutine SRTDM2
