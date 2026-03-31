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

subroutine SPIND2(ISYOP,MS2OP,IORBTAB,ISSTAB,IFSBTAB1,IFSBTAB4,PSI1,PSI4,SPD2)

use rassi_global_arrays, only: FSBANN1, FSBANN2, FSBANN3, FSBANN4
use Symmetry_Info, only: MUL
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: ISYOP, MS2OP, IORBTAB(*), ISSTAB(*), IFSBTAB1(*), IFSBTAB4(*)
real(kind=wp) :: PSI1(*), PSI4(*), SPD2(*)
integer(kind=iwp) :: IJKL, IMODE, ISMLAB, ISORB, ISPLAB, JISORB, JSMLAB, JSORB, JSPLAB, KLMS2, KLSORB, KLSYM, KOINFO, KSMLAB, &
                     KSORB, KSPLAB, LSMLAB, LSORB, LSPLAB, NASGEM, NASORB, ND1, ND2, ND3, ND4
real(kind=wp) :: COEFF, OVLP
real(kind=wp), allocatable :: ANN1(:), ANN2(:), ANN3(:), ANN4(:)
real(kind=wp), external :: OVERLAP_RASSI

! Nr of active spin-orbitals
NASORB = IORBTAB(4)
! Nr of active spin-orbital pairs:
NASGEM = (NASORB*(NASORB-1))/2
KOINFO = 19

do ISORB=2,NASORB
  ! Symmetry properties:
  ISMLAB = IORBTAB(KOINFO+1+8*(ISORB-1))
  ISPLAB = IORBTAB(KOINFO+3+8*(ISORB-1))
  ! Annihilate a single spin orbital, ISORB:
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
  do JSORB=1,ISORB-1
    ! Symmetry properties:
    JSMLAB = IORBTAB(KOINFO+1+8*(JSORB-1))
    JSPLAB = IORBTAB(KOINFO+3+8*(JSORB-1))
    ! Pair index:
    JISORB = ((ISORB-1)*(ISORB-2))/2+JSORB
    ! Annihilate once more, the spin orbital JSORB:
    IMODE = -1
    call FSBOP(IMODE,JSORB,IORBTAB,ISSTAB,FSBANN1,2)
    ND2 = FSBANN2(5)
    call mma_allocate(ANN2,ND2,Label='ANN2')
    ANN2(:) = Zero
    call PRIMSGM(IMODE,JSORB,IORBTAB,ISSTAB,FSBANN2,FSBANN1,COEFF,ANN2,ANN1)
    !TEST write(u6,*) ' The ANN2 wave function, with JSORB=',JSORB
    !TEST PRTHR = 0.01_wp
    !TEST call PRWVF(IORBTAB,ISSTAB,FSBANN2,PRTHR,ANN2)

    KLSYM = MUL(MUL(ISMLAB,JSMLAB),ISYOP)
    KLMS2 = MS2OP+ISPLAB+JSPLAB
    do LSORB=2,NASORB
      ! Symmetry properties:
      LSMLAB = IORBTAB(KOINFO+1+8*(LSORB-1))
      LSPLAB = IORBTAB(KOINFO+3+8*(LSORB-1))
      ! Annihilate a single spin orbital, LSORB:
      IMODE = -1
      call FSBOP(IMODE,LSORB,IORBTAB,ISSTAB,IFSBTAB4,4)
      ND4 = FSBANN4(5)
      COEFF = One
      call mma_allocate(ANN4,ND4,Label='ANN4')
      ANN4(:) = Zero
      call PRIMSGM(IMODE,LSORB,IORBTAB,ISSTAB,FSBANN4,IFSBTAB4,COEFF,ANN4,PSI4)
      !TEST write(u6,*) ' The ANN4 wave function, with LSORB=',LSORB
      !TEST PRTHR = 0.01_wp
      !TEST call PRWVF(IORBTAB,ISSTAB,FSBANN4,PRTHR,ANN4)
      do KSORB=1,LSORB-1
        OVLP = Zero
        ! Symmetry properties:
        KSMLAB = IORBTAB(KOINFO+1+8*(KSORB-1))
        KSPLAB = IORBTAB(KOINFO+3+8*(KSORB-1))
        if (MUL(KSMLAB,LSMLAB) /= KLSYM) cycle
        if (KSPLAB+LSPLAB /= KLMS2) cycle
        ! Pair index:
        KLSORB = ((LSORB-1)*(LSORB-2))/2+KSORB
        ! Annihilate once more, the spin orbital KSORB:
        IMODE = -1
        call FSBOP(IMODE,KSORB,IORBTAB,ISSTAB,FSBANN4,3)
        ND3 = FSBANN3(5)
        call mma_allocate(ANN3,ND3,Label='ANN3')
        ANN3(:) = Zero
        call PRIMSGM(IMODE,KSORB,IORBTAB,ISSTAB,FSBANN3,FSBANN4,COEFF,ANN3,ANN4)
        !TEST write(u6,*) ' The ANN3 wave function, with KSORB=',KSORB
        !TEST PRTHR = 0.01_wp
        !TEST call PRWVF(IORBTAB,ISSTAB,FSBANN3,PRTHR,ANN3)
        ! Compute the spin transition density matrix element:
        OVLP = OVERLAP_RASSI(FSBANN2,FSBANN3,ANN2,ANN3)
        !TEST write(u6,*) ' Their overlap:',OVLP
        IJKL = JISORB+NASGEM*(KLSORB-1)
        SPD2(IJKL) = SPD2(IJKL)+OVLP
        call mma_deallocate(ANN3)
        call mma_deallocate(FSBANN3)
      end do
      call mma_deallocate(ANN4)
      call mma_deallocate(FSBANN4)
    end do
    call mma_deallocate(ANN2)
    call mma_deallocate(FSBANN2)
  end do
  call mma_deallocate(ANN1)
  call mma_deallocate(FSBANN1)
end do

end subroutine SPIND2
