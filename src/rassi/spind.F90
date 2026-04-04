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

subroutine SPIND(ISYOP,MS2OP,IORBTAB,ISSTAB,IFSBTAB1,IFSBTAB2,PSI1,PSI2,SPD12)

use rassi_global_arrays, only: FSBANN1, FSBANN2
use Symmetry_Info, only: MUL
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: ISYOP, MS2OP, IORBTAB(*), ISSTAB(*), IFSBTAB1(*), IFSBTAB2(*)
real(kind=wp) :: PSI1(*), PSI2(*), SPD12(*)
integer(kind=iwp) :: IMODE, ISMLAB, ISORB, ISPLAB, JMS2, JSMLAB, JSORB, JSPLAB, JSYM, KOINFO, NASORB, NDETS1, NDETS2
real(kind=wp) :: COEFF, OVLP
real(kind=wp), allocatable :: ANN1(:), ANN2(:)
real(kind=wp), external :: OVERLAP_RASSI

!TEST write(u6,*) ' Test prints in SPIND.'
! Nr of active spin-orbitals
NASORB = IORBTAB(4)
! Loop over pairs of spin orbitals ISORB,JSORB:
KOINFO = 19
do ISORB=1,NASORB
  ISMLAB = IORBTAB(KOINFO+1+8*(ISORB-1))
  !UNUSED ISOIND = IORBTAB(KOINFO+2+8*(ISORB-1))
  ISPLAB = IORBTAB(KOINFO+3+8*(ISORB-1))
  ! Annihilate a single orbital:
  COEFF = One
  IMODE = -1
  call FSBOP(IMODE,ISORB,IORBTAB,ISSTAB,IFSBTAB1,1)
  NDETS1 = FSBANN1(5)
  call mma_allocate(ANN1,NDETS1,Label='ANN1')
  ANN1(:) = Zero
  call PRIMSGM(IMODE,ISORB,IORBTAB,ISSTAB,FSBANN1,IFSBTAB1,COEFF,ANN1,PSI1)
  !TEST write(u6,*) ' Prior to call to PRIMSGM.'
  !TEST write(u6,*) ' FS block structure at FSBANN1:'
  !TEST call PRFSBTAB(FSBANN1)
  !TEST write(u6,*) ' Wave function ANN1 after PRIMSGM:'
  !TEST PRTHR = 0.01_wp
  !TEST call PRWVF(IORBTAB,ISSTAB,FSBANN1,PRTHR,ANN1)
  ! Compute those ANN2 wave functions, which have the correct properties:
  JSYM = MUL(ISMLAB,ISYOP)
  JMS2 = MS2OP+ISPLAB
  do JSORB=1,NASORB
    OVLP = Zero
    JSMLAB = IORBTAB(KOINFO+1+8*(JSORB-1))
    !UNUSED JSOIND = IORBTAB(KOINFO+2+8*(JSORB-1))
    JSPLAB = IORBTAB(KOINFO+3+8*(JSORB-1))
    if ((JSMLAB == JSYM) .and. (JSPLAB == JMS2)) then
      COEFF = One
      IMODE = -1
      call FSBOP(IMODE,JSORB,IORBTAB,ISSTAB,IFSBTAB2,2)
      NDETS2 = FSBANN2(5)
      call mma_allocate(ANN2,NDETS2,Label='ANN2')
      ANN2(:) = Zero
      call PRIMSGM(IMODE,JSORB,IORBTAB,ISSTAB,FSBANN2,IFSBTAB2,COEFF,ANN2,PSI2)

      ! Compute the spin transition density matrix element:
      OVLP = OVERLAP_RASSI(FSBANN1,FSBANN2,ANN1,ANN2)
      !TEST write(u6,*) ' Their overlap:',OVLP
      call mma_deallocate(ANN2)
      call mma_deallocate(FSBANN2)
    end if
    SPD12(ISORB+NASORB*(JSORB-1)) = OVLP
  end do
  call mma_deallocate(ANN1)
  call mma_deallocate(FSBANN1)
end do

end subroutine SPIND
