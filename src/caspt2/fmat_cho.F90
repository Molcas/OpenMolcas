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

subroutine FMAT_CHO(CMO,NCMO,FIAO,FAAO,HONE,NHONE,FIMO,NFIMO,FIFA,NFIFA)
! THIS ROUTINE IS USED IF THE TWO-ELECTRON INTEGRALS ARE
! REPRESENTED BY CHOLESKY VECTORS:
! TRANSFORM FOCK MATRICES COMPUTED BY TRACHO
! TO MO BASIS FOR USE IN CASPT2.

use constants, only: Zero, One
use stdalloc, only: mma_allocate, mma_deallocate
use caspt2_module, only: NBTRI, notri, NSYM, NBAS, NORB, NFRO
use definitions, only: iwp, wp
#ifdef _DEBUGPRINT_
use definitions, only: u6
#endif

implicit none
integer(kind=iwp), intent(in) :: NCMO, NHONE, NFIMO, NFIFA
real(kind=wp), intent(in) :: CMO(NCMO)
real(kind=wp), intent(in) :: FIAO(NBTRI), FAAO(NBTRI)
real(kind=wp), intent(in) :: HONE(NHONE)
real(kind=wp), intent(out) :: FIMO(NFIMO), FIFA(NFIFA)
real(kind=wp), allocatable :: SCR1(:), SCR2(:), SCR3(:)
integer(kind=iwp) I, IFAO, IJ, IOFMO, ISYM, J, LSC, LSCI, NB, NBBMX, NBBT, NBOMX, NF, NO, NO_X, NOOMX
#ifdef _DEBUGPRINT_
integer(kind=iwp) ISTLT
real(kind=wp), allocatable :: FAMO(:)
#endif

FIMO(:) = Zero
FIFA(:) = Zero ! initially used as FAMO

NBBT = 0
NBBMX = 0
NBOMX = 0
NOOMX = 0
do ISYM=1,NSYM
  NB = NBAS(ISYM)
  NO = NORB(ISYM)
  NBBT = NBBT+(NB*(NB+1))/2
  NBBMX = max(NBBMX,NB*NB)
  NBOMX = max(NBOMX,NB*NO)
  NOOMX = max(NOOMX,NO*NO)
end do

call mma_allocate(SCR1,NBBMX,LABEL='SCR1')
call mma_allocate(SCR2,NBOMX,LABEL='SCR2')
call mma_allocate(SCR3,NOOMX,LABEL='SCR3')

IFAO = 1
IOFMO = 0
LSC = 1
do ISYM=1,NSYM
  NB = NBAS(ISYM)
  if (NB == 0) cycle
  NO = NORB(ISYM)
  NO_X = max(1,NO)
  NF = NFRO(ISYM)
  LSCI = LSC+NF*NB
  ! The inactive Fock matrix:
  call SQUARE(FIAO(IFAO),SCR1,NB,1,NB)
  call DGEMM_('N','N',NB,NO,NB,One,SCR1,NB,CMO(LSCI),NB,Zero,SCR2,NB)
  call DGEMM_('T','N',NO,NO,NB,One,CMO(LSCI),NB,SCR2,NB,Zero,SCR3,NO_X)
  IJ = 0
  do I=1,NO
    do J=1,I
      IJ = IJ+1
      FIMO(IOFMO+IJ) = SCR3(I+NO*(J-1))
    end do
  end do
  ! The active Fock matrix:
  call SQUARE(FAAO(IFAO),SCR1,NB,1,NB)
  call DGEMM_('N','N',NB,NO,NB,One,SCR1,NB,CMO(LSCI),NB,Zero,SCR2,NB)
  call DGEMM_('T','N',NO,NO,NB,One,CMO(LSCI),NB,SCR2,NB,Zero,SCR3,NO_X)
  IJ = 0
  do I=1,NO
    do J=1,I
      IJ = IJ+1
      FIFA(IOFMO+IJ) = SCR3(I+NO*(J-1))
    end do
  end do
  IFAO = IFAO+(NB*(NB+1))/2
  IOFMO = IOFMO+(NO*(NO+1))/2
  LSC = LSC+NB**2
end do

call mma_deallocate(SCR1)
call mma_deallocate(SCR2)
call mma_deallocate(SCR3)

FIMO(1:NoTri) = FIMO(:)+HONE(:)
FIFA(1:NoTri) = FIMO(:)+FIFA(:)

#ifdef _DEBUGPRINT_
write(6,*) '      INACTIVE FOCK MATRIX IN MO BASIS'
ISTLT = 1
do ISYM=1,NSYM
  NO = NORB(ISYM)
  if (NO > 0) then
    write(6,'(6X,A,I2)') ' SYMMETRY SPECIES:',ISYM
    call TRIPRT(' ',' ',FIMO(ISTLT),NO)
    ISTLT = ISTLT+(NO*(NO+1))/2
  end if
end do

call mma_allocate(FAMO,nFIMO,Label='FAMO')
FAMO(:) = FIFA(:)-FIMO(:)
write(6,*) '        ACTIVE FOCK MATRIX IN MO BASIS'
ISTLT = 1
do ISYM=1,NSYM
  NO = NORB(ISYM)
  if (NO > 0) then
    write(6,'(6X,A,I2)') ' SYMMETRY SPECIES:',ISYM
    call TRIPRT(' ',' ',FAMO(ISTLT),NO)
    ISTLT = ISTLT+(NO*(NO+1))/2
  end if
end do
call mma_deallocate(FAMO)

write(6,*) '      TOTAL FOCK MATRIX IN MO BASIS'
ISTLT = 1
do ISYM=1,NSYM
  if (NORB(ISYM) > 0) then
    write(6,'(6X,A,I2)') ' SYMMETRY SPECIES:',ISYM
    call TRIPRT(' ',' ',FIFA(ISTLT),NORB(ISYM))
    ISTLT = ISTLT+NORB(ISYM)*(NORB(ISYM)+1)/2
  end if
end do

#endif

end subroutine FMAT_CHO
