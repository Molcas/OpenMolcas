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

use Index_Functions, only: nTri_Elem
use caspt2_module, only: NBAS, NBTRI, NFRO, NORB, notri, NSYM
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

implicit none
integer(kind=iwp), intent(in) :: NCMO, NHONE, NFIMO, NFIFA
real(kind=wp), intent(in) :: CMO(NCMO), FIAO(NBTRI), FAAO(NBTRI), HONE(NHONE)
real(kind=wp), intent(out) :: FIMO(NFIMO), FIFA(NFIFA)
integer(kind=iwp) :: I, IFAO, IJ, IOFMO, ISYM, LSC, LSCI, NB, NBBMX, NBOMX, NF, NO, NO_X, NOOMX
real(kind=wp), allocatable :: SCR1(:), SCR2(:), SCR3(:)
#ifdef _DEBUGPRINT_
integer(kind=iwp) :: ISTLT
real(kind=wp), allocatable :: FAMO(:)
#endif

FIMO(:) = Zero
FIFA(:) = Zero ! initially used as FAMO

NBBMX = maxval(NBAS(1:NSYM)**2)
NBOMX = maxval(NBAS(1:NSYM)*NORB(1:NSYM))
NOOMX = maxval(NORB(1:NSYM)**2)

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
    FIMO(IOFMO+IJ+1:IOFMO+IJ+I) = SCR3(I:I+NO*(I-1):NO)
    IJ = IJ+I
  end do
  ! The active Fock matrix:
  call SQUARE(FAAO(IFAO),SCR1,NB,1,NB)
  call DGEMM_('N','N',NB,NO,NB,One,SCR1,NB,CMO(LSCI),NB,Zero,SCR2,NB)
  call DGEMM_('T','N',NO,NO,NB,One,CMO(LSCI),NB,SCR2,NB,Zero,SCR3,NO_X)
  IJ = 0
  do I=1,NO
    FIFA(IOFMO+IJ+1:IOFMO+IJ+I) = SCR3(I:I+NO*(I-1):NO)
    IJ = IJ+I
  end do
  IFAO = IFAO+nTri_Elem(NB)
  IOFMO = IOFMO+nTri_Elem(NO)
  LSC = LSC+NB**2
end do

call mma_deallocate(SCR1)
call mma_deallocate(SCR2)
call mma_deallocate(SCR3)

FIMO(1:NoTri) = FIMO(:)+HONE(:)
FIFA(1:NoTri) = FIMO(:)+FIFA(:)

#ifdef _DEBUGPRINT_
write(u6,*) '      INACTIVE FOCK MATRIX IN MO BASIS'
ISTLT = 1
do ISYM=1,NSYM
  NO = NORB(ISYM)
  if (NO > 0) then
    write(u6,'(6X,A,I2)') ' SYMMETRY SPECIES:',ISYM
    call TRIPRT(' ',' ',FIMO(ISTLT),NO)
    ISTLT = ISTLT+nTri_Elem(NO)
  end if
end do

call mma_allocate(FAMO,nFIMO,Label='FAMO')
FAMO(:) = FIFA(:)-FIMO(:)
write(u6,*) '        ACTIVE FOCK MATRIX IN MO BASIS'
ISTLT = 1
do ISYM=1,NSYM
  NO = NORB(ISYM)
  if (NO > 0) then
    write(u6,'(6X,A,I2)') ' SYMMETRY SPECIES:',ISYM
    call TRIPRT(' ',' ',FAMO(ISTLT),NO)
    ISTLT = ISTLT+nTri_Elem(NO)
  end if
end do
call mma_deallocate(FAMO)

write(u6,*) '      TOTAL FOCK MATRIX IN MO BASIS'
ISTLT = 1
do ISYM=1,NSYM
  if (NORB(ISYM) > 0) then
    write(u6,'(6X,A,I2)') ' SYMMETRY SPECIES:',ISYM
    call TRIPRT(' ',' ',FIFA(ISTLT),NORB(ISYM))
    ISTLT = ISTLT+nTri_Elem(NORB(ISYM))
  end if
end do

#endif

end subroutine FMAT_CHO
