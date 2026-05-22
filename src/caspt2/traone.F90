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

subroutine TRAONE(CMO,NCMO,HONE,nHONE)
! Objective: Transformation of one-electron integrals
! (effective one electron Hamiltonian) for CASPT2.

use OneDat, only: sNoNuc, sNoOri
use PrintLevel, only: VERBOSE
use caspt2_global, only: iPrGlb
use caspt2_module, only: ERFSELF, nBas, nBMX, nBSqT, nBTri, nDel, nFro, nFroT, nOrb, nOTri, nSym, PotNuc, RFPert
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Two, Half
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: NCMO, nHONE
real(kind=wp), intent(in) :: CMO(NCMO)
real(kind=wp), intent(inout) :: HONE(nHONE)
integer(kind=iwp) :: I, iAO, IB, ICMO, ICOMP, IERR, IFTEST, IJ, IMO, IOFF, IOPT, IRC, ISTLT, ISTMO, ISTSQ, ISYLBL, ISYM, Keep(8), &
                     NB, nBasXX(8), NF, NSYMXX, nTemp, NWTMP
real(kind=wp) :: ECORE, EONE, ETWO, ExFac
logical(kind=iwp) :: Found, iSquar
character(len=8) :: Label
real(kind=wp), allocatable :: Temp(:), WDLT(:), WDSQ(:), WFLT(:), WFMO(:), WTMP(:)
real(kind=wp), external :: DDot_
#include "warnings.h"

#ifdef _DEBUGPRINT_
IFTEST = 1
#else
IFTEST = 0
#endif

call GetOrd(IRC,iSquar,nSymXX,nBasXX,Keep)
if (IPRGLB >= VERBOSE) then
  if (iSquar) then
    write(u6,*) 'TRAONE OrdInt status: squared'
  else
    write(u6,*) 'TRAONE OrdInt status: non-squared'
  end if
end if
if (any(NBAS(1:NSYM) /= NBASXX(1:NSYM))) then
  IERR = 1
else
  IERR = 0
end if
if (IERR /= 0) then
  write(u6,*) '     *** ERROR IN SUBROUTINE TRAONE ***'
  write(u6,*) '          INCOMPATIBLE BASIS DATA'
  write(u6,*)
  write(u6,*) ' JOBIPH NR OF SYMM:',NSYM
  write(u6,*) ' JOBIPH NR OF BASIS FUNCTIONS/SYMM:'
  write(u6,'(1x,8I5)') (NBAS(I),I=1,NSYM)
  write(u6,*)
  write(u6,*) ' ORDINT NR OF SYMM:',NSYMXX
  write(u6,*) ' ORDINT NR OF BASIS FUNCTIONS/SYMM:'
  write(u6,'(1x,8I5)') (NBASXX(I),I=1,NSYMXX)
  call ABEND()
end if
! Allocate FLT,DLT, and DSQ.
call mma_allocate(WFLT,NBTRI,Label='WFLT')
! Read nuclear repulsion energy:
IRC = -1
IOPT = 0
ICOMP = 0
ISYLBL = 1
if (IFTEST /= 0) write(u6,*) ' GET POTNUC FROM RUNFILE'
call Get_dScalar('PotNuc',PotNuc)
if (IFTEST /= 0) write(u6,*) ' POTNUC:',POTNUC
! Read one-electron hamiltonian matrix into FLT.
IRC = -1
IOPT = ibset(ibset(0,sNoOri),sNoNuc)
ICOMP = 1
ISYLBL = 1
Label = 'OneHam'
if (IFTEST /= 0) write(u6,*) ' CALLING RDONE (ONEHAM)'
call RDONE(IRC,IOPT,Label,ICOMP,WFLT,ISYLBL)
if (IFTEST /= 0) write(u6,*) ' BACK FROM RDONE'
if (IRC /= 0) then
  write(u6,*) 'TRAONE Error: RDONE failed reading OneHam.'
  call Quit(_RC_IO_ERROR_READ_)
end if

if (IFTEST /= 0) then
  write(u6,*) '     TEST PRINTS FROM TRAONE.'
  write(u6,*) '     NAKED 1-EL HAMILTONIAN IN AO BASIS'
  ISTLT = 1
  do ISYM=1,NSYM
    if (NBAS(ISYM) > 0) then
      write(u6,'(6X,A,I2)') ' SYMMETRY SPECIES:',ISYM
      call TRIPRT(' ',' ',WFLT(ISTLT),NBAS(ISYM))
      ISTLT = ISTLT+NBAS(ISYM)*(NBAS(ISYM)+1)/2
    end if
  end do
end if

! If this is a perturbative reaction field calculation then
! modifiy the one-electron Hamiltonian by the reaction field and
! the nuclear attraction by the cavity self-energy

if (RFpert) then
  nTemp = sum(nBas(1:nSym)*(nBas(1:nSym)+1)/2)
  call mma_allocate(Temp,nTemp,Label='Temp')

  call f_Inquire('RUNOLD',Found)
  if (Found) call NameRun('RUNOLD')
  call Get_dScalar('RF Self Energy',ERFSelf)
  call Get_dArray('Reaction field',Temp,nTemp)
  if (Found) call NameRun('#Pop')
  PotNuc = PotNuc+ERFself
  WFLT(1:nTemp) = WFLT(1:nTemp)+Temp(:)

  call mma_deallocate(Temp)
  if (IFTEST /= 0) then
    write(u6,*) ' 1-EL HAMILTONIAN INCLUDING REACTION FIELD'
    ISTLT = 1
    do ISYM=1,NSYM
      if (NBAS(ISYM) > 0) then
        write(u6,'(6X,A,I2)') ' SYMMETRY SPECIES:',ISYM
        call TRIPRT(' ',' ',WFLT(ISTLT),NBAS(ISYM))
        ISTLT = ISTLT+NBAS(ISYM)*(NBAS(ISYM)+1)/2
      end if
    end do
  end if
end if

EONE = Zero
ETWO = Zero
! The following section is needed for frozen orbitals:
if (NFROT /= 0) then
  call mma_allocate(WDLT,NBTRI,LABEL='WDLT')
  call mma_allocate(WDSQ,NBSQT,LABEL='WDSQ')
  ! Compute the density matrix of the frozen orbitals
  ! The DLT matrix contains the same data as DSQ, but
  ! with symmetry blocks in lower triangular format, and
  ! with non-diagonal elements doubled.
  WDLT(:) = Zero
  WDSQ(:) = Zero
  ISTMO = 1
  ISTSQ = 1
  ISTLT = 1
  do ISYM=1,NSYM
    NF = NFRO(ISYM)
    NB = NBAS(ISYM)
    if (NB*NF > 0) then
      call DGEMM_('N','T',NB,NB,NF,Two,CMO(ISTMO),NB,CMO(ISTMO),NB,Zero,WDSQ(ISTSQ),NB)
      IJ = ISTLT-1
      do IB=1,NB
        WDLT(IJ+1:IJ+IB) = Two*WDSQ(ISTSQ+(IB-1)*NB:ISTSQ+(IB-1)*NB+IB-1)
        IJ = IJ+IB
        WDLT(IJ) = Half*WDLT(IJ)
      end do
    end if
    ISTMO = ISTMO+NB*NB
    ISTSQ = ISTSQ+NB*NB
    ISTLT = ISTLT+NB*(NB+1)/2
  end do

  ! One-electron contribution to the core energy.
  ! Note that FLT still contains only the naked
  !  one-electron hamiltonian.
  EONE = DDOT_(NBTRI,WDLT,1,WFLT,1)
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! Generate Fock-matrix for frozen orbitals
  ! and compute the total core energy
  ! Look out-- we temporarily allocate all available memory.

  ExFac = One
  call FTwo_Drv(nSym,nBas,nFro,KEEP,WDLT,WDSQ,WFLT,NBTRI,ExFac,nBMX,CMO)

  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! Compute the two-electron contribution to the core energy
  ETWO = Half*(DDOT_(NBTRI,WDLT,1,WFLT,1)-EONE)
  call mma_deallocate(WDSQ)
  call mma_deallocate(WDLT)
  ! Previous section was bypassed if NFROT == 0.
end if

ECORE = POTNUC+EONE+ETWO
if (IFTEST /= 0) then
  write(u6,'(6X,A,ES20.10)') 'NUCLEAR REPULSION ENERGY:',POTNUC
  write(u6,'(6X,A,ES20.10)') 'ONE-ELECTRON CORE ENERGY:',EONE
  write(u6,'(6X,A,ES20.10)') 'TWO-ELECTRON CORE ENERGY:',ETWO
  write(u6,'(6X,A,ES20.10)') '       TOTAL CORE ENERGY:',ECORE
end if

! Allocate FMO, TMP:
NWTMP = 2*NBMX**2
call mma_allocate(WFMO,notri,LABEL='WFMO')
call mma_allocate(WTMP,NWTMP,LABEL='WTMP')

! Transform one-electron effective Hamiltonian:
WFMO(:) = Zero
WTMP(:) = Zero
ICMO = 1
IAO = 1
IMO = 1
do ISYM=1,NSYM
  ICMO = ICMO+NBAS(ISYM)*NFRO(ISYM)
  IOFF = 1+NBAS(ISYM)*NBAS(ISYM)
  if (NORB(ISYM) > 0) then
    call SQUARE(WFLT(IAO),WTMP,1,NBAS(ISYM),NBAS(ISYM))

    call DGEMM_('T','N',NORB(ISYM),NBAS(ISYM),NBAS(ISYM),One,CMO(ICMO),NBAS(ISYM),WTMP,NBAS(ISYM),Zero,WTMP(IOFF),NORB(ISYM))

    call DGEMM_Tri('N','N',NORB(ISYM),NORB(ISYM),NBAS(ISYM),One,WTMP(IOFF),NORB(ISYM),CMO(ICMO),NBAS(ISYM),Zero,WFMO(IMO), &
                   NORB(ISYM))
  end if
  ICMO = ICMO+NBAS(ISYM)*(NORB(ISYM)+NDEL(ISYM))
  IAO = IAO+NBAS(ISYM)*(NBAS(ISYM)+1)/2
  IMO = IMO+NORB(ISYM)*(NORB(ISYM)+1)/2
end do

if (IFTEST /= 0) then
  write(u6,*) '      EFFECTIVE 1-EL HAMILTONIAN IN MO BASIS'
  ISTLT = 1
  do ISYM=1,NSYM
    if (NORB(ISYM) > 0) then
      write(u6,'(6X,A,I2)') ' SYMMETRY SPECIES:',ISYM
      call TRIPRT(' ',' ',WFMO(ISTLT),NORB(ISYM))
      ISTLT = ISTLT+NORB(ISYM)*(NORB(ISYM)+1)/2
    end if
  end do
end if

HONE(:) = Zero
HONE(1:NoTri) = WFMO(:)

call mma_deallocate(WTMP)
call mma_deallocate(WFMO)
call mma_deallocate(WFLT)

end subroutine TRAONE
