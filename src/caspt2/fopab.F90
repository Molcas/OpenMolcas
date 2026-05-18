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

subroutine FOPAB(FIFA,NFIFA,IBRA,IKET,FOPEL)

use constants, only: Zero, One, Two
use sguga, only: SGS, L2ACT, EXS, CIS
use caspt2_global, only: LUCIEX, IDCIEX
use stdalloc, only: mma_allocate, mma_deallocate
use caspt2_module, only: NSYM, NORB, NISH, ISCF, NCONF, STSYM, NASH, NAES
use definitions, only: iwp, wp, u6

implicit none
integer(kind=iwp), intent(in) :: NFIFA, IBRA, IKET
real(kind=wp), intent(in) :: FIFA(NFIFA)
real(kind=wp), intent(out) :: FOPEL
integer(kind=iwp) IOFF(8)
integer(kind=iwp) :: nLev
real(kind=wp), allocatable :: BRA(:), KET(:), SGM(:)
integer(kind=iwp) IOF, ISYM, IFTEST, IJ, I, ID, II, ISCR, IST, ISU, IT, ITABS, ITTOT, ITUTOT, IU, IUABS, IUTOT, J, LEVT, LEVU, NI
real(kind=wp) ESUM, OCC, EINACT, FTU, TRC
real(kind=wp), external :: DDot_

nLev = SGS%nLev

! Procedure for computing one matrix element of the Fock matrix in the
! basis of the CASSCF states: <BRA|FOP|KET>
! In: The (possibly average) Fock matrix, active indices only, over the
! original CASSCF orbitals and the indices of the two states

! Offset table for accessing FIFA array:
IOF = 0
do ISYM=1,NSYM
  IOFF(ISYM) = IOF
  IOF = IOF+(NORB(ISYM)*(NORB(ISYM)+1))/2
end do

IFTEST = 0
if (IFTEST > 0) then
  write(u6,*) ' The FIFA array:'
  do ISYM=1,NSYM
    IJ = IOFF(ISYM)+1
    do I=1,NORB(ISYM)
      write(u6,'(1x,5F16.8)') (FIFA(IJ+J),J=0,I-1)
      IJ = IJ+I
    end do
  end do
end if

! Specialized code for Closed-shell or Hi-spin HF:
! Sum up diagonal elements of FIFA times occ. number
! FIXME: This only works for diagonal elements, thus
! in the case of XMS this will not work...
if ((ISCF == 1) .or. (ISCF == 2)) then
  ESUM = Zero
  if (IBRA == IKET) then
    do ISYM=1,NSYM
      OCC = Two
      do I=1,NISH(ISYM)
        ESUM = ESUM+OCC*FIFA(IOFF(ISYM)+(I*(I+1))/2)
      end do
      if (ISCF == 2) OCC = One
      do J=1,NASH(ISYM)
        I = NISH(ISYM)+J
        ESUM = ESUM+OCC*FIFA(IOFF(ISYM)+(I*(I+1))/2)
      end do
    end do
  else
    write(u6,*) ' Warning: neglecting the off-diagonal entries'
    write(u6,*) ' of H0, XMS will be equal to MS!'
  end if
  FOPEL = ESUM
  return
end if

! General CASSCF or RASSCF case:
! Sum up trace of FIFA over inactive orbitals only:
TRC = Zero
do ISYM=1,NSYM
  do I=1,NISH(ISYM)
    II = IOFF(ISYM)+(I*(I+1))/2
    TRC = TRC+FIFA(II)
  end do
end do
! Contribution from inactive orbitals:
EINACT = Two*TRC

if (IFTEST > 0) write(u6,*) ' Energy contrib from inactive orbitals:',EINACT

! Allocate arrays for ket and bra wave functions
call mma_allocate(BRA,NCONF,Label='BRA')
call mma_allocate(KET,NCONF,Label='KET')
! Allocate array for sigma = Fock operator acting on ket:
call mma_allocate(SGM,NCONF,LABEL='SGM')

! Load ket wave function
ID = IDCIEX(IKET)
call DDAFILE(LUCIEX,2,KET,NCONF,ID)

if (IFTEST > 0) then
  write(u6,*) ' IKET:',IKET
  write(u6,*) ' Ket CI array:'
  ISCR = min(NCONF,20)
  write(u6,'(1x,5F16.8)') (KET(I),I=1,ISCR)
end if

! Compute (lowering part of) FIFA operator acting on
! the ket wave function.
call DCOPY_(NCONF,[Zero],0,SGM,1)
do LEVU=1,NLEV
  IUABS = L2ACT(LEVU)
  ISU = SGS%ISM(LEVU)
  IU = IUABS-NAES(ISU)
  NI = NISH(ISU)
  IUTOT = NI+IU
  do LEVT=1,LEVU
    if (SGS%ISM(LEVT) /= ISU) cycle
    ITABS = L2ACT(LEVT)
    IST = ISU
    IT = ITABS-NAES(IST)
    ITTOT = NI+IT
    ITUTOT = (IUTOT*(IUTOT-1))/2+ITTOT
    if (ITTOT > IUTOT) ITUTOT = (ITTOT*(ITTOT-1))/2+IUTOT
    FTU = FIFA(IOFF(ISU)+ITUTOT)
    if (abs(FTU) < 1.0e-16_wp) cycle
    call SG_Epq_Psi(SGS,CIS,EXS,LEVT,LEVU,FTU,STSYM,KET,SGM)
  end do
end do
! Add contribution from inactive part:
call DAXPY_(NCONF,EINACT,KET,1,SGM,1)

if (IFTEST > 0) then
  write(u6,*) ' SGM array from (lowering F)|KET>:'
  write(u6,'(1x,5F16.8)') (SGM(I),I=1,NCONF)
end if

! Load bra wave function
ID = IDCIEX(IBRA)
call DDAFILE(LUCIEX,2,BRA,NCONF,ID)

! Put matrix element into FOPEL:
FOPEL = DDOT_(NCONF,BRA,1,SGM,1)

if (IFTEST > 0) then
  write(u6,*) ' FOPEL is now:'
  write(u6,'(1x,5f16.8)') FOPEL
end if

! Compute (strictly lowering part of) FIFA operator acting on |BRA>.
! We are computing contributions <KET|Etu|BRA> with t<u, then
! using them as <BRA|Eut|KET>
! Note that I already have BRA in memory
call DCOPY_(NCONF,[Zero],0,SGM,1)
do LEVU=2,NLEV
  IUABS = L2ACT(LEVU)
  ISU = SGS%ISM(LEVU)
  IU = IUABS-NAES(ISU)
  NI = NISH(ISU)
  IUTOT = NI+IU
  do LEVT=1,LEVU-1
    if (SGS%ISM(LEVT) /= ISU) cycle
    ITABS = L2ACT(LEVT)
    IST = ISU
    IT = ITABS-NAES(IST)
    ITTOT = NI+IT
    ITUTOT = (IUTOT*(IUTOT-1))/2+ITTOT
    if (ITTOT > IUTOT) ITUTOT = (ITTOT*(ITTOT-1))/2+IUTOT
    FTU = FIFA(IOFF(ISU)+ITUTOT)
    if (abs(FTU) < 1.0e-16_wp) cycle
    call SG_Epq_Psi(SGS,CIS,EXS,LEVT,LEVU,FTU,STSYM,BRA,SGM)
  end do
end do

if (IFTEST > 0) then
  write(u6,*) ' SGM array from (strictly lowering F)|BRA>:'
  write(u6,'(1x,5F16.8)') (SGM(I),I=1,NCONF)
end if

! Load ket wave function
ID = IDCIEX(IKET)
call DDAFILE(LUCIEX,2,KET,NCONF,ID)

! Add contribution to matrix element FOPEL
FOPEL = FOPEL+DDOT_(NCONF,KET,1,SGM,1)

if (IFTEST > 0) then
  write(u6,*) ' FOPEL is now:'
  write(u6,'(1x,5f16.8)') FOPEL
end if

call mma_deallocate(SGM)
call mma_deallocate(BRA)
call mma_deallocate(KET)

end subroutine FOPAB

