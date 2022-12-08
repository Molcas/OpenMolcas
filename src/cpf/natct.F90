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
! Copyright (C) 1986, Per E. M. Siegbahn                               *
!               1986, Margareta R. A. Blomberg                         *
!***********************************************************************

subroutine NATCT(C,FC)

use cpf_global, only: BNAME, DETOT, ETOT, ICASE, ICPF, INCPF, INDX, ISDCI, ITOC17, JSY, Lu_CPFORB, Lu_TraOne, NBAS, NORB, NPFRO, &
                      NSYM
use OneDat, only: sNoNuc, sNoOri
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, Two
use Definitions, only: wp, iwp, u6

implicit none
real(kind=wp), intent(inout) :: C(*)
real(kind=wp), intent(in) :: FC(*)
integer(kind=iwp) :: ICMO, iComp, IDISK, iDum, iDummy(7,8), iiRC, IOCC, iOpt, iSYM, M, N2SUM, n2Tri, nbMax, NSUM
real(kind=wp) :: dum, Dummy(1), EREL, ErelDC, ErelMV
character(len=72) :: Header
character(len=8) :: Label
real(kind=wp), allocatable :: CAO(:), CMO(:), CMO2(:), D(:), DSYM(:), OCC(:), OP(:), S(:)

NSUM = 0
N2SUM = 0
n2Tri = 0
nbMax = 0
do ISYM=1,NSYM
  nbMax = max(nbMax,nBas(iSym))
  NSUM = NSUM+NBAS(ISYM)
  N2SUM = N2SUM+NBAS(ISYM)**2
  n2Tri = n2Tri+nBas(iSym)*(nBas(iSym)+1)/2
end do

! Read MO coefficients
IDISK = ITOC17(1)
call mma_allocate(CMO,N2SUM,label='CMO')
call dDAFILE(Lu_TraOne,2,CMO,N2SUM,IDISK)

! Loop over irreps and compute natural orbitals

call mma_allocate(OCC,NSUM,label='OCC')
call mma_allocate(CMO2,nbMax**2,label='CMO2')
call mma_allocate(DSYM,nbMax**2,label='DSYM')
call mma_allocate(CAO,nbMax**2,label='CAO')
IOCC = 1
ICMO = 1
do M=1,NSYM
  ! set occupation number of orbitals prefrozen in MOTRA
  OCC(IOCC:IOCC+NBAS(M)-1) = Zero
  ! skip orbitals prefrozen in MOTRA
  OCC(IOCC:IOCC+NPFRO(M)-1) = Two
  call NATORB_CPF(FC,CMO(ICMO+NBAS(M)*NPFRO(M)),CMO2,DSYM,CAO,OCC(IOCC+NPFRO(M)),M)
  call DCOPY_(NORB(M)*NBAS(M),CAO,1,CMO(ICMO+NBAS(M)*NPFRO(M)),1)
  ICMO = ICMO+NBAS(M)**2
  IOCC = IOCC+NBAS(M)
end do
call mma_deallocate(CMO2)
call mma_deallocate(DSYM)
call mma_deallocate(CAO)

call mma_allocate(D,n2Sum,label='D')
call mma_allocate(OP,n2Tri,label='OP')
call RelEne(ErelMV,ErelDC,nSym,nBas,CMO,OCC,D,OP)
call mma_deallocate(D)
call mma_deallocate(OP)

EREL = ERELMV+ERELDC
write(u6,'(/,5X,A)') 'FIRST ORDER RELATIVISTIC CORRECTIONS'
write(u6,'(5X,A,F17.8)') 'MASS-VELOCITY        ',ErelMV
write(u6,'(5X,A,F17.8)') '1-EL DARWIN CONTACT  ',ErelDC
write(u6,'(5X,A,F17.8)') 'TOTAL REL. CORRECTION',Erel
if (ISDCI == 1) then
  write(u6,'(5X,A,F17.8)') 'REL. CI ENERGY       ',ETOT+Erel
  write(u6,'(5X,A,F17.8)') 'REL. CI+Q ENERGY     ',DETOT+Erel
else
  write(u6,'(5X,A,F17.8)') 'TOTAL REL. ENERGY    ',ETOT+Erel
end if
call PRWF_CPF(ICASE,JSY,INDX,C)
if (iCPF == 1) then
  Header = ' CPF natural orbitals'
else if (iSDCI == 1) then
  Header = ' SDCI natural orbitals'
else if (iNCPF == 1) then
  Header = ' ACPF natural orbitals'
else
  Header = ' MCPF natural orbitals'
end if
call Primo(Header,.true.,.false.,1.0e-4_wp,dum,nSym,nBas,nBas,BName,Dummy,OCC,CMO,-1)

! Read the overlap matrix in ao basis
iiRC = -1
iOpt = ibset(ibset(0,sNoOri),sNoNuc)
Label = 'MLTPL  0'
iComp = 1
call mma_allocate(S,n2Tri,label='S')
call RdOne(iiRC,iOpt,Label,iComp,S,iDum)
if (iiRC /= 0) then
  write(u6,*) 'Natct: Error reading overlap matrix!'
  call Abend()
end if
call Charge(nSym,nBas,BName,CMO,OCC,S,2,.true.,.true.)
call Prpt_old(nSym,nBas,nSum,n2Sum,CMO,OCC)
call mma_deallocate(S)

if (iCPF == 1) then
  Header = '* CPF NO COEFS'
else if (iSDCI == 1) then
  Header = '* SDCI NO COEFS'
else if (iNCPF == 1) then
  Header = '* ACPF NO COEFS'
else
  Header = '* MCPF NO COEFS'
end if
call WrVec('CPFORB',Lu_CPFORB,'CO',nSym,nBas,nBas,CMO,OCC,Dummy,iDummy,Header)
call mma_deallocate(CMO)
call mma_deallocate(OCC)

return

end subroutine NATCT
