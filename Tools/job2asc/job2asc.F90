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

program JOB2ASC
! SVC 20071004
! convert JOBIPH to a formatted file FMTIPH
! with some additional explanation.

use rasscf_global, only: BName, Header, IADR15, iPT2, iRoot, lRoots, NACPAR, NACPR2, NORBT, nRoots, NTOT3, PotNuc, Title, Weight
use general_data, only: ispin, jobiph, nactel, nash, nbas, nconf, ndel, ndel, nelec3, nfro, nhole1, nish, norb, nrs1, nrs2, nrs3, &
                        nsym, ntot, ntot2
use Molcas, only: LenIn, MxOrb, MxRoot, MxSym
use RASDim, only: MxIter, MxTit
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: FMTIPH, I, IAD15, ISYM, LSYM, NASHT, NFOCK, nHeader, nName, nTitle
logical(kind=iwp) :: Found
real(kind=wp), allocatable :: ADR1(:), ADR2(:), ADR(:)
integer(kind=iwp), external :: isFreeUnit

call IniMem()
call Init_LinAlg()
call PrgmInit('Job2Asc')

call f_inquire('JOBIPH',Found)
if (.not. Found) then
  call WarningMessage(2,'JOBIPH not found')
  stop
end if

JOBIPH = isFreeUnit(15)
call DANAME(JOBIPH,'JOBIPH')

FMTIPH = isFreeUnit(JOBIPH+1)
call MOLCAS_OPEN(FMTIPH,'FMTIPH')

IAD15 = 0
call IDAFILE(JOBIPH,2,IADR15,30,IAD15)
write(FMTIPH,*) 'print out of IADR15: '
write(FMTIPH,'(15I10)') IADR15(1:30)

nName = (LenIn+8)*mxOrb
nHeader = 144
nTitle = 4*18*mxTit

IAD15 = IADR15(1)
call WR_RASSCF_Info(JOBIPH,2,IAD15,NACTEL,ISPIN,NSYM,LSYM,NFRO,NISH,NASH,NDEL,NBAS,MxSym,BName,nName,NCONF,HEADER,nHeader,TITLE, &
                    nTitle,POTNUC,LROOTS,NROOTS,IROOT,MxRoot,NRS1,NRS2,NRS3,NHOLE1,NELEC3,IPT2,WEIGHT)

write(FMTIPH,*) 'number of active electrons'
write(FMTIPH,100) nActEl
write(FMTIPH,*) 'spin of the wave function'
write(FMTIPH,100) iSpin
write(FMTIPH,*) 'number of irreps'
write(FMTIPH,100) nSym
write(FMTIPH,*) 'irrep of the wave function'
write(FMTIPH,100) lSym
write(FMTIPH,*) 'number of frozen orbitals'
write(FMTIPH,200) nFro(1:MxSym)
write(FMTIPH,*) 'number of inactive orbitals'
write(FMTIPH,200) nISh(1:MxSym)
write(FMTIPH,*) 'number of active orbitals'
write(FMTIPH,200) nASh(1:MxSym)
write(FMTIPH,*) 'number of deleted orbitals'
write(FMTIPH,200) nDel(1:MxSym)
write(FMTIPH,*) 'number of basis functions'
write(FMTIPH,200) nBas(1:MxSym)
write(FMTIPH,*) 'number of configurations'
write(FMTIPH,100) nConf
write(FMTIPH,*) 'Header'
write(FMTIPH,'(36A)') Header(1:72)
write(FMTIPH,*) 'Title'
write(FMTIPH,'(80A)') Title(1:9)
write(FMTIPH,*) 'nuclear potential energy'
write(FMTIPH,400) PotNuc
write(FMTIPH,*) 'number of roots in the small CI'
write(FMTIPH,100) lRoots
write(FMTIPH,*) 'number of roots in the averaging'
write(FMTIPH,100) nRoots
write(FMTIPH,*) 'rootnumbers'
write(FMTIPH,200) iRoot(1:MxRoot)
write(FMTIPH,*) 'number of RAS1 orbitals'
write(FMTIPH,200) nRS1(1:MxSym)
write(FMTIPH,*) 'number of RAS2 orbitals'
write(FMTIPH,200) nRS2(1:MxSym)
write(FMTIPH,*) 'number of RAS3 orbitals'
write(FMTIPH,200) nRS3(1:MxSym)
write(FMTIPH,*) 'max number of holes in RAS1'
write(FMTIPH,100) nHole1
write(FMTIPH,*) 'max number of electrons in RAS3'
write(FMTIPH,100) nElec3
write(FMTIPH,*) 'iPT2'
write(FMTIPH,100) iPT2
write(FMTIPH,*) 'weight of the roots in the averaging'
write(FMTIPH,500) Weight(1:MxRoot)

NTOT = 0
NTOT2 = 0
NASHT = 0
NTOT3 = 0
NFOCK = 0
do ISYM=1,NSYM
  NTOT = NTOT+nBas(ISYM)
  NTOT2 = NTOT2+nBas(ISYM)**2
  NASHT = NASHT+nASh(ISYM)
  nOrb(iSym) = nBas(iSym)-nFro(iSym)-nDel(iSym)
  NTOT3 = NTOT3+nOrb(iSym)*(nOrb(iSym)+1)/2
  NFOCK = NFOCK+(nISh(ISYM)+nASh(ISYM))**2
end do
NACPAR = NASHT*(NASHT+1)/2
NACPR2 = NACPAR*(NACPAR+1)/2

write(FMTIPH,*) 'basis function labels and type'
write(FMTIPH,300) BName(1:NTOT)

IAD15 = IADR15(2)
call mma_allocate(ADR1,NTOT2,Label='ADR1')
call mma_allocate(ADR2,NTOT,Label='ADR2')
call DDAFILE(JOBIPH,2,ADR1,NTOT2,IAD15)
call DDAFILE(JOBIPH,2,ADR2,NTOT,IAD15)
write(FMTIPH,*) 'average orbitals'
write(FMTIPH,*) 'MO coefficients'
write(FMTIPH,*) 'length = ',NTOT2
write(FMTIPH,500) ADR1(:)
write(FMTIPH,*) 'occupation numbers'
write(FMTIPH,*) 'length = ',NTOT
write(FMTIPH,500) ADR2(:)
call mma_deallocate(ADR1)
call mma_deallocate(ADR2)

IAD15 = IADR15(3)
call mma_allocate(ADR1,NACPAR,Label='ADR1')
call mma_allocate(ADR2,NACPR2,Label='ADR2')
write(FMTIPH,*) 'density matrices for the active orbitals'
do I=1,LROOTS
  write(FMTIPH,*) 'root = ',I
  write(FMTIPH,*) 'D  : one-body density matrix'
  write(FMTIPH,*) 'length = ',NACPAR
  call DDafIle(JOBIPH,2,ADR1,NACPAR,IAD15)
  write(FMTIPH,500) ADR1(:)
  write(FMTIPH,*) 'DS : spin density matrix'
  write(FMTIPH,*) 'length = ',NACPAR
  call DDafIle(JOBIPH,2,ADR1,NACPAR,IAD15)
  write(FMTIPH,500) ADR1(:)
  write(FMTIPH,*) 'P  : symmetric two-body density matrix'
  write(FMTIPH,*) 'length = ',NACPR2
  call DDafIle(JOBIPH,2,ADR2,NACPR2,IAD15)
  write(FMTIPH,500) ADR2(:)
  write(FMTIPH,*) 'PA : antisymmetric two-body density matrix'
  write(FMTIPH,*) 'length = ',NACPR2
  call DDafIle(JOBIPH,2,ADR2,NACPR2,IAD15)
  write(FMTIPH,500) ADR2(:)
end do
call mma_deallocate(ADR1)
call mma_deallocate(ADR2)

IAD15 = IADR15(4)
call mma_allocate(ADR,NCONF,Label='ADR')
write(FMTIPH,*) 'CI coefficients'
do I=1,LROOTS
  write(FMTIPH,*) 'root = ',I
  write(FMTIPH,*) 'length = ',NCONF
  call DDafIle(JOBIPH,2,ADR,NCONF,IAD15)
  write(FMTIPH,500) ADR(:)
end do
call mma_deallocate(ADR)

IAD15 = IADR15(5)
call mma_allocate(ADR,NFOCK,Label='ADR')
call DDAFILE(JOBIPH,2,ADR,NFOCK,IAD15)
write(FMTIPH,*) 'Fock matrix for the occupied orbitals'
write(FMTIPH,*) 'length = ',NFOCK
write(FMTIPH,500) ADR(:)
call mma_deallocate(ADR)

IAD15 = IADR15(6)
call mma_allocate(ADR,mxRoot*mxIter,Label='ADR')
call DDAFILE(JOBIPH,2,ADR,mxRoot*mxIter,IAD15)
write(FMTIPH,*) 'energies from array ENER(mxRoot,mxIter)'
write(FMTIPH,*) 'length = ',mxRoot*mxIter
write(FMTIPH,500) ADR(:)
call mma_deallocate(ADR)

IAD15 = IADR15(7)
call mma_allocate(ADR,6*mxIter,Label='ADR')
call DDAFILE(JOBIPH,2,ADR,6*mxIter,IAD15)
write(FMTIPH,*) 'convergence parameters from array CONV(6,mxIter)'
write(FMTIPH,*) 'length = ',6*mxIter
write(FMTIPH,510) ADR(:)
call mma_deallocate(ADR)

IAD15 = IADR15(9)
call mma_allocate(ADR,NTOT2,Label='ADR')
call DDAFILE(JOBIPH,2,ADR,NTOT2,IAD15)
write(FMTIPH,*) 'final canonical MOs written in FCKPT2'
write(FMTIPH,*) 'length = ',NTOT2
write(FMTIPH,500) ADR(:)
call mma_deallocate(ADR)

IAD15 = IADR15(10)
call mma_allocate(ADR,NTOT3,Label='ADR')
call DDAFILE(JOBIPH,2,ADR,NTOT3,IAD15)
write(FMTIPH,*) 'inactive Fock matrix FI for CASPT2 in MO basis'
write(FMTIPH,*) 'length = ',NTOT3
write(FMTIPH,500) ADR(:)
call DDAFILE(JOBIPH,2,ADR,NTOT3,IAD15)
write(FMTIPH,*) 'CASPT2 Fock matrix FP'
write(FMTIPH,*) 'length = ',NTOT3
write(FMTIPH,500) ADR(:)
call mma_deallocate(ADR)

IAD15 = IADR15(11)
call mma_allocate(ADR,NORBT,Label='ADR')
call DDAFILE(JOBIPH,2,ADR,NORBT,IAD15)
write(FMTIPH,*) 'diagonal of the CASPT2 Fock matrix FP'
write(FMTIPH,*) 'length = ',NORBT
write(FMTIPH,500) ADR(:)
call mma_deallocate(ADR)

IAD15 = IADR15(12)
call mma_allocate(ADR1,NTOT2,Label='ADR1')
call mma_allocate(ADR2,NTOT,Label='ADR2')
write(FMTIPH,*) 'natural orbitals for the final wave functions'
do I=1,LROOTS
  write(FMTIPH,*) 'root = ',I
  write(FMTIPH,*) 'MO coefficients'
  write(FMTIPH,*) 'length = ',NTOT2
  call DDAFILE(JOBIPH,2,ADR1,NTOT2,IAD15)
  write(FMTIPH,500) ADR1(:)
  write(FMTIPH,*) 'occupation numbers'
  write(FMTIPH,*) 'length = ',NTOT
  call DDAFILE(JOBIPH,2,ADR2,NTOT,IAD15)
  write(FMTIPH,500) ADR2(:)
end do
call mma_deallocate(ADR1)
call mma_deallocate(ADR2)

IAD15 = IADR15(14)
call mma_allocate(ADR1,NTOT2,Label='ADR1')
call mma_allocate(ADR2,NTOT,Label='ADR2')
write(FMTIPH,*) 'spin orbitals for the final wave functions'
do I=1,LROOTS
  write(FMTIPH,*) 'root = ',I
  write(FMTIPH,*) 'MO coefficients'
  write(FMTIPH,*) 'length = ',NTOT2
  call DDAFILE(JOBIPH,2,ADR1,NTOT2,IAD15)
  write(FMTIPH,500) ADR1(:)
  write(FMTIPH,*) 'occupation numbers'
  write(FMTIPH,*) 'length = ',NTOT
  call DDAFILE(JOBIPH,2,ADR2,NTOT,IAD15)
  write(FMTIPH,500) ADR2(:)
end do
call mma_deallocate(ADR1)
call mma_deallocate(ADR2)

IAD15 = IADR15(17)
call mma_allocate(ADR,LROOTS**2,Label='ADR')
write(FMTIPH,*) 'effective hamiltonian from MS-CASPT2'
write(FMTIPH,*) 'length = ',LROOTS**2
call DDAFILE(JOBIPH,2,ADR,LROOTS**2,IAD15)
write(FMTIPH,500) ADR(:)
call mma_deallocate(ADR)

call DACLOS(JOBIPH)
close(FMTIPH)

call GetMem('Finish','LIST','REAL',I,0)
call GetMem('Finish','TERM','REAL',I,0)

100 format(i10)
200 format(2x,10i10)
300 format(a)
400 format(es20.12)
500 format(2x,5es20.12)
510 format(2x,6es20.12)

end program JOB2ASC
