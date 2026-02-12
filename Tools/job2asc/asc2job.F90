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

program ASC2JOB

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
call PrgmInit('Asc2Job')

call f_inquire('FMTIPH',Found)
if (.not. Found) then
  call WarningMessage(2,'FMTIPH not found')
  stop
end if

JOBIPH = isFreeUnit(15)
call DANAME(JOBIPH,'JOBIPH')

FMTIPH = isFreeUnit(JOBIPH+1)
call MOLCAS_OPEN(FMTIPH,'FMTIPH')

read(FMTIPH,*)
read(FMTIPH,'(15I10)') IADR15(1:30)
IAD15 = 0
call IDAFILE(JOBIPH,1,IADR15,30,IAD15)

nName = (LenIn+8)*mxOrb
nHeader = 144
nTitle = 4*18*mxTit

read(FMTIPH,*)
read(FMTIPH,100) nActEl
read(FMTIPH,*)
read(FMTIPH,100) iSpin
read(FMTIPH,*)
read(FMTIPH,100) nSym
read(FMTIPH,*)
read(FMTIPH,100) lSym
read(FMTIPH,*)
read(FMTIPH,200) nFro(1:MxSym)
read(FMTIPH,*)
read(FMTIPH,200) nISh(1:MxSym)
read(FMTIPH,*)
read(FMTIPH,200) nASh(1:MxSym)
read(FMTIPH,*)
read(FMTIPH,200) nDel(1:MxSym)
read(FMTIPH,*)
read(FMTIPH,200) nBas(1:MxSym)
read(FMTIPH,*)
read(FMTIPH,100) nConf
read(FMTIPH,*)
read(FMTIPH,'(36A)') Header(1:72)
read(FMTIPH,*)
read(FMTIPH,'(80A)') Title(1:9)
read(FMTIPH,*)
read(FMTIPH,400) PotNuc
read(FMTIPH,*)
read(FMTIPH,100) lRoots
read(FMTIPH,*)
read(FMTIPH,100) nRoots
read(FMTIPH,*)
read(FMTIPH,200) iRoot(1:MxRoot)
read(FMTIPH,*)
read(FMTIPH,200) nRS1(1:MxSym)
read(FMTIPH,*)
read(FMTIPH,200) nRS2(1:MxSym)
read(FMTIPH,*)
read(FMTIPH,200) nRS3(1:MxSym)
read(FMTIPH,*)
read(FMTIPH,100) nHole1
read(FMTIPH,*)
read(FMTIPH,100) nElec3
read(FMTIPH,*)
read(FMTIPH,100) iPT2
read(FMTIPH,*)
read(FMTIPH,500) Weight(1:MxRoot)

IAD15 = IADR15(1)
call WR_RASSCF_Info(JOBIPH,1,IAD15,NACTEL,ISPIN,NSYM,LSYM,NFRO,NISH,NASH,NDEL,NBAS,MxSym,BName,nName,NCONF,HEADER,nHeader,TITLE, &
                    nTitle,POTNUC,LROOTS,NROOTS,IROOT,MxRoot,NRS1,NRS2,NRS3,NHOLE1,NELEC3,IPT2,WEIGHT)

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

read(FMTIPH,*)
read(FMTIPH,300) BName(1:NTOT)

IAD15 = IADR15(2)
call mma_allocate(ADR1,NTOT2,Label='ADR1')
call mma_allocate(ADR2,NTOT,Label='ADR2')
read(FMTIPH,*)
read(FMTIPH,*)
read(FMTIPH,*)
read(FMTIPH,500) ADR1(:)
read(FMTIPH,*)
read(FMTIPH,*)
read(FMTIPH,500) ADR2(:)
call DDAFILE(JOBIPH,1,ADR1,NTOT2,IAD15)
call DDAFILE(JOBIPH,1,ADR2,NTOT,IAD15)
call mma_deallocate(ADR1)
call mma_deallocate(ADR2)

IAD15 = IADR15(3)
call mma_allocate(ADR1,NACPAR,Label='ADR1')
call mma_allocate(ADR2,NACPR2,Label='ADR2')
read(FMTIPH,*)
do I=1,LROOTS
  read(FMTIPH,*)
  read(FMTIPH,*)
  read(FMTIPH,*)
  read(FMTIPH,500) ADR1(:)
  call DDafIle(JOBIPH,1,ADR1,NACPAR,IAD15)
  read(FMTIPH,*)
  read(FMTIPH,*)
  read(FMTIPH,500) ADR1(:)
  call DDafIle(JOBIPH,1,ADR1,NACPAR,IAD15)
  read(FMTIPH,*)
  read(FMTIPH,*)
  read(FMTIPH,500) ADR2(:)
  call DDafIle(JOBIPH,1,ADR2,NACPR2,IAD15)
  read(FMTIPH,*)
  read(FMTIPH,*)
  read(FMTIPH,500) ADR2(:)
  call DDafIle(JOBIPH,1,ADR2,NACPR2,IAD15)
end do
call mma_deallocate(ADR1)
call mma_deallocate(ADR2)

IAD15 = IADR15(4)
call mma_allocate(ADR,NCONF,Label='ADR')
read(FMTIPH,*)
do I=1,LROOTS
  read(FMTIPH,*)
  read(FMTIPH,*)
  read(FMTIPH,500) ADR(:)
  call DDafIle(JOBIPH,1,ADR,NCONF,IAD15)
end do
call mma_deallocate(ADR)

IAD15 = IADR15(5)
call mma_allocate(ADR,NFOCK,Label='ADR')
read(FMTIPH,*)
read(FMTIPH,*)
read(FMTIPH,500) ADR(:)
call DDAFILE(JOBIPH,1,ADR,NFOCK,IAD15)
call mma_deallocate(ADR)

IAD15 = IADR15(6)
call mma_allocate(ADR,mxRoot*mxIter,Label='ADR')
read(FMTIPH,*)
read(FMTIPH,*)
read(FMTIPH,500) ADR(:)
call DDAFILE(JOBIPH,1,ADR,mxRoot*mxIter,IAD15)
call mma_deallocate(ADR)

IAD15 = IADR15(7)
call mma_allocate(ADR,6*mxIter,Label='ADR')
read(FMTIPH,*)
read(FMTIPH,*)
read(FMTIPH,510) ADR(:)
call DDAFILE(JOBIPH,1,ADR,6*mxIter,IAD15)
call mma_deallocate(ADR)

IAD15 = IADR15(9)
call mma_allocate(ADR,NTOT2,Label='ADR')
read(FMTIPH,*)
read(FMTIPH,*)
read(FMTIPH,500) ADR(:)
call DDAFILE(JOBIPH,1,ADR,NTOT2,IAD15)
call mma_deallocate(ADR)

IAD15 = IADR15(10)
call mma_allocate(ADR,NTOT3,Label='ADR')
read(FMTIPH,*)
read(FMTIPH,*)
read(FMTIPH,500) ADR(:)
call DDAFILE(JOBIPH,1,ADR,NTOT3,IAD15)
read(FMTIPH,*)
read(FMTIPH,*)
read(FMTIPH,500) ADR(:)
call DDAFILE(JOBIPH,1,ADR,NTOT3,IAD15)
call mma_deallocate(ADR)

IAD15 = IADR15(11)
call mma_allocate(ADR,NORBT,Label='ADR')
read(FMTIPH,*)
read(FMTIPH,*)
read(FMTIPH,500) ADR(:)
call DDAFILE(JOBIPH,1,ADR,NORBT,IAD15)
call mma_deallocate(ADR)

IAD15 = IADR15(12)
call mma_allocate(ADR1,NTOT2,Label='ADR1')
call mma_allocate(ADR2,NTOT,Label='ADR2')
read(FMTIPH,*)
do I=1,LROOTS
  read(FMTIPH,*)
  read(FMTIPH,*)
  read(FMTIPH,*)
  read(FMTIPH,500) ADR1(:)
  call DDAFILE(JOBIPH,1,ADR1,NTOT2,IAD15)
  read(FMTIPH,*)
  read(FMTIPH,*)
  read(FMTIPH,500) ADR2(:)
  call DDAFILE(JOBIPH,1,ADR2,NTOT,IAD15)
end do
call mma_deallocate(ADR1)
call mma_deallocate(ADR2)

IAD15 = IADR15(14)
call mma_allocate(ADR1,NTOT2,Label='ADR1')
call mma_allocate(ADR2,NTOT,Label='ADR2')
read(FMTIPH,*)
do I=1,LROOTS
  read(FMTIPH,*)
  read(FMTIPH,*)
  read(FMTIPH,*)
  read(FMTIPH,500) ADR1(:)
  call DDAFILE(JOBIPH,1,ADR1,NTOT2,IAD15)
  read(FMTIPH,*)
  read(FMTIPH,*)
  read(FMTIPH,500) ADR2(:)
  call DDAFILE(JOBIPH,1,ADR2,NTOT,IAD15)
end do
call mma_deallocate(ADR1)
call mma_deallocate(ADR2)

IAD15 = IADR15(17)
call mma_allocate(ADR,LROOTS**2,Label='ADR')
read(FMTIPH,*)
read(FMTIPH,*)
read(FMTIPH,500) ADR(:)
call DDAFILE(JOBIPH,1,ADR,LROOTS**2,IAD15)
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

end program ASC2JOB
