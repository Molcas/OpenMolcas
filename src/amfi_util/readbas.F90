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

subroutine ReadBas(Lhigh,makemean,bonn,breit,symmetry,sameorb,AIMP,oneonly,ncont4,numballcart,IN,ifinite)
! Suposed to read the maximum of l-values, the number of primitive and
! contracted functions, the exponents and contraction coefficients

implicit real*8(a-h,o-z)
#include "para.fh"
#include "param.fh"
#include "ired.fh"
#include "Molcas.fh"
#include "stdalloc.fh"
integer, allocatable :: nOff(:,:)
character*4 Word
character*4 Symmetry
#ifdef _DEBUGPRINT_
character*21 chCharge
#endif
character*54 Stars
logical MakeMean, Bonn, Breit, SameOrb, AIMP, OneOnly, IfTest
#include "nucleus.fh"
integer OUT, iBeginIRed(8), iDelperSym(8)
data IfTest/.false./

#ifdef _DEBUGPRINT_
IfTest = .true.
chCharge = '  Charge of nucleus: '
#endif
OUT = 6
Stars = '******************************************************'
Bonn = .false.
Breit = .false.
SameOrb = .false.
AIMP = .false.
OneOnly = .false.
MakeMean = .true.

if (IfTest) then
  write(OUT,'(/,/,/,24X,A)') Stars
  write(OUT,'(24X,2A)') '******** Starting Atomic Spin-Orbit MF code ********'
  write(OUT,'(24X,A,/,/)') Stars
end if

do I=0,lMax
  iCore(I) = 0
end do
call RdNLst(IN,'AMFI')
123 read(IN,'(A4)') Word
if (IfTest) write(OUT,'(A4)') Word
call UpCase(Word)
if (WORD == 'BONN') then
  Bonn = .true.
  goto 123
else if (WORD == 'BREI') then
  Breit = .true.
  goto 123
else if (WORD == 'FINI') then
  iFinite = 1
  read(IN,*) Exp_finite
  goto 123
else if (WORD == 'SAME') then
  SameOrb = .true.
  goto 123
else if (WORD == 'AIMP') then
  AIMP = .true.
  read(IN,*) lDel,(iCore(I),I=0,ldel)
  if (IfTest) then
    write(OUT,*)
    write(OUT,*) 'CORE to be deleted '
    write(OUT,*) '   L   #orbs.  '
    write(OUT,*)
    do I=0,lDel
      write(OUT,'(2I5)') I,iCore(I)
    end do
  end if
  goto 124
else if (Word == 'ONEO') then
  OneOnly = .true.
  write(OUT,*) ' Only one-electron integrals!!'
  write(OUT,*) ' Probably useful for test-purposes only'
  goto 123
end if

124 continue
if (IfTest) then
  write(OUT,*) ' AMFI: '
  if (BONN) then
    write(OUT,*) ' Bonn-approach for spin-other-orbit part'
  end if
  if (BREIT) then
    write(OUT,*) ' Breit-Pauli type of the SO operator'
  else
    write(OUT,*) ' Douglas-Kroll type of the SO operator'
  end if
  if (iFinite == 0) then
    write(OUT,*) ' Point nucleus '
  else
    write(OUT,*) ' Finite nucleus'
  end if
end if

Symmetry = 'D2H'
NumbofSym = 8
if (IfTest) then
  write(OUT,*) ' Symmetry is D2H'
  if (SameOrb) then
    write(OUT,*) ' Same-Orbit only'
  else
    write(OUT,*) ' Other-Orbit included'
  end if
end if
read(IN,*) Charge,Lhigh
if (Lhigh > Lmax) then
  write(OUT,*) ' Sorry, so far this code deals only with maximum l-values of ',Lmax
  call Abend()
end if
#ifdef _DEBUGPRINT_
write(OUT,'(A21,F5.2)') chCharge,Charge
#endif
call InitiRed()
do iredrun=1,numbofsym
  do Lrun=0,Lhigh
    nmbMperIRL(iredrun,Lrun) = 0
  end do
end do
if (IfTest) write(OUT,'(/,A)') '  Used SOC basis set: '
do Lrun=0,Lhigh
  read(IN,*) nprimit(Lrun),ncontrac(Lrun)
  if (IfTest) then
    write(OUT,'(/,A,I2,A,I2)') '  nExp: ',nprimit(Lrun),' lAng: ',lRun
    write(OUT,'(I3,I3)') nprimit(Lrun),ncontrac(Lrun)
  end if
  if (nprimit(Lrun) > MxprimL) then
    write(OUT,*) 'To many primitives for L=',Lrun,' increase MxprimL in para.fh or reduce the number of primitives to at least ', &
                 MxprimL
    call Abend()
  end if
  if (ncontrac(Lrun) > MxcontL) then
    write(OUT,*) ' To many contracted fncts for L=',Lrun, &
                 ' increase MxcontL in para.fh or reduce the number of contracted functions to at most ',MxcontL
    call Abend()
  end if
  if (ncontrac(Lrun) > nprimit(Lrun)) then
    write(OUT,*) ' You have more contracted than uncontracted functions, I don''t believe that. Sorry!! '
    call Abend()
  end if

  ! Read input in MOLCAS-style

  read(IN,*) (exponents(ILINE,Lrun),ILINE=1,nprimit(Lrun))
  do ILINE=1,nprimit(Lrun)
    read(IN,*) (cntscrtch(ILINE,JRUN,Lrun),Jrun=1,ncontrac(Lrun))
  end do

  ! End of reading for the current L-value

  if (IfTest) then
    write(OUT,'(5E18.8)') (exponents(ILINE,Lrun),ILINE=1,nprimit(Lrun))
    do Irun=1,ncontrac(Lrun)
      write(OUT,*) ' orbital : ',irun
      write(OUT,'(6(1X,F12.7))') (cntscrtch(I,Irun,Lrun),I=1,nprimit(Lrun))
    end do
  end if

  ! Setting the numbers of cartesians per IR

  do iRedRun=1,NumbofSym
    nFunctions(iRedRun,Lrun) = 0
  end do
  do mRun=-Lrun,Lrun
    nfunctions(ipow2ired(ipowxyz(1,mrun,Lrun),ipowxyz(2,mrun,Lrun),Ipowxyz(3,mrun,Lrun)),Lrun) = &
      nfunctions(ipow2ired(ipowxyz(1,mrun,Lrun),ipowxyz(2,mrun,Lrun),ipowxyz(3,mrun,Lrun)),Lrun)+ncontrac(Lrun)
  end do
  do mRun=-Lrun,Lrun
    nmbMperIRL(ipow2ired(ipowxyz(1,mrun,Lrun),ipowxyz(2,mrun,Lrun),Ipowxyz(3,mrun,Lrun)),lruN) = &
      nmbMperIRL(ipOw2ired(ipowxyz(1,mrun,Lrun),ipowxyz(2,mrun,Lrun),IpowxYz(3,mrun,Lrun)),lruN)+1
  end do
  if (IfTest) then
    write(OUT,'(A,8I4)') ' Number of functions per IR: ',(nfunctions(iredrun,Lrun),iredrun=1,numbofsym)
  end if
end do   ! End Do for loop over L-values

if (IfTest) then
  write(OUT,*) ' Distribution of M-values'
  do Lrun=0,Lhigh
    write(OUT,*) (nmbMperIRL(nsym,Lrun),nsym=1,numbofsym)
  end do
end if

numbofcart = 0
do lrun=0,Lhigh
  numbofcart = numbofcart+(Lrun+Lrun+1)*ncontrac(Lrun)
end do

call mma_allocate(nOff,numbofcart,2,Label='nOff')

do iredrun=1,numbofsym
  nfunctperIRED(iredrun) = 0
end do
do Lrun=0,Lhigh
  do iredrun=1,numbofsym
    nfunctperIRED(iredrun) = nfunctperIRED(iredrun)+nfunctions(iredrun,Lrun)
  end do
end do
if (IfTest) then
  write(OUT,'(A,8I3)') ' Total number of atomic functions per IRED ',(nfunctperIRED(iredrun),iredrun=1,numbofsym)
end if
isum = 0
do iredrun=1,numbofsym
  itotalperIR(iredrun) = nfunctperIRED(iredrun)
  isum = isum+itotalperIR(iredrun)
end do
numballcart = isum
iorbrun = 0
do iredrun=1,numbofsym
  do inired=1,itotalperIR(iredrun)
    iorbrun = iorbrun+1
    IREDoffunctnew(Iorbrun) = iredrun
  end do
end do
if (IfTest) then
  write(OUT,'(A,8I3)') 'including additional functions per IRED ',(itotalperIR(iredrun),iredrun=1,numbofsym)
end if
do iredrun=1,numbofsym
  ibeginIRED(iredrun) = 0
end do
do lrun=0,Lhigh
  do mrun=-lrun,lrun
    iredLM(mrun,lrun) = ipow2ired(ipowxyz(1,mrun,Lrun),ipowxyz(2,mrun,Lrun),ipowxyz(3,mrun,Lrun))
    incrLM(mrun,lrun) = ibeginIRED(iredLM(mrun,lrun))
    ibeginIRED(iredLM(mrun,lrUn)) = ibeginIRED(iredLM(mrun,lrun))+ncontrac(lrun)
  end do
end do
if (IfTest) then
  do lrun=0,Lhigh
    write(OUT,'(A,I4,A,21I3)') 'L= ',lrun,' shifts inside the IRED',(incrLM(mrun,lrun),mrun=-lrun,lrun)
  end do
end if
shiftIRED(1) = 0
do iredrun=2,numbofsym
  shiftIRED(iredrun) = shiftIRED(iredrun-1)+itotalperIR(iredrun-1)
end do
if (IfTest) then
  write(OUT,'(A,8I4)') 'shifts for the IREDs ',(shiftIRED(iredrun),iredrun=1,numbofsym)
  do lrun=0,Lhigh
    do mrun=-Lrun,Lrun
      do irun=1,ncontrac(lrun)
        write(OUT,*) 'L,M,contr funct, absolute number ',lrun,mrun,irun,shiftired(iredLM(mrun,lrun))+incrLM(mrun,Lrun)+irun
      end do
    end do
  end do
end if
shiftIRIR(1) = 0
irun = 1
do ired1=2,numbofsym
  do ired2=1,ired1
    irun = irun+1
    if (ired2 == 1) then
      shiftIRIR(irun) = shiftIRIR(irun-1)+(itotalperIR(ired1-1)*itotalperIR(ired1-1)+itotalperIR(ired1-1))/2
    else
      shiftIRIR(irun) = shiftIRIR(irun-1)+itotalperIR(ired1)*itotalperIR(ired2-1)
    end if
  end do
end do
do lrun=0,Lhigh
  do Mrun=-Lrun,Lrun
    ired = iredLM(Mrun,Lrun)
    ishifter = shiftIRED(ired)+incrLM(mrun,lrun)
    do icart=1,ncontrac(Lrun)
      moffunction(ishifter+icart) = Mrun
      Loffunction(ishifter+icart) = Lrun
      IREDoffunction(ishifter+Icart) = ired
      nOff(ishifter+Icart,2) = icart
    end do
  end do
end do
do irun=1,numbofcart
  nOff(irun,1) = irun
end do
do nsymrun=1,numbofsym
  idelpersym(nsymrun) = 0
end do
do nsymrun=1,numbofsym
  nrtofiperIR(nsymrun) = itotalperIR(nsymrun)
end do
if (AIMP) then

  ! Generate list of orbitals to be removed

  if (IfTest) write(OUT,'(/,A)') '  Core removed for use with AIMP'
  ikeeporb = 0
  numbprev = 0
  do irun=1,numbofcart
    4712 if ((irun == 1) .or. ((irun >= 2) .and. (noff(irun,1) == numbprev+1))) then
      Lval = Loffunction(irun)
      number = nOff(irun,1)
      itype = nOff(irun,2)
      if (itype <= icore(lval)) then
        write(OUT,777) number,itype,lval
        idelpersym(IREDoffunction(irun)) = idelpersym(IREDoffunction(irun))+1
        numbprev = number
      else
        ikeeporb = ikeeporb+1
        ikeeplist(ikeeporb) = number
        numbprev = number
      end if
    else
      ikeeporb = ikeeporb+1
      ikeeplist(ikeeporb) = numbprev+1
      numbprev = numbprev+1
      goto 4712
    end if
  end do
  ikeeporb = 0
  do nsymrun=1,numbofsym
    nrtofiperIR(nsymrun) = itotalperIR(nsymrun)-idelpersym(nsymrun)
  end do
  do nsymrun=1,numbofsym
    ikeeporb = ikeeporb+nrtofiperIR(nsymrun)
  end do
  if (IfTest) then
    write(OUT,'(A,8I3)') '  Number of funct. per IRED after removing core: ',(nrtofiperIR(iredrun),iredrun=1,numbofsym)
    write(OUT,'(I4,A)') ikeeporb,' orbitals left after deleting core'
  end if
end if
nmax = max(6,ncontrac(0))
do lrun=1,Lhigh
  nmax = max(nmax,ncontrac(lrun))
end do
ncont4 = nmax*nmax*nmax*nmax

call mma_deallocate(nOff)

return

777 format('  Orbital number ',I4,' is the ',I3,'th of L-value ',I3,' it will be removed !!!')

end subroutine ReadBas
