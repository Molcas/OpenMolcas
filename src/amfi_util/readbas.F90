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

subroutine ReadBas(Lhigh,makemean,bonn,breit,symmetry,sameorb,AIMP,oneonly,ncont4,numballcart,LUIN,ifinite)
! Supposed to read the maximum of l-values, the number of primitive and
! contracted functions, the exponents and contraction coefficients

use AMFI_global, only: charge, cntscrtch, Exp_finite, exponents, icore, ikeeplist, ikeeporb, incrLM, ipow2ired, ipowxyz, iredLM, &
                       iredoffunctnew, itotalperIR, Lmax, Loffunction, Moffunction, MxcontL, MxprimL, ncontrac, nprimit, &
                       nrtofiperIR, numbofsym, shiftIRED, shiftIRIR
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: iwp, u6

implicit none
integer(kind=iwp), intent(out) :: Lhigh, ncont4, numballcart, ifinite
logical(kind=iwp), intent(out) :: MakeMean, Bonn, Breit, SameOrb, AIMP, OneOnly
character(len=4), intent(out) :: Symmetry
integer(kind=iwp), intent(in) :: LUIN
integer(kind=iwp) :: I, iBeginIRed(8), icart, iDelperSym(8), ILINE, inired, iorbrun, ired, ired1, ired2, iredrun, Irun, ishifter, &
                     isum, itype, JRUN, lDel, Lrun, Lval, mRun, nfunctperIRED(8), nmax, nsymrun, numbofcart, numbprev, numbr
character(len=54) :: Stars
character(len=4) :: Word
integer(kind=iwp), allocatable :: IREDoffunction(:), nfunctions(:,:), nmbMperIRL(:,:), nOff(:,:)
#ifdef _DEBUGPRINT_
character(len=21) :: chCharge
#define _TEST_ .true.
#else
#define _TEST_ .false.
#endif
logical(kind=iwp), parameter :: IfTest = _TEST_

#ifdef _DEBUGPRINT_
chCharge = '  Charge of nucleus: '
#endif
Stars = '******************************************************'
Bonn = .false.
Breit = .false.
SameOrb = .false.
AIMP = .false.
OneOnly = .false.
MakeMean = .true.
ifinite = 0

if (IfTest) then
  write(u6,'(/,/,/,24X,A)') Stars
  write(u6,'(24X,2A)') '******** Starting Atomic Spin-Orbit MF code ********'
  write(u6,'(24X,A,/,/)') Stars
end if

iCore(:) = 0
call RdNLst(LUIN,'AMFI')
do
  read(LUIN,'(A4)') Word
  if (IfTest) write(u6,'(A4)') Word
  call UpCase(Word)
  select case (WORD)
    case ('BONN')
      Bonn = .true.
    case ('BREI')
      Breit = .true.
    case ('FINI')
      iFinite = 1
      read(LUIN,*) Exp_finite
    case ('SAME')
      SameOrb = .true.
    case ('AIMP')
      AIMP = .true.
      read(LUIN,*) lDel,(iCore(I),I=0,ldel)
      if (IfTest) then
        write(u6,*)
        write(u6,*) 'CORE to be deleted '
        write(u6,*) '   L   #orbs.  '
        write(u6,*)
        do I=0,lDel
          write(u6,'(2I5)') I,iCore(I)
        end do
      end if
      exit
    case ('ONEO')
      OneOnly = .true.
      write(u6,*) ' Only one-electron integrals!!'
      write(u6,*) ' Probably useful for test-purposes only'
    case default
      exit
  end select
end do

if (IfTest) then
  write(u6,*) ' AMFI: '
  if (BONN) then
    write(u6,*) ' Bonn-approach for spin-other-orbit part'
  end if
  if (BREIT) then
    write(u6,*) ' Breit-Pauli type of the SO operator'
  else
    write(u6,*) ' Douglas-Kroll type of the SO operator'
  end if
  if (iFinite == 0) then
    write(u6,*) ' Point nucleus '
  else
    write(u6,*) ' Finite nucleus'
  end if
end if

Symmetry = 'D2H'
NumbofSym = 8
if (IfTest) then
  write(u6,*) ' Symmetry is D2H'
  if (SameOrb) then
    write(u6,*) ' Same-Orbit only'
  else
    write(u6,*) ' Other-Orbit included'
  end if
end if
read(LUIN,*) Charge,Lhigh
if (Lhigh > Lmax) then
  write(u6,*) ' Sorry, so far this code deals only with maximum l-values of ',Lmax
  call Abend()
end if
#ifdef _DEBUGPRINT_
write(u6,'(A21,F5.2)') chCharge,Charge
#endif
call InitiRed(Symmetry)
call mma_allocate(nfunctions,[1,numbofsym],[0,Lhigh],label='nfunctions')
call mma_allocate(nmbMperIRL,[1,numbofsym],[0,Lhigh],label='nmbMperIRL')
nmbMperIRL(:,:) = 0
if (IfTest) write(u6,'(/,A)') '  Used SOC basis set: '
do Lrun=0,Lhigh
  read(LUIN,*) nprimit(Lrun),ncontrac(Lrun)
  if (IfTest) then
    write(u6,'(/,A,I2,A,I2)') '  nExp: ',nprimit(Lrun),' lAng: ',lRun
    write(u6,'(I3,I3)') nprimit(Lrun),ncontrac(Lrun)
  end if
  if (nprimit(Lrun) > MxprimL) then
    write(u6,*) 'Too many primitives for L=',Lrun, &
                ' increase MxprimL in amfi_global or reduce the number of primitives to at least ',MxprimL
    call Abend()
  end if
  if (ncontrac(Lrun) > MxcontL) then
    write(u6,*) ' Too many contracted functions for L=',Lrun, &
                ' increase MxcontL in amfi_global or reduce the number of contracted functions to at most ',MxcontL
    call Abend()
  end if
  if (ncontrac(Lrun) > nprimit(Lrun)) then
    write(u6,*) ' You have more contracted than uncontracted functions, I do not believe that. Sorry! '
    call Abend()
  end if

  ! Read input in MOLCAS-style

  read(LUIN,*) (exponents(ILINE,Lrun),ILINE=1,nprimit(Lrun))
  do ILINE=1,nprimit(Lrun)
    read(LUIN,*) (cntscrtch(ILINE,JRUN,Lrun),Jrun=1,ncontrac(Lrun))
  end do

  ! End of reading for the current L-value

  if (IfTest) then
    write(u6,'(5ES18.8)') (exponents(ILINE,Lrun),ILINE=1,nprimit(Lrun))
    do Irun=1,ncontrac(Lrun)
      write(u6,*) ' orbital : ',irun
      write(u6,'(6(1X,F12.7))') (cntscrtch(I,Irun,Lrun),I=1,nprimit(Lrun))
    end do
  end if

  ! Setting the numbers of cartesians per IR

  do iRedRun=1,NumbofSym
    nFunctions(iRedRun,Lrun) = 0
  end do
  do mRun=-Lrun,Lrun
    nfunctions(ipow2ired(ipowxyz(1,mrun,Lrun),ipowxyz(2,mrun,Lrun),ipowxyz(3,mrun,Lrun)),Lrun) = &
      nfunctions(ipow2ired(ipowxyz(1,mrun,Lrun),ipowxyz(2,mrun,Lrun),ipowxyz(3,mrun,Lrun)),Lrun)+ncontrac(Lrun)
  end do
  do mRun=-Lrun,Lrun
    nmbMperIRL(ipow2ired(ipowxyz(1,mrun,Lrun),ipowxyz(2,mrun,Lrun),ipowxyz(3,mrun,Lrun)),Lrun) = &
      nmbMperIRL(ipOw2ired(ipowxyz(1,mrun,Lrun),ipowxyz(2,mrun,Lrun),ipowxYz(3,mrun,Lrun)),Lrun)+1
  end do
  if (IfTest) write(u6,'(A,8I4)') ' Number of functions per IR: ',(nfunctions(iredrun,Lrun),iredrun=1,numbofsym)
end do   ! End Do for loop over L-values

if (IfTest) then
  write(u6,*) ' Distribution of M-values'
  do Lrun=0,Lhigh
    write(u6,*) nmbMperIRL(:,Lrun)
  end do
end if

numbofcart = 0
do lrun=0,Lhigh
  numbofcart = numbofcart+(Lrun+Lrun+1)*ncontrac(Lrun)
end do

call mma_allocate(nOff,numbofcart,2,Label='nOff')

nfunctperIRED(1:numbofsym) = 0
do Lrun=0,Lhigh
  nfunctperIRED(1:numbofsym) = nfunctperIRED(1:numbofsym)+nfunctions(1:numbofsym,Lrun)
end do
call mma_deallocate(nfunctions)
call mma_deallocate(nmbMperIRL)
if (IfTest) write(u6,'(A,8I3)') ' Total number of atomic functions per IRED ',(nfunctperIRED(iredrun),iredrun=1,numbofsym)
itotalperIR(1:numbofsym) = nfunctperIRED(1:numbofsym)
isum = 0
do iredrun=1,numbofsym
  isum = isum+itotalperIR(iredrun)
end do
numballcart = isum
iorbrun = 0
do iredrun=1,numbofsym
  do inired=1,itotalperIR(iredrun)
    iorbrun = iorbrun+1
    IREDoffunctnew(iorbrun) = iredrun
  end do
end do
if (IfTest) then
  write(u6,'(A,8I3)') 'including additional functions per IRED ',(itotalperIR(iredrun),iredrun=1,numbofsym)
end if
ibeginIRED(1:numbofsym) = 0
do lrun=0,Lhigh
  do mrun=-lrun,lrun
    iredLM(mrun,lrun) = ipow2ired(ipowxyz(1,mrun,Lrun),ipowxyz(2,mrun,Lrun),ipowxyz(3,mrun,Lrun))
    incrLM(mrun,lrun) = ibeginIRED(iredLM(mrun,lrun))
    ibeginIRED(iredLM(mrun,lrUn)) = ibeginIRED(iredLM(mrun,lrun))+ncontrac(lrun)
  end do
end do
if (IfTest) then
  do lrun=0,Lhigh
    write(u6,'(A,I4,A,21I3)') 'L= ',lrun,' shifts inside the IRED',(incrLM(mrun,lrun),mrun=-lrun,lrun)
  end do
end if
shiftIRED(1) = 0
do iredrun=1,numbofsym-1
  shiftIRED(iredrun+1) = shiftIRED(iredrun)+itotalperIR(iredrun)
end do
if (IfTest) then
  write(u6,'(A,8I4)') 'shifts for the IREDs ',(shiftIRED(iredrun),iredrun=1,numbofsym)
  do lrun=0,Lhigh
    do mrun=-Lrun,Lrun
      do irun=1,ncontrac(lrun)
        write(u6,*) 'L,M,contr funct, absolute number ',lrun,mrun,irun,shiftired(iredLM(mrun,lrun))+incrLM(mrun,Lrun)+irun
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
call mma_allocate(IREDoffunction,numbofcart,label='IREDoffunction')
do lrun=0,Lhigh
  do Mrun=-Lrun,Lrun
    ired = iredLM(Mrun,Lrun)
    ishifter = shiftIRED(ired)+incrLM(mrun,lrun)
    do icart=1,ncontrac(Lrun)
      Moffunction(ishifter+icart) = Mrun
      Loffunction(ishifter+icart) = Lrun
      IREDoffunction(ishifter+Icart) = ired
      nOff(ishifter+Icart,2) = icart
    end do
  end do
end do
do irun=1,numbofcart
  nOff(irun,1) = irun
end do
idelpersym(1:numbofsym) = 0
nrtofiperIR(1:numbofsym) = itotalperIR(1:numbofsym)
if (AIMP) then

  ! Generate list of orbitals to be removed

  if (IfTest) write(u6,'(/,A)') '  Core removed for use with AIMP'
  ikeeporb = 0
  numbprev = 0
  do irun=1,numbofcart
    do
      if ((irun == 1) .or. ((irun >= 2) .and. (noff(irun,1) == numbprev+1))) then
        Lval = Loffunction(irun)
        numbr = nOff(irun,1)
        itype = nOff(irun,2)
        if (itype <= iCore(lval)) then
          write(u6,777) numbr,itype,lval
          idelpersym(IREDoffunction(irun)) = idelpersym(IREDoffunction(irun))+1
          numbprev = numbr
        else
          ikeeporb = ikeeporb+1
          ikeeplist(ikeeporb) = numbr
          numbprev = numbr
        end if
        exit
      end if
      ikeeporb = ikeeporb+1
      ikeeplist(ikeeporb) = numbprev+1
      numbprev = numbprev+1
    end do
  end do
  ikeeporb = 0
  nrtofiperIR(1:numbofsym) = itotalperIR(1:numbofsym)-idelpersym(1:numbofsym)
  do nsymrun=1,numbofsym
    ikeeporb = ikeeporb+nrtofiperIR(nsymrun)
  end do
  if (IfTest) then
    write(u6,'(A,8I3)') '  Number of funct. per IRED after removing core: ',(nrtofiperIR(iredrun),iredrun=1,numbofsym)
    write(u6,'(I4,A)') ikeeporb,' orbitals left after deleting core'
  end if
end if
call mma_deallocate(IREDoffunction)
nmax = max(6,ncontrac(0))
do lrun=1,Lhigh
  nmax = max(nmax,ncontrac(lrun))
end do
ncont4 = nmax*nmax*nmax*nmax

call mma_deallocate(nOff)

return

777 format('  Orbital number ',I4,' is the ',I3,'th of L-value ',I3,' it will be removed !!!')

end subroutine ReadBas
