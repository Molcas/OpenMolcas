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
! Copyright (C) 1993, Markus P. Fuelscher                              *
!               1993, Per Ake Malmqvist                                *
!               Pavel Neogrady                                         *
!***********************************************************************

subroutine RdInpPN(run_triples,run_sort)
!***********************************************************************
!                                                                      *
!     purpose:                                                         *
!     - read input                                                     *
!     - set defaults                                                   *
!                                                                      *
!     calling parameters: none                                         *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
!     written by:                                                      *
!     M.P. Fuelscher and P.-AA. Malmqvist                              *
!     University of Lund, Sweden, 1993                                 *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
!     history: none                                                    *
!     reduced by P.N.                                                  *
!                                                                      *
!***********************************************************************

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp, u6

implicit none
logical(kind=iwp) :: run_triples, run_sort
#include "ccsort.fh"
#include "reorg.fh"
#include "motra.fh"
integer(kind=iwp) :: i, IAD15, iCmd, istatus, isym, J, jCmd, LROOTS, LuSpool, M, N, nhelp, ntAsh
character(len=72) :: Line
character(len=4) :: Command
real(kind=wp), allocatable :: Weights(:)
character(len=LenIn8), allocatable :: CName(:)
character(len=4), parameter :: Cmd(20) = ['TITL','END ','CCSD','CCT ','CLOS','OPEN','FROZ','DELE','PRIN','NOOP','IOKE','ZROF', &
                                          'DENO','SHIF','ACCU','ADAP','EXTR','TRIP','NOSO','ITER']

!---  Initialize -------------------------------------------------------*
LROOT = 0
MAXIT = 0
CONV = 1.0e-6_wp
THRSHN = 1.0e-10_wp
THRSHS = 1.0e-8_wp
THRSHF = 0.05_wp
ORBIN = 'DEFAULT '
THRENE = 1.5_wp
ORBIT = 'DEFAULT '
THROCC = Zero
FOCKTYPE = 'STANDARD'
HZERO = 'STANDARD'
METHOD = 'CONJ'
IFJAC = 0
RFpert = .false.
NTIT = 0

lunsta = 21
luna1 = 22
luna2 = 23
luna3 = 24
luna4 = 25
lunab = 50
lunt3 = 26
lunda1 = 9
lunda2 = 10
lunpublic = 29

!---  Open JOBIPH and LUONEM files ------------------------------------*

! Job interface
JOBIPH = 15
! Job interface
call DANAME(JOBIPH,'JOBIPH')

!---  Read input from JOBIPH file -------------------------------------*
IAD15 = 0
call iDAFILE(JOBIPH,2,IADR15,15,IAD15)
!DIVNUO
IAD15 = IADR15(1)
call mma_allocate(CName,mxOrb,label='CName')
call mma_allocate(Weights,mxRoot,label='Weights')
call WR_RASSCF_Info(JOBIPH,2,iAd15,NACTEL,ISPIN,NSYM,LSYM,NFRO,NISH,NASH,NDEL,NBAS,8,CName,LenIn8*MxOrb,NCONF,HEADER,2*72,TITLE, &
                    4*18*MXTIT,POTNUC,LROOTS,NROOTS,IROOT,MXROOT,NRAS1,NRAS2,NRAS3,NHOLE1,NELE3,IPT2,Weights)
call mma_deallocate(CName)
call mma_deallocate(Weights)

! define defaults for REORG

ntAsh = 0
do isym=1,nSym
  ntAsh = ntAsh+nAsh(iSym)
end do
cckey = 1
t3key = 1
clopkey = 1
if (ntAsh == 0) clopkey = 2
ndelr(1:nsym) = NDEL(1:nsym)
nfror(1:nsym) = NFRO(1:nsym)
!GG fullprint = 0
noop = 0
iokey = 1
zrkey = 1
run_triples = .true.
run_sort = .true.

!---  Read input from TRAONE file -------------------------------------*
call RdTraOne()
nFror(1:nSym) = nFroX(1:nSym)
nDelr(1:nSym) = nDelX(1:nSym)

!---  Read input from LuSpool -----------------------------------------*
LuSpool = 17
call SpoolInp(LuSpool)
rewind(LuSpool)
Command = '&REO'
call RdNlst(LuSpool,'CCSDT')
do
  read(LuSpool,'(A)',iostat=istatus) Line
  if (istatus < 0) call Error(1)
  if ((Line(1:1) /= '*') .and. (Line /= '')) exit
end do
Command = Line(1:4)
call UpCase(Command)
jCmd = 0
do iCmd=1,size(Cmd)
  if (Command == Cmd(iCmd)) jCmd = iCmd
end do
if (jCmd == 0) call Error(2)
do
  select case (jCmd)
    case default !(1) !TITL
      !---  process TITLE command -------------------------------------*
      do
        read(LuSpool,'(A)',iostat=istatus) Line
        if (istatus < 0) call Error(1)
        Command = Line(1:4)
        call UpCase(Command)
        if (Command(1:1) == '*') cycle
        jCmd = 0
        do iCmd=1,size(Cmd)
          if (Command == Cmd(iCmd)) jCmd = iCmd
        end do
        if (jCmd /= 0) exit
        if (nTit >= mxTit) cycle
        nTit = nTit+1
        if (nTit <= mxTit) read(Line,'(18A4)') (Title(nTit,i),i=1,18)
      end do
    case (2) !END
      exit
    case (3) !CCSD
      !---  CCSD  command ---------------------------------------------*
      cckey = 1
      t3key = 0
      run_triples = .false.
      jCmd = 1
    case (4) !CCT
      !---  CCT  command ----------------------------------------------*
      cckey = 1
      t3key = 1
      run_triples = .true.
      jCmd = 1
    case (5) !CLOS
      !---  CLOSed  command -------------------------------------------*
      clopkey = 2
      jCmd = 1
    case (6) !OPEN
      !---  OPEN  command ---------------------------------------------*
      clopkey = 1
      jCmd = 1
    case (7) !FROZ
      !---  FROZen  command -------------------------------------------*
      read(LuSpool,*) (nfror(nhelp),nhelp=1,nsym)
      jCmd = 1
    case (8) !DELE
      !---  DELEte  command -------------------------------------------*
      read(LuSpool,*) (ndelr(nhelp),nhelp=1,nsym)
      jCmd = 1
    case (9) !PRIN
      !---  PRINt  command --------------------------------------------*
      read(LuSpool,*) fullprint
      jCmd = 1
    case (10) !NOOP
      !---  NOOPeration  command --------------------------------------*
      noop = 1
      jCmd = 1
    case (11) !IOKE
      !---  IOKEy  command --------------------------------------------*
      read(LuSpool,*) iokey
      if ((iokey < 1) .or. (iokey > 2)) iokey = 2
      jCmd = 1
    case (12) !ZROF
      !---  ZROFf  command --------------------------------------------*
      zrkey = 0
      jCmd = 1
    case (19) !NOSO
      !---  NOSOrt  command -------------------------------------------*
      run_sort = .false.
      jCmd = 1
    case (13:18,20) !DENO, SHIF, ACCU, ADAP, EXTR, TRIP, ITER
      !---  DENO  command ---------------------------------------------*
      !---  SHIF  command ---------------------------------------------*
      !---  ACCU  command ---------------------------------------------*
      !---  ADAP  command ---------------------------------------------*
      !---  EXTR  command ---------------------------------------------*
      !---  TRIP  command ---------------------------------------------*
      !---  ITER  command ---------------------------------------------*
      jCmd = 1
  end select
end do

!---  The end of the input section, complete input processing ---------*
NFROT = 0
NISHT = 0
NASHT = 0
NRAS1T = 0
NRAS2T = 0
NRAS3T = 0
NOSHT = 0
NSSHT = 0
NDELT = 0
NORBT = 0
NBAST = 0
NBAS2 = 0
NORB1 = 0
do ISYM=1,NSYM
  NIES(ISYM) = NISHT
  NAES(ISYM) = NASHT
  NSES(ISYM) = NSSHT
  NOSH(ISYM) = NISH(ISYM)+NASH(ISYM)
  NSSH(ISYM) = NBAS(ISYM)-NFRO(ISYM)-NOSH(ISYM)-NDEL(ISYM)
  NORB(ISYM) = NOSH(ISYM)+NSSH(ISYM)
  NORBT = NORBT+NORB(ISYM)
  NBAS2 = NBAS2+NBAS(ISYM)**2
  NORB1 = NORB1+(NORB(ISYM)**2+NORB(ISYM))/2
  NFROT = NFROT+NFRO(ISYM)
  NISHT = NISHT+NISH(ISYM)
  NASHT = NASHT+NASH(ISYM)
  NOSHT = NOSHT+NOSH(ISYM)
  NRAS1T = NRAS1T+NRAS1(ISYM)
  NRAS2T = NRAS2T+NRAS2(ISYM)
  NRAS3T = NRAS3T+NRAS3(ISYM)
  NSSHT = NSSHT+NSSH(ISYM)
  NDELT = NDELT+NDEL(ISYM)
  NBAST = NBAST+NBAS(ISYM)
end do

!---  Identify the wave function type ---------------------------------*
ISCF = 0
if (NASHT == 0) ISCF = 1
if (NACTEL == 2*NASHT) ISCF = 1
if ((iSpin > 1) .and. (ISPIN == NACTEL+1) .and. (NACTEL == NASHT)) ISCF = 2

! test agreement between REORG input and JOBIPH
!if (clopkey /= (3-ISCF)) then
!  write(u6,*) ' Diference in closed/open specification'
!  write(u6,*) ' Plaese, correct REORG input file'
!  write(u6,*) clopkey,ISCF
!  call Quit(16)
!end if
!
!  Should not be necessary anymore, let's just hardwire it in.
!  This will save a keyword
clopkey = 3-ISCF
!---  Identify the reference function ---------------------------------*
if (ISCF > 0) then
  LROOT = 1
  NROOTS = 1
  IROOT(1) = 1
end if
if ((LROOT == 0) .and. (NROOTS == 1)) LROOT = IROOT(1)

!---  Create the symmetry multiplication table ------------------------*
MUL(1,1) = 1
M = 1
do N=1,3
  do I=1,M
    do J=1,M
      MUL(I+M,J) = M+MUL(I,J)
      MUL(I,J+M) = MUL(I+M,J)
      MUL(I+M,J+M) = MUL(I,J)
    end do
  end do
  M = 2*M
end do
call Close_LuSpool(LuSpool)

!---  Exit ------------------------------------------------------------*
return

contains

!---  Error exits -----------------------------------------------------*
subroutine Error(code)

  integer :: code

  write(u6,*)
  select case (code)
    case (1)
      write(u6,*) ' *** input error ***'
      write(u6,*) ' hitting end of file mark'
    case (2)
      write(u6,*) ' *** input error ***'
      write(u6,*) ' unknown input'
      write(u6,*) ' line: ',Line
  end select
  write(u6,*)
  call Quit_OnUserError()

end subroutine Error

end subroutine RdInpPN
