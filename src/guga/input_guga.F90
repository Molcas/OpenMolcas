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
!***********************************************************************

subroutine INPUT_GUGA(L0,L1,L2,L3,ISPAC)

use guga_global, only: IADD10, ICASE, ICH, IFIRST, ILIM, IPRINT, ISPIN, JRC, LN, LNP, Lu_10, MXVERT, N, NIORB, NSM, NSYM, S
use guga_util_global, only: IAD10, nIOCR
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Half
use Definitions, only: iwp, u5, u6

implicit none
integer(kind=iwp), intent(out) :: L0(4*MXVERT), L1(4*MXVERT), L2(4*MXVERT), L3(4*MXVERT), ISPAC
#include "warnings.h"
integer(kind=iwp), parameter :: mxTit = 10, nCmd = 18
integer(kind=iwp) :: I, ICIALL, iCmd, ICOR(55), IFCORE, IN_, IN1, IN2, IN3, INTNUM, IOM, IONE(8), iOpt, IR, IR1, IR2, iRef, &
                     istatus, ISUM, IVER, J, jCmd, jEnd, JJS(18), JONE(8), jStart, LN1, LN2, LSYM, LV, MN, MX, NACTEL, NCOR(8), &
                     NCORI, NFREF, NISH(8), NISHI, NISHT, nJJS, nJRC, NN, NO, NONE_, NREF, NRLN1, NSH(8), NSHI, nTit, NVAL(8), NVALI
logical(kind=iwp) :: skip
character(len=132) :: ModLine
character(len=72) :: Line, Title(mxTit)
character(len=4) :: Command
integer(kind=iwp), allocatable :: IOCR(:), JREFX(:), JSY(:)
integer(kind=iwp), parameter :: MLL(64) = [1,2,3,4,5,6,7,8,2,1,4,3,6,5,8,7,3,4,1,2,7,8,5,6,4,3,2,1,8,7,6,5,5,6,7,8,1,2,3,4,6,5,8, &
                                           7,2,1,4,3,7,8,5,6,3,4,1,2,8,7,6,5,4,3,2,1]
character(len=4), parameter :: Cmd(nCmd) = ['TITL','ELEC','SPIN','SYMM','ACTI','PRIN','REFE','FIRS','INAC','CIAL','VALE','INTE', &
                                            'NOCO','ONEO','EXTR','NONI','NACT','END ']

!---  Initialize data and variables -----------------------------------*
IOM = 55
IVER = MXVERT
IFIRST = 0
IPRINT = 0
ICIALL = 0
ILIM = 4
N = -1
NACTEL = -1
NISHT = 0
!NSYM = 1
call Get_iScalar('nSym',NSYM)
NREF = 0
!PAM97 New default: Interacting space.
INTNUM = 1
!PAM97 IFCORE /= 0 means core-polarization orbitals (NOCO keyword).
IFCORE = 0
LSYM = 1
do I=1,8
  NISH(I) = 0
  NVAL(I) = 0
  NCOR(I) = 0
  NSH(I) = 0
  IONE(I) = 0
end do
do I=1,55
  ICOR(I) = 0
end do
nTit = 0
jCmd = 0
call mma_allocate(IOCR,nIOCR,label='IOCR')

!---  Read input from standard input ----------------------------------*
call RdNLst(u5,'GUGA')
skip = .false.
do
  if (skip) then
    skip = .false.
  else
    read(u5,'(A)',iostat=istatus) Line
    if (istatus < 0) call Error(1)
    Command = Line(1:8)
    call UpCase(Command)
    if (Command(1:1) == '*') cycle
    if (Command == ' ') cycle
    jCmd = 0
    do iCmd=1,nCmd
      if (Command == Cmd(iCmd)) jCmd = iCmd
    end do
  end if
  select case (jCmd)
    case default
      write(u6,*) 'Input: Illegal Keyword'
      write(u6,'(A,A)') 'Command=',Command
      call Quit(_RC_INPUT_ERROR_)
    case (1)
      !---  process TITLE    command ----------------------------------*
      do
        do
          read(u5,'(A)',iostat=istatus) Line
          if (istatus < 0) call Error(1)
          Command = Line(1:8)
          call UpCase(Command)
          if (Command(1:1) /= '*') exit
        end do
        jCmd = 0
        do iCmd=1,nCmd
          if (Command == Cmd(iCmd)) jCmd = iCmd
        end do
        if (jCmd /= 0) then
          skip = .true.
          exit
        end if
        nTit = nTit+1
        if (nTit <= mxTit) Title(nTit) = Line
      end do
    case (2)
      !---  process ELECTRON command ----------------------------------*
      do
        read(u5,'(A)',iostat=istatus) Line
        if (istatus < 0) call Error(1)
        if (Line(1:1) /= '*') exit
      end do
      read(Line,*,iostat=istatus) N
      if (istatus > 0) call Error(2)
    case (3)
      !---  process SPIN     command ----------------------------------*
      do
        read(u5,'(A)',iostat=istatus) Line
        if (istatus < 0) call Error(1)
        if (Line(1:1) /= '*') exit
      end do
      read(Line,*,iostat=istatus) ISPIN
      if (istatus > 0) call Error(2)
    case (4)
      !---  process SYMMETRY command ----------------------------------*
      write(u6,*) 'Input_GUGA: keyword SYMMETRY is obsolete and ignored!'
      do
        read(u5,'(A)',iostat=istatus) Line
        if (istatus < 0) call Error(1)
        if (Line(1:1) /= '*') exit
      end do
      read(Line,*,iostat=istatus) I
      if (istatus > 0) call Error(2)
    case (5)
      !---  process ACTIVE   command ----------------------------------*
      do
        read(u5,'(A)',iostat=istatus) Line
        if (istatus < 0) call Error(1)
        if (Line(1:1) /= '*') exit
      end do
      ModLine = Line//' 0 0 0 0 0 0 0 0'
      read(ModLine,*,iostat=istatus) (NSH(I),I=1,8)
      if (istatus > 0) call Error(2)
    case (6)
      !---  process PRINT    command ----------------------------------*
      do
        read(u5,'(A)',iostat=istatus) Line
        if (istatus < 0) call Error(1)
        if (Line(1:1) /= '*') exit
      end do
      read(Line,*,iostat=istatus) IPRINT
      if (istatus > 0) call Error(2)
    case (7)
      !---  process REFERENC command ----------------------------------*
      do
        read(u5,'(A)',iostat=istatus) Line
        if (istatus < 0) call Error(1)
        if (Line(1:1) /= '*') exit
      end do
      read(Line,*,iostat=istatus) nRef,LN1
      if (istatus > 0) call Error(2)
      if (LN1 /= 0) then
        jEnd = 0
        do iRef=1,nRef
          jStart = jEnd+1
          jEnd = jEnd+LN1
          read(u5,'(80I1)',iostat=istatus) (IOCR(j),j=jStart,jEnd)
          if (istatus < 0) call Error(1)
          if (istatus > 0) call Error(2)
        end do
      end if
    case (8)
      !---  process FIRST    command ----------------------------------*
      IFIRST = 1
      ILIM = 2
    case (9)
      !---  process INACTIVE command ----------------------------------*
      do
        read(u5,'(A)',iostat=istatus) Line
        if (istatus < 0) call Error(1)
        if (Line(1:1) /= '*') exit
      end do
      ModLine = Line//' 0 0 0 0 0 0 0 0'
      read(ModLine,*,iostat=istatus) (NISH(I),I=1,8)
      if (istatus > 0) call Error(2)
      NISHT = 0
      do I=1,8
        NISHT = NISHT+NISH(I)
      end do
    case (10)
      !---  process CIALL    command ----------------------------------*
      do
        read(u5,'(A)',iostat=istatus) Line
        if (istatus < 0) call Error(1)
        if (Line(1:1) /= '*') exit
      end do
      read(Line,*,iostat=istatus) LSYM
      if (istatus > 0) call Error(2)
      ICIALL = 1
    case (11)
      !---  process VALENCE  command ----------------------------------*
      do
        read(u5,'(A)',iostat=istatus) Line
        if (istatus < 0) call Error(1)
        if (Line(1:1) /= '*') exit
      end do
      ModLine = Line//' 0 0 0 0 0 0 0 0'
      read(ModLine,*,iostat=istatus) (NVAL(I),I=1,8)
      if (istatus > 0) call Error(2)
    case (12)
      !---  process INTERACT command ----------------------------------*
      INTNUM = 1
    case (13)
      !---  process NOCORR   command ----------------------------------*
      do
        read(u5,'(A)',iostat=istatus) Line
        if (istatus < 0) call Error(1)
        if (Line(1:1) /= '*') exit
      end do
      ModLine = Line//' 0 0 0 0 0 0 0 0'
      read(ModLine,*,iostat=istatus) (NCOR(I),I=1,8)
      if (istatus > 0) call Error(2)
      !PAM97 IFCORE was not set -- assume bug. Following line inserted:
      IFCORE = 1
    case (14)
      !---  process ONEOCC   command ----------------------------------*
      do
        read(u5,'(A)',iostat=istatus) Line
        if (istatus < 0) call Error(1)
        if (Line(1:1) /= '*') exit
      end do
      ModLine = Line//' 0 0 0 0 0 0 0 0'
      read(ModLine,*,iostat=istatus) (IONE(I),I=1,8)
      if (istatus > 0) call Error(2)
    case (15)
      !---  process EXTRACT  command ----------------------------------*
      write(u6,*) 'Input: EXTRACT option is redundant and is ignored!'
    case (16)
      !---  process NON-INTERACT command ------------------------------*
      INTNUM = 0
    case (17)
      !---  process NACTEL       command ------------------------------*
      do
        read(u5,'(A)',iostat=istatus) Line
        if (istatus < 0) call Error(1)
        if (Line(1:1) /= '*') exit
      end do
      read(Line,*,iostat=istatus) NACTEL
      if (istatus > 0) call Error(2)
    case (18)
      exit
  end select
end do
!---  The end of the input is reached, print the title ----------------*
if (ntit == 0) then
  ntit = 1
  title(1) = ' (No title was given)'
end if

! Nr of correlated electrons:
if ((N == -1) .and. (NACTEL == -1)) then
  write(u6,*) ' Neither the number of correlated electrons (Keyword ELECTRONS)'
  write(u6,*) ' nor the nr of active electrons in the reference space (NACTEL) has been specified.'
  write(u6,*) ' The number of active electrons are set to zero.'
  NACTEL = 0
else if ((N > -1) .and. (NACTEL > -1)) then
  write(u6,*) ' Both the number of correlated electrons (Keyword ELECTRONS)'
  write(u6,*) ' and the nr of active electrons in the reference space (NACTEL) have been specified.'
  if (N /= 2*NISHT+NACTEL) then
    N = 2*NISHT+NACTEL
    write(u6,*) ' Number of correlated electrons is recomputed,=',N
  end if
end if
if (N == -1) N = 2*NISHT+NACTEL
if (NACTEL == -1) NACTEL = N-2*NISHT

write(u6,*)
write(u6,'(6X,A)') repeat('*',120)
write(u6,'(6X,A,A,A)') '*',repeat(' ',118),'*'
write(u6,'(6X,A,A,A6,A,A)') '*',repeat(' ',56),'Title:',repeat(' ',56),'*'
do i=1,nTit
  call Center_Text(Title(i))
  write(u6,'(6X,A,A,A72,A,A)') '*',repeat(' ',23),Title(i),repeat(' ',23),'*'
end do
write(u6,'(6X,A,A,A)') '*',repeat(' ',118),'*'
write(u6,'(6X,A)') repeat('*',120)
write(u6,*)
S = (ISPIN-1)*Half
if (IFIRST == 0) write(u6,2)
if (IFIRST /= 0) write(u6,1)
write(u6,110) N,S
write(u6,109) (I,I=1,NSYM)
write(u6,106) (NISH(I),I=1,NSYM)
write(u6,108) (NSH(I),I=1,NSYM)
write(u6,208) (NVAL(I),I=1,NSYM)
write(u6,206) (NCOR(I),I=1,NSYM)
write(u6,209) (IONE(I),I=1,NSYM)
LN = 0
LV = 0
NIORB = 0
do I=1,NSYM
  NIORB = NIORB+NISH(I)
  LN = LN+NSH(I)+NISH(I)
  LV = LV+NVAL(I)
end do
NONE_ = 0
IN_ = LV+NIORB
do I=1,NSYM
  NN = IONE(I)
  do NO=1,NN
    NONE_ = NONE_+1
    IN_ = IN_+1
    JONE(NONE_) = IN_
  end do
  IN_ = IN_+NSH(I)-NN
end do
LN = LN+LV
if (ICIALL == 1) LN1 = LN-LV-NIORB
if (LN /= LN1+LV+NIORB) then
  write(u6,*) 'Input: LN /= LN1+LV+NIORB'
  write(u6,*) 'LN,LN1,LV,NIORB=',LN,LN1,LV,NIORB
  call Quit(_RC_INPUT_ERROR_)
end if
LNP = LN*(LN+1)/2
IN_ = 0
IN3 = 0
IN1 = LV
IN2 = NIORB+LV
do I=1,NSYM
  NISHI = NISH(I)
  NSHI = NSH(I)
  NVALI = NVAL(I)
  NCORI = NCOR(I)
  do J=1,NVALI
    IN3 = IN3+1
    IN_ = IN_+1
    NSM(IN3) = I
    ICH(IN_) = IN3
  end do
  do J=1,NISHI
    IN1 = IN1+1
    IN_ = IN_+1
    NSM(IN1) = I
    ICH(IN_) = IN1
    if (J <= NCORI) ICOR(IN1) = 1
  end do
  do J=1,NSHI
    IN2 = IN2+1
    IN_ = IN_+1
    NSM(IN2) = I
    ICH(IN_) = IN2
  end do
end do
if (LN > IOM) then
  write(u6,*) 'Input: LN > IOM'
  write(u6,*) 'LN,IOM=',LN,IOM
  call Quit(_RC_INPUT_ERROR_)
end if
call TAB2F(IVER-1,LV)
call TAB2(NREF,IOCR,nIOCR,L0,L1,L2,L3,INTNUM,LV,LSYM,ICIALL,IFCORE,ICOR,NONE_,JONE)
LN2 = LN1
if (LN1 > 8) LN2 = 16
if (LN1 == 0) then
  write(u6,55)
else
  write(u6,107) (I,I=1,LN2)
  NO = N-2*NIORB
  MX = 0
  do IREF=1,NREF
    MN = MX+1
    MX = MX+LN1
    write(u6,112) IREF,(IOCR(J),J=MN,MX)

    ! Sum up the occupation numbers of the first reference:
    ISUM = 0
    do I=1,LN1
      ISUM = ISUM+IOCR(MN+I-1)
    end do
    if (ISUM /= NO) then
      write(u6,*) ' Summed occupation nums of this reference does'
      write(u6,*) ' not match nr of electrons.'
      write(u6,*) ' In closed shells: 2*NIORB=',2*NIORB
      write(u6,*) ' Summed occupation nums   =',NO
      write(u6,*) ' Sum total is             =',2*NIORB+NO
      write(u6,*) ' But input says nr of elec=',N
      call Quit(_RC_INPUT_ERROR_)
    end if
  end do
end if

! Here with ILIM=2 (FIRST command) or 4 (normal, default).
call mma_allocate(JSY,3000,label='JSY')
call mma_allocate(JREFX,9000,label='JREFX')
call CONFIG(NREF,IOCR,nIOCR,L0,L1,L2,L3,JSY,INTNUM,LSYM,JJS,LV,IFCORE,ICOR,NONE_,JONE,JREFX,NFREF)
IR = JRC(ILIM)
ISPAC = IR*LNP
if (IPRINT >= 2) write(u6,9) ISPAC
NRLN1 = NREF*LN1
if (LN1 == 0) NRLN1 = 1
!PAM97 IR1 = (LN*IR+29)/30
IR1 = (LN*IR+14)/15
IR2 = (IR+9)/10

iOpt = 1
nJJS = size(JJS)
nJRC = size(JRC)
call WR_GUGA(Lu_10,iOpt,IADD10,NFREF,S,N,LN,NSYM,IR1,IR2,IFIRST,INTNUM,LSYM,NREF,LN1,NRLN1,NSH,NISH,8,JRC,nJRC,JJS,nJJS,NVAL,IOCR, &
             nIOCR)
call iDAFILE(Lu_10,1,ICASE,IR1,IADD10)
call iDAFILE(Lu_10,1,JSY,IR2,IADD10)
IAD10(2) = IADD10
call iDAFILE(Lu_10,1,JREFX,JRC(1),IADD10)
call mma_deallocate(IOCR)
call mma_deallocate(JREFX)
call mma_deallocate(JSY)

return

1 format(//,6X,'ONLY SINGLE REPLACEMENTS INCLUDED')
2 format(//,6X,'ALL SINGLE AND DOUBLE REPLACEMENTS')
9 format(//,6X,'ELEMENTS TO BE SORTED',I7)
55 format(//,6X,'ONE CLOSED SHELL REFERENCE STATE')
106 format(6X,'INACTIVE',8I5)
107 format(//,6X,'OCCUPATION OF REFERENCE STATES',//,6X,'REF.STATE',2X,'ORB:',I2,15I4)
108 format(6X,'ACTIVE  ',8I5)
109 format(//,14X,'ORBITALS PER SYMMETRY',/,14X,8I5)
110 format(/,6X,'NUMBER OF ELECTRONS IN CI',I10,/,6X,'TOTAL SPIN QUANTUM NUMBER',F10.2)
112 format(6X,I5,8X,16I4)
206 format(6X,'CORE    ',8I5)
208 format(6X,'VALENCE ',8I5)
209 format(6X,'ONEOCC  ',8I5)

contains

subroutine Error(code)

  integer(kind=iwp), intent(in) :: code

  if (code == 1) then
    write(u6,*) 'Input: End of input file encountered'
  else if (code == 2) then
    write(u6,*) 'Input: Error while reading input!'
  end if
  write(u6,'(A,A)') 'Last Command: ',Command
  call Quit(_RC_INPUT_ERROR_)

end subroutine Error

end subroutine INPUT_GUGA
