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

subroutine INPUT_GUGA(ISO,JSYM,JSY,L0,L1,L2,L3,ISPAC)

implicit real*8(A-H,O-Z)
#include "niocr.fh"
#include "SysDef.fh"
#include "files_addr.fh"
dimension ISO(*), JSYM(*), JSY(*), L0(*), L1(*), L2(*), L3(*)
#include "real_guga.fh"
#include "integ.fh"
#include "warnings.h"
dimension MLL(64), IOCR(nIOCR), JREFX(9000), NISH(8), JJS(18), NVAL(8), NCOR(8), ICOR(55), IONE(8), JONE(8)
parameter(nCmd=18)
parameter(mxTit=10)
character*4 Command, Cmd(nCmd)
character*72 Line, Title(mxTit)
character*132 ModLine
data MLL/1,2,3,4,5,6,7,8,2,1,4,3,6,5,8,7,3,4,1,2,7,8,5,6,4,3,2,1,8,7,6,5,5,6,7,8,1,2,3,4,6,5,8,7,2,1,4,3,7,8,5,6,3,4,1,2,8,7,6,5, &
         4,3,2,1/
data Cmd/'TITL','ELEC','SPIN','SYMM','ACTI','PRIN','REFE','FIRS','INAC','CIAL','VALE','INTE','NOCO','ONEO','EXTR','NONI','NACT', &
         'END '/

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
IN = 0
do I=1,8
  NISH(I) = 0
  NVAL(I) = 0
  NCOR(I) = 0
  NSH(I) = 0
  IONE(I) = 0
  do J=1,8
    IN = IN+1
    MUL(I,J) = MLL(IN)
  end do
end do
do I=1,55
  ICOR(I) = 0
end do
nTit = 0

!---  Read input from standard input ----------------------------------*
call RdNLst(5,'GUGA')
10 read(5,'(A)',end=991) Line
Command = Line(1:8)
call UpCase(Command)
if (Command(1:1) == '*') goto 10
if (Command == ' ') goto 10
jCmd = 0
do iCmd=1,nCmd
  if (Command == Cmd(iCmd)) jCmd = iCmd
end do
20 goto(100,200,300,400,500,600,700,800,900,1000,1100,1200,1300,1400,1500,1600,1700,1800) jCmd
write(6,*) 'Input: Illegal Keyword'
write(6,'(A,A)') 'Command=',Command
call Quit(_RC_INPUT_ERROR_)

!---  process TITLE    command ----------------------------------------*
100 continue
read(5,'(A)',end=991) Line
Command = Line(1:8)
call UpCase(Command)
if (Command(1:1) == '*') goto 100
jCmd = 0
do iCmd=1,nCmd
  if (Command == Cmd(iCmd)) jCmd = iCmd
end do
if (jCmd /= 0) goto 20
nTit = nTit+1
if (nTit <= mxTit) Title(nTit) = Line
goto 100

!---  process ELECTRON command ----------------------------------------*
200 continue
read(5,'(A)',end=991) Line
if (Line(1:1) == '*') goto 200
read(Line,*,Err=992) N
goto 10

!---  process SPIN     command ----------------------------------------*
300 continue
read(5,'(A)',end=991) Line
if (Line(1:1) == '*') goto 300
read(Line,*,Err=992) ISPIN
goto 10

!---  process SYMMETRY command ----------------------------------------*
400 continue
write(6,*) 'Input_GUGA: keyword SYMMETRY is obsolete and ignored!'
read(5,'(A)',end=991) Line
if (Line(1:1) == '*') goto 400
read(Line,*,Err=992) I
goto 10

!---  process ACTIVE   command ----------------------------------------*
500 continue
read(5,'(A)',end=991) Line
if (Line(1:1) == '*') goto 500
ModLine = Line//' 0 0 0 0 0 0 0 0'
read(ModLine,*,Err=992) (NSH(I),I=1,8)
goto 10

!---  process PRINT    command ----------------------------------------*
600 continue
read(5,'(A)',end=991) Line
if (Line(1:1) == '*') goto 600
read(Line,*,Err=992) IPRINT
goto 10

!---  process REFERENC command ----------------------------------------*
700 continue
read(5,'(A)',end=991) Line
if (Line(1:1) == '*') goto 700
read(Line,*,Err=992) nRef,LN1
if (LN1 == 0) goto 10
jEnd = 0
do iRef=1,nRef
  jStart = jEnd+1
  jEnd = jEnd+LN1
  read(5,'(80I1)',end=991,Err=992) (IOCR(j),j=jStart,jEnd)
end do
goto 10

!---  process FIRST    command ----------------------------------------*
800 continue
IFIRST = 1
ILIM = 2
goto 10

!---  process INACTIVE command ----------------------------------------*
900 continue
read(5,'(A)',end=991) Line
if (Line(1:1) == '*') goto 900
ModLine = Line//' 0 0 0 0 0 0 0 0'
read(ModLine,*,Err=992) (NISH(I),I=1,8)
NISHT = 0
do I=1,8
  NISHT = NISHT+NISH(I)
end do
goto 10

!---  process CIALL    command ----------------------------------------*
1000 continue
read(5,'(A)',end=991) Line
if (Line(1:1) == '*') goto 1000
read(Line,*,Err=992) LSYM
ICIALL = 1
goto 10

!---  process VALENCE  command ----------------------------------------*
1100 continue
read(5,'(A)',end=991) Line
if (Line(1:1) == '*') goto 1100
ModLine = Line//' 0 0 0 0 0 0 0 0'
read(ModLine,*,Err=992) (NVAL(I),I=1,8)
goto 10

!---  process INTERACT command ----------------------------------------*
1200 continue
INTNUM = 1
goto 10

!---  process NOCORR   command ----------------------------------------*
1300 continue
read(5,'(A)',end=991) Line
if (Line(1:1) == '*') goto 1300
ModLine = Line//' 0 0 0 0 0 0 0 0'
read(ModLine,*,Err=992) (NCOR(I),I=1,8)
!PAM97 IFCORE was not set -- assume bug. Following line inserted:
IFCORE = 1
goto 10

!---  process ONEOCC   command ----------------------------------------*
1400 continue
read(5,'(A)',end=991) Line
if (Line(1:1) == '*') goto 1400
ModLine = Line//' 0 0 0 0 0 0 0 0'
read(ModLine,*,Err=992) (IONE(I),I=1,8)
goto 10

!---  process EXTRACT  command ----------------------------------------*
1500 write(6,*) 'Input: EXTRACT option is redundant and is ignored!'
goto 10

!---  process NON-INTERACT command ------------------------------------*
1600 continue
INTNUM = 0
goto 10

!---  process NACTEL       command ------------------------------------*
1700 continue
read(5,'(A)',end=991) Line
if (Line(1:1) == '*') goto 200
read(Line,*,Err=992) NACTEL
goto 10

!---  The end of the input is reached, print the title ----------------*
1800 continue
if (ntit == 0) then
  ntit = 1
  title(1) = ' (No title was given)'
end if

! Nr of correlated electrons:
if ((N == -1) .and. (NACTEL == -1)) then
  write(6,*) ' Neither the number of correlated electrons (Keyword ELECTRONS)'
  write(6,*) ' nor the nr of active electrons in the reference space (NACTEL) has been specified.'
  write(6,*) ' The number of active electrons are set to zero.'
  NACTEL = 0
else if ((N > -1) .and. (NACTEL > -1)) then
  write(6,*) ' Both the number of correlated electrons (Keyword ELECTRONS)'
  write(6,*) ' and the nr of active electrons in the reference space (NACTEL) have been specified.'
  if (N /= 2*NISHT+NACTEL) then
    N = 2*NISHT+NACTEL
    write(6,*) ' Number of correlated electrons is recomputed,=',N
  end if
end if
if (N == -1) N = 2*NISHT+NACTEL
if (NACTEL == -1) NACTEL = N-2*NISHT

write(6,*)
write(6,'(6X,120A1)') ('*',i=1,120)
write(6,'(6X,120A1)') '*',(' ',i=1,118),'*'
write(6,'(6X,57A1,A6,57A1)') '*',(' ',i=1,56),'Title:',(' ',i=1,56),'*'
do i=1,nTit
  call Center_Text(Title(i))
  write(6,'(6X,24A1,A72,24A1)') '*',(' ',j=1,23),Title(i),(' ',j=1,23),'*'
end do
write(6,'(6X,120A1)') '*',(' ',i=1,118),'*'
write(6,'(6X,120A1)') ('*',i=1,120)
write(6,*)
S = (ISPIN-1)*0.5d0
if (IFIRST == 0) write(IW,2)
2 format(//,6X,'ALL SINGLE AND DOUBLE REPLACEMENTS')
if (IFIRST /= 0) write(IW,1)
1 format(//,6X,'ONLY SINGLE REPLACEMENTS INCLUDED')
write(IW,110) N,S
110 format(/,6X,'NUMBER OF ELECTRONS IN CI',I10,/,6X,'TOTAL SPIN QUANTUM NUMBER',F10.2)
write(IW,109) (I,I=1,NSYM)
109 format(//,14X,'ORBITALS PER SYMMETRY',/,14X,8I5)
write(IW,106) (NISH(I),I=1,NSYM)
106 format(6X,'INACTIVE',8I5)
write(IW,108) (NSH(I),I=1,NSYM)
108 format(6X,'ACTIVE  ',8I5)
write(IW,208) (NVAL(I),I=1,NSYM)
208 format(6X,'VALENCE ',8I5)
write(IW,206) (NCOR(I),I=1,NSYM)
206 format(6X,'CORE    ',8I5)
write(IW,209) (IONE(I),I=1,NSYM)
209 format(6X,'ONEOCC  ',8I5)
LN = 0
LV = 0
NIORB = 0
do I=1,NSYM
  NIORB = NIORB+NISH(I)
  LN = LN+NSH(I)+NISH(I)
  LV = LV+NVAL(I)
end do
NONE = 0
IN = LV+NIORB
do I=1,NSYM
  NN = IONE(I)
  if (NN == 0) GO TO 905
  do NO=1,NN
    NONE = NONE+1
    IN = IN+1
    JONE(NONE) = IN
  end do
905 IN = IN+NSH(I)-NN
end do
LN = LN+LV
if (ICIALL == 1) LN1 = LN-LV-NIORB
if (LN /= LN1+LV+NIORB) then
  write(6,*) 'Input: LN.NE.LN1+LV+NIORB'
  write(6,*) 'LN,LN1,LV,NIORB=',LN,LN1,LV,NIORB
  call Quit(_RC_INPUT_ERROR_)
end if
LNP = LN*(LN+1)/2
IN = 0
IN3 = 0
IN1 = LV
IN2 = NIORB+LV
do I=1,NSYM
  NISHI = NISH(I)
  NSHI = NSH(I)
  NVALI = NVAL(I)
  NCORI = NCOR(I)
  if (NVALI == 0) GO TO 806
  do J=1,NVALI
    IN3 = IN3+1
    IN = IN+1
    NSM(IN3) = I
    ICH(IN) = IN3
  end do
806 if (NISHI == 0) GO TO 804
  do J=1,NISHI
    IN1 = IN1+1
    IN = IN+1
    NSM(IN1) = I
    ICH(IN) = IN1
    if (J > NCORI) GO TO 805
    ICOR(IN1) = 1
805 continue
  end do
804 if (NSHI == 0) GO TO 801
  do J=1,NSHI
    IN2 = IN2+1
    IN = IN+1
    NSM(IN2) = I
    ICH(IN) = IN2
  end do
801 continue
end do
if (LN > IOM) then
  write(6,*) 'Input: LN > IOM'
  write(6,*) 'LN,IOM=',LN,IOM
  call Quit(_RC_INPUT_ERROR_)
end if
call TAB2F(IVER-1,LV)
call TAB2(NREF,IOCR,nIOCR,L0,L1,L2,L3,INTNUM,LV,LSYM,ICIALL,IFCORE,ICOR,NONE,JONE)
LN2 = LN1
if (LN1 > 8) LN2 = 16
if (LN1 /= 0) GO TO 50
write(IW,55)
55 format(//,6X,'ONE CLOSED SHELL REFERENCE STATE')
GO TO 75

50 write(IW,107) (I,I=1,LN2)
107 format(//,6X,'OCCUPATION OF REFERENCE STATES',//,6X,'REF.STATE',2X,'ORB:',I2,15I4)
NO = N-2*NIORB
MAX = 0
do IREF=1,NREF
  MIN = MAX+1
  MAX = MAX+LN1
  write(IW,112) IREF,(IOCR(J),J=MIN,MAX)
112 format(6X,I5,8X,16I4)

  ! Sum up the occupation numbers of the first reference:
  ISUM = 0
  do I=1,LN1
    ISUM = ISUM+IOCR(MIN+I-1)
  end do
  if (ISUM /= NO) then
    write(6,*) ' Summed occupation nums of this reference does'
    write(6,*) ' not match nr of electrons.'
    write(6,*) ' In closed shells: 2*NIORB=',2*NIORB
    write(6,*) ' Summed occupation nums   =',NO
    write(6,*) ' Sum total is             =',2*NIORB+NO
    write(6,*) ' But input says nr of elec=',N
    call Quit(_RC_INPUT_ERROR_)
  end if
end do

75 continue
! Here with ILIM=2 (FIRST command) or 4 (normal, default).
call CONFIG(NREF,IOCR,nIOCR,L0,L1,L2,L3,JSYM,JSY,INTNUM,LSYM,JJS,ISO,LV,IFCORE,ICOR,NONE,JONE,JREFX,NFREF)
IR = JRC(ILIM)
ISPAC = IR*LNP
if (IPRINT >= 2) write(IW,9) ISPAC,ISPA
9 format(//,6X,'ELEMENTS TO BE SORTED',I7,/6X,'SORTING AREA',I16)
NRLN1 = NREF*LN1
if (LN1 == 0) NRLN1 = 1
!PAM97 IR1 = (LN*IR+29)/30
IR1 = (LN*IR+14)/15
IR2 = (IR+9)/10

iOpt = 1
nMUL = 64
nJJS = 18
nJRC = 4
call WR_GUGA(Lu_10,iOpt,IADD10,NFREF,S,N,LN,NSYM,IR1,IR2,IFIRST,INTNUM,LSYM,NREF,LN1,NRLN1,MUL,nMUL,NSH,NISH,8,JRC,nJRC,JJS,nJJS, &
             NVAL,IOCR,nIOCR)
call iDAFILE(Lu_10,1,ICASE,IR1,IADD10)
call iDAFILE(Lu_10,1,JSY,IR2,IADD10)
IAD10(2) = IADD10
call iDAFILE(Lu_10,1,JREFX,JRC(1),IADD10)

return
991 write(6,*) 'Input: End of input file encountered'
write(6,'(A,A)') 'Last Command: ',Command
call Quit(_RC_INPUT_ERROR_)
992 write(6,*) 'Input: Error while reading input!'
write(6,'(A,A)') 'Last Command: ',Command
call Quit(_RC_INPUT_ERROR_)

end subroutine INPUT_GUGA
