************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
*                                                                      *
* Copyright (C) 1986, Per E. M. Siegbahn                               *
************************************************************************
      SUBROUTINE INPUT_GUGA(ISO,JSYM,JSY,L0,L1,L2,L3,ISPAC)
      IMPLICIT REAL*8 (A-H,O-Z)
#include "niocr.fh"
#include "SysDef.fh"
#include "files_guga.fh"
      DIMENSION ISO(*),JSYM(*),JSY(*),L0(*),L1(*),L2(*),L3(*)
#include "real_guga.fh"
#include "integ.fh"
      COMMON/CNSTS/D0,D1,D2
#include "addr_guga.fh"
#include "warnings.fh"
      DIMENSION MLL(64),IOCR(nIOCR),JREFX(9000),
     *NISH(8),JJS(18),NVAL(8),NCOR(8),ICOR(55),IONE(8),JONE(8)
*
      Parameter ( nCmd=18 )
      Parameter ( mxTit=10 )
      Character*4 Command,Cmd(nCmd)
      Character*72  Line,Title(mxTit)
      Character*132 ModLine
      DATA MLL/1,2,3,4,5,6,7,8, 2,1,4,3,6,5,8,7, 3,4,1,2,7,8,5,6,
     * 4,3,2,1,8,7,6,5, 5,6,7,8,1,2,3,4, 6,5,8,7,2,1,4,3,
     * 7,8,5,6,3,4,1,2, 8,7,6,5,4,3,2,1/
      Data Cmd /'TITL','ELEC','SPIN','SYMM','ACTI',
     &          'PRIN','REFE','FIRS','INAC','CIAL',
     &          'VALE','INTE','NOCO','ONEO','EXTR',
     &          'NONI','NACT','END '/
*
*---  Initialize data and variables -----------------------------------*
      IOM=55
      IVER=MXVERT
      IFIRST=0
      IPRINT=0
      ICIALL=0
      ILIM=4
      N=-1
      NACTEL=-1
      IS=1
      NISHT=0
*     NSYM=1
      Call Get_iScalar('nSym',NSYM)
      NREF=0
CPAM97 New default: Interacting space.
      INTNUM=1
CPAM97 IFCORE.ne.0 means core-polarization orbitals (NOCO keyword).
      IFCORE=0
      LSYM=1
      IN=0
      DO 5 I=1,8
        NISH(I)=0
        NVAL(I)=0
        NCOR(I)=0
        NSH(I)=0
        IONE(I)=0
        DO 6 J=1,8
          IN=IN+1
          MUL(I,J)=MLL(IN)
6       CONTINUE
5     CONTINUE
      DO 207 I=1,55
        ICOR(I)=0
207   CONTINUE
      nTit=0
*
*---  Read input from standard input ----------------------------------*
      Call RdNLst(5,'GUGA')
10    Read(5,'(A)',End=991) Line
      Command=Line(1:8)
      Call UpCase(Command)
      If ( Command(1:1).eq.'*' ) Goto 10
      If (Command.eq.' ') Goto 10
      jCmd=0
      Do iCmd=1,nCmd
         If ( Command.eq.Cmd(iCmd) ) jCmd=iCmd
      End Do
20    Goto ( 100, 200, 300, 400, 500, 600, 700 ,800, 900,1000,
     &      1100,1200,1300,1400,1500,1600,1700,1800           ) jCmd
      Write (6,*) 'Input: Illegal Keyword'
      Write (6,'(A,A)') 'Command=',Command
      Call Quit(_RC_INPUT_ERROR_)
*
*---  process TITLE    command ----------------------------------------*
 100  Continue
      Read(5,'(A)',End=991) Line
      Command=Line(1:8)
      Call UpCase(Command)
      If ( Command(1:1).eq.'*' ) Goto 100
      jCmd=0
      Do iCmd=1,nCmd
         If ( Command.eq.Cmd(iCmd) ) jCmd=iCmd
      End Do
      If ( jCmd.ne.0 ) Goto 20
      nTit=nTit+1
      If ( nTit.le.mxTit ) Title(nTit)=Line
      Goto 100
*
*---  process ELECTRON command ----------------------------------------*
 200  Continue
      Read(5,'(A)',End=991) Line
      If ( Line(1:1).eq.'*' ) Goto 200
      Read(Line,*,Err=992) N
      Goto 10
*
*---  process SPIN     command ----------------------------------------*
 300  Continue
      Read(5,'(A)',End=991) Line
      If ( Line(1:1).eq.'*' ) Goto 300
      Read(Line,*,Err=992) ISPIN
      Goto 10
*
*---  process SYMMETRY command ----------------------------------------*
 400  Continue
      Write (6,*)'Input_GUGA: keyword SYMMETRY is obsolete and ignored!'
      Read(5,'(A)',End=991) Line
      If ( Line(1:1).eq.'*' ) Goto 400
      Read(Line,*,Err=992) iDummy
      Goto 10
*
*---  process ACTIVE   command ----------------------------------------*
 500  Continue
      Read(5,'(A)',End=991) Line
      If ( Line(1:1).eq.'*' ) Goto 500
      ModLine=Line//' 0 0 0 0 0 0 0 0'
      Read(ModLine,*,Err=992) (NSH(I),I=1,8)
      Goto 10
*
*---  process PRINT    command ----------------------------------------*
 600  Continue
      Read(5,'(A)',End=991) Line
      If ( Line(1:1).eq.'*' ) Goto 600
      Read(Line,*,Err=992) IPRINT
      Goto 10
*
*---  process REFERENC command ----------------------------------------*
 700  Continue
      Read(5,'(A)',End=991) Line
      If ( Line(1:1).eq.'*' ) Goto 700
      Read(Line,*,Err=992) nRef,LN1
      If ( LN1.eq.0 ) Goto 10
      jEnd=0
      Do 710 iRef=1,nRef
         jStart=jEnd+1
         jEnd=jEnd+LN1
         Read(5,'(80I1)',End=991,Err=992) (IOCR(j),j=jStart,jEnd)
 710  Continue
      Goto 10
*
*---  process FIRST    command ----------------------------------------*
 800  Continue
      IFIRST=1
      ILIM=2
      Goto 10
*
*---  process INACTIVE command ----------------------------------------*
 900  Continue
      Read(5,'(A)',End=991) Line
      If ( Line(1:1).eq.'*' ) Goto 900
      ModLine=Line//' 0 0 0 0 0 0 0 0'
      Read(ModLine,*,Err=992) (NISH(I),I=1,8)
      NISHT=0
      DO I=1,8
       NISHT=NISHT+NISH(I)
      END DO
      Goto 10
*
*---  process CIALL    command ----------------------------------------*
1000  Continue
      Read(5,'(A)',End=991) Line
      If ( Line(1:1).eq.'*' ) Goto 1000
      Read(Line,*,Err=992) LSYM
      ICIALL=1
      Goto 10
*
*---  process VALENCE  command ----------------------------------------*
1100  Continue
      Read(5,'(A)',End=991) Line
      If ( Line(1:1).eq.'*' ) Goto 1100
      ModLine=Line//' 0 0 0 0 0 0 0 0'
      Read(ModLine,*,Err=992) (NVAL(I),I=1,8)
      Goto 10
*
*---  process INTERACT command ----------------------------------------*
1200  Continue
      INTNUM=1
      Goto 10
*
*---  process NOCORR   command ----------------------------------------*
1300  Continue
      Read(5,'(A)',End=991) Line
      If ( Line(1:1).eq.'*' ) Goto 1300
      ModLine=Line//' 0 0 0 0 0 0 0 0'
      Read(ModLine,*,Err=992) (NCOR(I),I=1,8)
CPAM97 IFCORE was not set -- assume bug. Following line inserted:
      IFCORE=1
      Goto 10
*
*---  process ONEOCC   command ----------------------------------------*
1400  Continue
      Read(5,'(A)',End=991) Line
      If ( Line(1:1).eq.'*' ) Goto 1400
      ModLine=Line//' 0 0 0 0 0 0 0 0'
      Read(ModLine,*,Err=992) (IONE(I),I=1,8)
      Goto 10
*
*---  process EXTRACT  command ----------------------------------------*
1500  Write (6,*) 'Input: EXTRACT option is redundant and is ignored!'
      Goto 10
*
*---  process NON-INTERACT command ----------------------------------------*
1600  Continue
      INTNUM=0
      Goto 10
*
*---  process NACTEL       command ----------------------------------------*
1700  Continue
      Read(5,'(A)',End=991) Line
      If ( Line(1:1).eq.'*' ) Goto 200
      Read(Line,*,Err=992) NACTEL
      Goto 10
*
*---  The end of the input is reached, print the title ----------------*
1800  Continue
      if(ntit.eq.0) then
        ntit=1
        title(1)=' (No title was given)'
      end if

* Nr of correlated electrons:
      IF(N.eq.-1 .and. NACTEL.eq.-1) THEN
        Write(6,*)' Neither the number of correlated electrons '//
     &            '(Keyword ELECTRONS)'
        Write(6,*)' nor the nr of active electrons in the reference '//
     &            'space (NACTEL) has been specified.'
        Write(6,*)' The number of active electrons are set to zero.'
        NACTEL=0
      ELSE IF(N.gt.-1 .and. NACTEL.gt.-1) THEN
        Write(6,*)' Both the number of correlated electrons '//
     &            '(Keyword ELECTRONS)'
        Write(6,*)' and the nr of active electrons in the reference '//
     &            'space (NACTEL) have been specified.'
        IF(N.ne.2*NISHT+NACTEL) THEN
         N=2*NISHT+NACTEL
         Write(6,*)' Number of correlated electrons is recomputed,=',N
        END IF
      END IF
      IF(N.eq.-1) N=2*NISHT+NACTEL
      IF(NACTEL.eq.-1) NACTEL=N-2*NISHT

      Write(6,*)
      Write(6,'(6X,120A1)') ('*',i=1,120)
      Write(6,'(6X,120A1)') '*',(' ',i=1,118),'*'
      Write(6,'(6X,57A1,A6,57A1)')'*',(' ',i=1,56),'Title:',
     &                                (' ',i=1,56),'*'
      Do i=1,nTit
         Call Center(Title(i))
         Write(6,'(6X,24A1,A72,24A1)')
     &        '*',(' ',j=1,23),Title(i),(' ',j=1,23),'*'
      End Do
      Write(6,'(6X,120A1)') '*',(' ',i=1,118),'*'
      Write(6,'(6X,120A1)') ('*',i=1,120)
      Write(6,*)
      S=(ISPIN-1)*0.5D0
      IF(IFIRST.EQ.0)WRITE(IW,2)
2     FORMAT(//,6X,'ALL SINGLE AND DOUBLE REPLACEMENTS')
      IF(IFIRST.NE.0)WRITE(IW,1)
1     FORMAT(//,6X,'ONLY SINGLE REPLACEMENTS INCLUDED')
      WRITE(IW,110)N,S
110   FORMAT(/,6X,'NUMBER OF ELECTRONS IN CI',I10,
     */,6X,'TOTAL SPIN QUANTUM NUMBER',F10.2)
      WRITE(IW,109)(I,I=1,NSYM)
109   FORMAT(//,14X,'ORBITALS PER SYMMETRY',/,14X,8I5)
      WRITE(IW,106)(NISH(I),I=1,NSYM)
106   FORMAT(6X,'INACTIVE',8I5)
      WRITE(IW,108)(NSH(I),I=1,NSYM)
108   FORMAT(6X,'ACTIVE  ',8I5)
      WRITE(IW,208)(NVAL(I),I=1,NSYM)
208   FORMAT(6X,'VALENCE ',8I5)
      WRITE(IW,206)(NCOR(I),I=1,NSYM)
206   FORMAT(6X,'CORE    ',8I5)
      WRITE(IW,209)(IONE(I),I=1,NSYM)
209   FORMAT(6X,'ONEOCC  ',8I5)
      LN=0
      LV=0
      NIORB=0
      DO 803 I=1,NSYM
      NIORB=NIORB+NISH(I)
      LN=LN+NSH(I)+NISH(I)
      LV=LV+NVAL(I)
803   CONTINUE
      NONE=0
      IN=LV+NIORB
      DO 904 I=1,NSYM
      NN=IONE(I)
      IF(NN.EQ.0)GO TO 905
      DO 906 NO=1,NN
      NONE=NONE+1
      IN=IN+1
      JONE(NONE)=IN
906   CONTINUE
905   IN=IN+NSH(I)-NN
904   CONTINUE
      LN=LN+LV
      IF(ICIALL.EQ.1)LN1=LN-LV-NIORB
      If (LN.NE.LN1+LV+NIORB) Then
         Write (6,*) 'Input: LN.NE.LN1+LV+NIORB'
         Write (6,*) 'LN,LN1,LV,NIORB=',LN,LN1,LV,NIORB
         Call Quit(_RC_INPUT_ERROR_)
      End If
      LNP=LN*(LN+1)/2
      IN=0
      IN3=0
      IN1=LV
      IN2=NIORB+LV
      DO 801 I=1,NSYM
      NISHI=NISH(I)
      NSHI=NSH(I)
      NVALI=NVAL(I)
      NCORI=NCOR(I)
      IF(NVALI.EQ.0)GO TO 806
      DO 807 J=1,NVALI
      IN3=IN3+1
      IN=IN+1
      NSM(IN3)=I
      ICH(IN)=IN3
807   CONTINUE
806   IF(NISHI.EQ.0)GO TO 804
      DO 805 J=1,NISHI
      IN1=IN1+1
      IN=IN+1
      NSM(IN1)=I
      ICH(IN)=IN1
      IF(J.GT.NCORI)GO TO 805
      ICOR(IN1)=1
805   CONTINUE
804   IF(NSHI.EQ.0)GO TO 801
      DO 802 J=1,NSHI
      IN2=IN2+1
      IN=IN+1
      NSM(IN2)=I
      ICH(IN)=IN2
802   CONTINUE
801   CONTINUE
      If (LN.GT.IOM) Then
         Write (6,*) 'Input: LN.GT.IOM'
         Write (6,*) 'LN,IOM=',LN,IOM
         Call Quit(_RC_INPUT_ERROR_)
      End If
      CALL TAB2F(IVER-1,LV)
      CALL TAB2(NREF,IOCR,nIOCR,L0,L1,L2,L3,INTNUM,LV,LSYM,ICIALL,
     &          IFCORE,ICOR,NONE,JONE)
      LN2=LN1
      IF(LN1.GT.8)LN2=16
      IF(LN1.NE.0)GO TO 50
      WRITE(IW,55)
55    FORMAT(//,6X,'ONE CLOSED SHELL REFERENCE STATE')
      GO TO 75

50    WRITE(IW,107)(I,I=1,LN2)
107   FORMAT(//,6X,'OCCUPATION OF REFERENCE STATES',
     &//,6X,'REF.STATE',2X,'ORB:',I2,15I4)
      NO=N-2*NIORB
      MAX=0
      DO 111 IREF=1,NREF
       MIN=MAX+1
       MAX=MAX+LN1
       WRITE(IW,112)IREF,(IOCR(J),J=MIN,MAX)
112    FORMAT(6X,I5,8X,16I4)

* Sum up the occupation numbers of the first reference:
       ISUM=0
       DO 113 I=1,LN1
         ISUM=ISUM+IOCR(MIN+I-1)
113    CONTINUE
       If (ISUM.NE.NO) Then
          Write (6,*) ' Summed occupation nums of this reference does'
          Write (6,*) ' not match nr of electrons.'
          Write (6,*) ' In closed shells: 2*NIORB=',2*NIORB
          Write (6,*) ' Summed occupation nums   =',NO
          Write (6,*) ' Sum total is             =',2*NIORB+NO
          Write (6,*) ' But input says nr of elec=',N
          Call Quit(_RC_INPUT_ERROR_)
       End If
111   CONTINUE

75    CONTINUE
* Here with ILIM=2 (FIRST command) or 4 (normal, default).
      CALL CONFIG(NREF,IOCR,nIOCR,L0,L1,L2,L3,JSYM,JSY,INTNUM,LSYM,
     &            JJS,ISO,LV,IFCORE,ICOR,NONE,JONE,JREFX,NFREF)
      IR=JRC(ILIM)
      ISPAC=IR*LNP
      IF(IPRINT.GE.2) WRITE(IW,9)ISPAC,ISPA
9     FORMAT(//,6X,'ELEMENTS TO BE SORTED',I7,/6X,'SORTING AREA',I16)
      NRLN1=NREF*LN1
      IF(LN1.EQ.0)NRLN1=1
CPAM97      IR1=(LN*IR+29)/30
      IR1=(LN*IR+14)/15
      IR2=(IR+9)/10
*
      iOpt=1
      nMUL=64
      nJJS=18
      nJRC=4
      Call WR_GUGA(Lu_10,iOpt,IADD10,
     &             NFREF,S,N,LN,NSYM,IR1,IR2,IFIRST,INTNUM,
     &             LSYM,NREF,LN1,NRLN1,MUL,nMUL,NSH,NISH,8,
     &             JRC,nJRC,JJS,nJJS,NVAL,IOCR,nIOCR)
      CALL iDAFILE(Lu_10,1,ICASE,IR1,IADD10)
      CALL iDAFILE(Lu_10,1,JSY,IR2,IADD10)
      IAD10(2)=IADD10
      CALL iDAFILE(Lu_10,1,JREFX,JRC(1),IADD10)
*
      RETURN
991   Write (6,*) 'Input: End of input file encountered'
      Write (6,'(A,A)') 'Last Command: ',Command
      Call Quit(_RC_INPUT_ERROR_)
992   Write (6,*) 'Input: Error while reading input!'
      Write (6,'(A,A)') 'Last Command: ',Command
      Call Quit(_RC_INPUT_ERROR_)
      END
