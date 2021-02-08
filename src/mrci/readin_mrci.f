************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
************************************************************************
*PAM04      SUBROUTINE READIN(HWork,iHWork)
      SUBROUTINE READIN_MRCI()
      IMPLICIT REAL*8 (A-H,O-Z)
#include "SysDef.fh"
#include "warnings.fh"
#include "mrci.fh"
#include "WrkSpc.fh"
#include "niocr.fh"
      DIMENSION IOCR(nIOCR),NOTOT(8)
*PAM04      DIMENSION HWork(*), iHWork(*)
*
      Parameter ( nCmd=20 )
      Parameter ( mxTit=10 )
      Character*4 Command,Cmd(nCmd)
      Character*72  Line,Title(mxTit)
      Character*88  ModLine
      Data Cmd /'TITL','THRP','PRIN','FROZ','DELE',
     *          'MAXI','ECON','REST','ROOT','ACPF',
     *          'SDCI','GVAL','PROR','REFC','SELE',
     *          'NRRO','MXVE','TRAN','EXTR','END '/
*
*---- convert a pointer in H to a pointer for iH
*     ipointer(i)=(i-1)*RtoI+1
*
*
*     Initialize data and set defaults
*
      IOM=MXORB
      KBUFF1=2*9600
      ETHRE=1.0D-08
      SQNLIM=1.0D-10
      CTRSH=0.05D00
      THRORB=1.0D-05
      ENP=1.0D00
      NRROOT=1
      NSEL=0
      IPRINT=1
      MAXIT=20
      IREST=0
      ICPF=0
      IREFCI=0
      ITRANS=0
      IGFAC=0
      MXVC=0
      DO 1 I=1,8
        NFRO(I)=0
        NDEL(I)=0
        NBAS(I)=0
        NORB(I)=0
1     CONTINUE
      DO 2 I=1,IOM+1
        IROW(I)=I*(I-1)/2
2     CONTINUE
      DO 4 I=1,12
        IROOT(I)=I
4     CONTINUE
      nTit=0
*
*     Read the header of the ONEINT file
*
      NAMSIZ=LENIN8*MXORB
      IDISK=0
      CALL WR_MOTRA_Info(LUONE,2,iDisk,
     &                   ITOC17,64, POTNUC,
     &                   NSYM, NBAS, NORB,NFMO,NDMO,8,NAME,NAMSIZ)
*
*---  Read input from standard input ----------------------------------*
      Call RdNLst(5,'MRCI')
10    Read(5,'(A)',End=991) Line
      Command=Line(1:4)
      Call UpCase(Command)
      If ( Command(1:1).eq.'*' ) Goto 10
      if (Command.eq.' ') Goto 10
      jCmd=0
      Do iCmd=1,nCmd
         If ( Command.eq.Cmd(iCmd) ) jCmd=iCmd
      End Do
20    Goto ( 100, 200, 300, 400, 500, 600, 700 ,800, 900,1000,
     &      1100,1200,1300,1400,1500,1600,1700,1800,1900,2000 ) jCmd
      WRITE(6,*)'READIN Error: Command not recognized.'
      WRITE(6,*)'The command is:'//''''//Command//''''
      CALL QUIT(_RC_INPUT_ERROR_)
*
*---  process TITL command --------------------------------------------*
 100  Continue
      Read(5,'(A)',End=991) Line
      Command=Line(1:4)
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
*---  process THRP command --------------------------------------------*
 200  Continue
      Read(5,'(A)',End=991) Line
      If ( Line(1:1).eq.'*' ) Goto 200
      Read(Line,*,Err=992) CTRSH
      Goto 10
*
*---  process PRIN command --------------------------------------------*
 300  Continue
      Read(5,'(A)',End=991) Line
      If ( Line(1:1).eq.'*' ) Goto 300
      Read(Line,*,Err=992) IPRINT
      Goto 10
*
*---  process FROZ command --------------------------------------------*
 400  Continue
      Read(5,'(A)',End=991) Line
      If ( Line(1:1).eq.'*' ) Goto 400
      ModLine=Line//' 0 0 0 0 0 0 0 0'
      Read(ModLine,*,Err=992) (NFRO(I),I=1,8)
      Goto 10
*
*---  process DELE command --------------------------------------------*
 500  Continue
      Read(5,'(A)',End=991) Line
      If ( Line(1:1).eq.'*' ) Goto 500
      ModLine=Line//' 0 0 0 0 0 0 0 0'
      Read(ModLine,*,Err=992) (NDEL(I),I=1,8)
      Goto 10
*
*---  process MAXI command --------------------------------------------*
 600  Continue
      Read(5,'(A)',End=991) Line
      If ( Line(1:1).eq.'*' ) Goto 600
      Read(Line,*,Err=992) MAXIT
      Goto 10
*
*---  process ECON command --------------------------------------------*
 700  Continue
      Read(5,'(A)',End=991) Line
      If ( Line(1:1).eq.'*' ) Goto 700
      Read(Line,*,Err=992) ETHRE
      Goto 10
*
*---  process REST command --------------------------------------------*
 800  Continue
      IREST=1
      Goto 10
*
*---  process ROOT command --------------------------------------------*
 900  Continue
      Read(5,'(A)',End=991) Line
      If ( Line(1:1).eq.'*' ) Goto 900
      Read(Line,*,Err=992) (IROOT(I),I=1,NRROOT)
      Goto 10
*
*---  process ACPF command --------------------------------------------*
1000  Continue
      ICPF=1
      Goto 10
*
*---  process SDCI command --------------------------------------------*
1100  Continue
      ICPF=0
      Goto 10
*
*---  process GVAL command --------------------------------------------*
1200  Continue
      Read(5,'(A)',End=991) Line
      If ( Line(1:1).eq.'*' ) Goto 1200
      Read(Line,*,Err=992) GFAC
      IGFAC=1
      Goto 10
*
*---  process PROR command --------------------------------------------*
1300  Continue
      Read(5,'(A)',End=991) Line
      If ( Line(1:1).eq.'*' ) Goto 1300
      Read(Line,*,Err=992) THRORB
      Goto 10
*
*---  process REFC command --------------------------------------------*
1400  Continue
      IREFCI=1
      Goto 10
*
*---  process SELE command --------------------------------------------*
1500  Continue
      Read(5,'(A)',End=991) Line
      If ( Line(1:1).eq.'*' ) Goto 1500
      Read(Line,*,Err=992) NSEL
      JJ=0
      Do 1510 I=1,NSEL
        Read(5,*,End=991,Err=992) NC,(CSEL(JJ+J),SSEL(JJ+J),J=1,NC)
        JJ=JJ+NC
        NCOMP(I)=NC
1510  CONTINUE
      Goto 10
*
*---  process NRRO command --------------------------------------------*
1600  Continue
      Read(5,'(A)',End=991) Line
      If ( Line(1:1).eq.'*' ) Goto 1600
      Read(Line,*,Err=992) NRROOT
      if(nrroot.gt.mxvec) then
       write(6,1610) nrroot,mxvec
1610   format('Too many roots,',i3,', max allowed is',i3)
       call quit(_RC_INPUT_ERROR_)
      endif
      DO  I=1,NRROOT
        IROOT(I)=I
      Enddo
      Goto 10
*
*---  process MXVE command --------------------------------------------*
1700  Continue
      Read(5,'(A)',End=991) Line
      If ( Line(1:1).eq.'*' ) Goto 1700
      Read(Line,*,Err=992) MXVC
       if(mxvc.gt.mxvec) then
       write(6,1710) mxvc,mxvec
1710   format('Too many vectors,',i3,', max allowed is',i3)
       call quit(_RC_INPUT_ERROR_)
      endif
      Goto 10
*
*---  process TRAN command --------------------------------------------*
1800  Continue
      ITRANS=1
      Goto 10
*
*---  process EXTR command --------------------------------------------*
1900  WRITE(6,*) 'The EXTRACT option is redundant and is ignored!'
      Goto 10
*
*---  The end of the input is reached, print the title ----------------*
2000  Continue
      if(ntit.eq.0) then
        ntit=1
        title(1)=' ( No title was given )'
      end if
      WRITE(6,*)
      CALL XFLUSH(6)
      WRITE(6,'(6X,120A1)') ('*',i=1,120)
      CALL XFLUSH(6)
      WRITE(6,'(6X,120A1)') '*',(' ',i=1,118),'*'
      CALL XFLUSH(6)
      WRITE(6,'(6X,57A1,A6,57A1)')
     &        '*',(' ',i=1,56),'Title:',(' ',i=1,56),'*'
      CALL XFLUSH(6)
      Do i=1,nTit
         Call Center(Title(i))
         WRITE(6,'(6X,24A1,A72,24A1)')
     &        '*',(' ',j=1,23),Title(i),(' ',j=1,23),'*'
      CALL XFLUSH(6)
      End Do
      WRITE(6,'(6X,120A1)') '*',(' ',i=1,118),'*'
      CALL XFLUSH(6)
      WRITE(6,'(6X,120A1)') ('*',i=1,120)
      CALL XFLUSH(6)
      WRITE(6,*)
*
*---  print the coordinates of the system -----------------------------*
      CALL XFLUSH(6)
      Call PrCoor
*
*---  read the header of CIGUGA ---------------------------------------*
*
*     Read the header of the CIGUGA file
*
      IADD10=0
      CALL iDAFILE(LUSYMB,2,IAD10,9,IADD10)
      iOpt=2
      nMUL=64
      nJJS=18
      nIRC=4
      Call WR_GUGA(LUSYMB,iOpt,IADD10,
     &             NREF,SPIN,NELEC,LN,NSYM,NCSPCK,NINTSY,IFIRST,INTNUM,
     &             LSYM,NRF,LN1,NRLN1,MUL,nMUL,NASH,NISH,8,
     &             IRC,nIRC,JJS,nJJS,NVAL,IOCR,nIOCR)
      IF(ICPF.EQ.1) THEN
        WRITE(6,*)'      THIS IS AN   A C P F   CALCULATION'
      ELSE
        WRITE(6,*)'      THIS IS AN   S D C I   CALCULATION'
        WRITE(6,*)'      (But an ACPF correction will be computed)'
      END IF
      IF(IGFAC.EQ.0) THEN
        GFAC=2.0D00/NELEC
        WRITE(6,*)'      USE THE DEFAULT ACPF G-VALUE GFAC=',GFAC
      ELSE
        WRITE(6,*)'      THE ACPF G-VALUE HAS BEEN SET TO GFAC=',GFAC
      END IF
      WRITE(6,*)
      IF(IREST.NE.0) WRITE(6,*)'      RESTARTED CALCULATION.'
      WRITE(6,*)'      A SMALL CI IS PERFORMED INVOLVING ONLY'
     &        //' THE REFERENCE STATES.'
      WRITE(6,*)'      THIS REFERENCE CI WILL USE THE FOLLOWING ROOT'
     &        //' SELECTION CRITERIA:'
      WRITE(6,*)
      CALL XFLUSH(6)
      WRITE(6,*)
      IF(MXVC.EQ.0) MXVC=MAX(NRROOT,10)
      IF(NSEL.EQ.0) THEN
        WRITE(6,*)'      ROOT SELECTION BY ENERGY ORDERING.'
        IF(NRROOT.EQ.1) THEN
          WRITE(6,'(A,I8)')'      ONE SINGLE ROOT, NUMBER ',IROOT(1)
        ELSE
          WRITE(6,*)'      THE FOLLOWING ROOTS WILL BE SELECTED:'
          WRITE(6,'(1X,/(1x,12I3))')  (IROOT(I),I=1,NRROOT)
        END IF
      ELSE
        WRITE(6,*)
     &       '      ROOT SELECTION BY PROJECTION: THE EIGENVECTORS OF'
        WRITE(6,*)
     &       '      THE REFERENCE CI ARE ORDERED BY DECREASING SIZE OF'
        WRITE(6,*)
     &       '      THEIR PROJECTIONS ONTO A SELECTION SPACE.'
        IF(NRROOT.EQ.1) THEN
          WRITE(6,*)
     &         '     SELECT THE EIGENVECTOR WITH LARGEST PROJECTION.'
        ELSE
          WRITE(6,'(A,I2,A)')'      SELECT THE ',NRROOT,
     *                       ' EIGENVECTORS WITH LARGEST PROJECTION.'
        END IF
        WRITE(6,*)
     &       '      THE SELECTION SPACE IS SPANNED BY THE FOLLOWING',
     *             'VECTORS (NONZERO COMPONENTS ONLY):'
        JJ=0
        DO 1234 I=1,NSEL
          WRITE(6,'(6X,A,I2)') ' VECTOR NR. ',I
          NC=NCOMP(I)
          WRITE(6,'(11X,I2,5X,A20,F12.8)')
     *          (J,SSEL(JJ+J),CSEL(JJ+J),J=1,NC)
          JJ=JJ+NC
1234  CONTINUE
      END IF
      WRITE(6,*)
      IF(IREFCI.EQ.0) THEN
        IF(IREST.EQ.0) THEN
      WRITE(6,*)'      THE REFERENCE CI IS FOLLOWED BY THE FULL SPACE'
      WRITE(6,*)'      CALCULATION, WHERE THE SELECTION CRITERION'
      WRITE(6,*)'      IS MAXIMUM OVERLAP WITH THE ROOT(S) SELECTED IN'
      WRITE(6,*)'      THE REFERENCE CI.'
        ELSE
      WRITE(6,*)'      THE REFERENCE CI IS FOLLOWED BY THE FULL SPACE'
      WRITE(6,*)'      CALCULATION, WITH ITERATIONS RESTARTED FROM'
      WRITE(6,*)'      CI VECTOR(S) READ FROM FILE. THE ROOT SELECTION'
      WRITE(6,*)'      CRITERION IS MAXIMUM OVERLAP WITH THE START'
      WRITE(6,*)'      VECTORS.'
        END IF
      ELSE
        WRITE(6,*)'      ONLY THE REFERENCE CI WAS REQUESTED.'
      END IF
      IF(LN.GT.IOM) THEN
        WRITE(6,*)'READIN Error: Too many orbitals.'
        WRITE(6,'(1X,A,2I5)')'actual,allowed:',LN,IOM
        CALL QUIT(_RC_INPUT_ERROR_)
      END IF
      NISHT=0
      LV=0
      DO 811 I=1,NSYM
        NISHT=NISHT+NISH(I)
        LV=LV+NVAL(I)
811   CONTINUE
      IN=0
      IR=0
      IVA=0
      IU=NISHT+LV
      IT=LV
      IV=LN
      NBAST=0
      NORBT=0
      NFMOT=0
      NFROT=0
      NASHT=0
      NVALT=0
      NVIRT=0
      NCSHT=0
      NDELT=0
      NDMOT=0
      DO 7 I=1,NSYM
        NORBI=NORB(I)
        NBASI=NBAS(I)
        NFMOI=NFMO(I)
        NFROI=NFRO(I)
        NISHI=NISH(I)
        NASHI=NASH(I)
        NVALI=NVAL(I)
        NDELI=NDEL(I)
        NDMOI=NDMO(I)
        NVIR(I)=NORBI-NFROI-NASHI-NISHI-NVALI-NDELI
        NVIRI=NVIR(I)
        NCSH(I)=NISHI+NASHI+NVALI+NVIRI
        NCSHI=NCSH(I)
        NBAST=NBAST+NBASI
        NORBT=NORBT+NORBI
        NFMOT=NFMOT+NFMOI
        NFROT=NFROT+NFROI
        NASHT=NASHT+NASHI
        NVALT=NVALT+NVALI
        NVIRT=NVIRT+NVIRI
        NCSHT=NCSHT+NCSHI
        NDELT=NDELT+NDELI
        NDMOT=NDMOT+NDMOI
        DO 8 J=1,NFROI
          IN=IN+1
          IR=IR-1
          ICH(IN)=IR
8       CONTINUE
        DO 9 J=1,NISHI
          IN=IN+1
          IT=IT+1
          ICH(IN)=IT
          NSM(IT)=I
9       CONTINUE
        DO 11 J=1,NASHI
          IN=IN+1
          IU=IU+1
          ICH(IN)=IU
          NSM(IU)=I
11      CONTINUE
        DO 12 J=1,NVALI
          IN=IN+1
          IVA=IVA+1
          ICH(IN)=IVA
          NSM(IVA)=I
12      CONTINUE
        DO 13 J=1,NVIRI
          IN=IN+1
          IV=IV+1
          ICH(IN)=IV
          NSM(IV)=I
13      CONTINUE
        DO 14 J=1,NDELI
          IN=IN+1
          ICH(IN)=0
14      CONTINUE
7     CONTINUE
      IORBS=0
      DO 721 ISYM=1,NSYM
        NOTOT(ISYM)=0
721   CONTINUE
      DO 723 ISYM=1,NSYM
        IO=NOTOT(ISYM)
        DO 722 I=1,NFMO(ISYM)+NFRO(ISYM)
          IO=IO+1
722     CONTINUE
        NOTOT(ISYM)=IO
723   CONTINUE
      DO 725 ISYM=1,NSYM
        IO=NOTOT(ISYM)
        DO 724 I=1,NISH(ISYM)
          IO=IO+1
          IORBS=IORBS+1
          IORB(IORBS)=IO
724     CONTINUE
        NOTOT(ISYM)=IO
725   CONTINUE
      DO 727 ISYM=1,NSYM
        IO=NOTOT(ISYM)
        DO 726 I=1,NASH(ISYM)
          IO=IO+1
          IORBS=IORBS+1
          IORB(IORBS)=IO
726     CONTINUE
        NOTOT(ISYM)=IO
727   CONTINUE
      DO 729 ISYM=1,NSYM
        IO=NOTOT(ISYM)
        DO 728 I=1,NVAL(ISYM)
          IO=IO+1
          IORBS=IORBS+1
          IORB(IORBS)=IO
728     CONTINUE
        NOTOT(ISYM)=IO
729   CONTINUE
      DO 731 ISYM=1,NSYM
        IO=NOTOT(ISYM)
        DO 730 I=1,NVIR(ISYM)
          IO=IO+1
          IORBS=IORBS+1
          IORB(IORBS)=IO
730     CONTINUE
        NOTOT(ISYM)=IO
731   CONTINUE
C NR OF VIRTUALS IN PREVIOUS SYMMETRIES:
      ISUM=0
      DO 732 I=1,NSYM
        NVIRP(I)=ISUM
        ISUM=ISUM+NVIR(I)
732   CONTINUE
      NCMO=0
      NBMAX=0
      DO 350 I=1,NSYM
        IF(NBAS(I).GT.NBMAX)NBMAX=NBAS(I)
        NCMO=NCMO+NBAS(I)**2
350   CONTINUE
      NBTRI=(NBAST*(NBAST+1))/2
      NVT=IROW(NVIRT+1)
      NVT2=IROW(NVIRT)
      WRITE(6,*)
      WRITE(6,'(A)')'      MALMQVIST DIAGONALIZATION'
      WRITE(6,*)
      WRITE(6,'(A,I8)')  '      PRINT LEVEL                   ',IPRINT
      WRITE(6,'(A,I12)') '      WORKSPACE WORDS, (Re*8)   '    ,MEMTOT
      WRITE(6,'(A,I8)')  '      MAXIMUM NR OF ORBITALS        ',IOM
      WRITE(6,'(A,I8)')  '      MAX NR OF STORED CI/SGM ARR.  ',MXVC
      WRITE(6,'(A,I8)')  '      MAX NR OF ITERATIONS          ',MAXIT
      WRITE(6,'(A,D9.2)')'      ENERGY CONVERGENCE THRESHOLD ' ,ETHRE
      WRITE(6,'(A,F8.1)')'      SPIN QUANTUM NUMBER           ',SPIN
      WRITE(6,'(A,I8)')  '      CORRELATED ELECTRONS          ',NELEC
      WRITE(6,'(A,I8)')  '      WAVE FUNCTION SYMMETRY LABEL  ',LSYM
      WRITE(6,'(A,I8)')  '      POINT GROUP ORDER             ',NSYM
      CALL XFLUSH(6)
      WRITE(6,*)
      WRITE(6,101)'SYMMETRY LABEL:',(I,I=1,NSYM)
      WRITE(6,101)'INACTIVE ORBITALS',(NISH(I),I=1,NSYM),NISHT
      WRITE(6,101)'ACTIVE ORBITALS',(NASH(I),I=1,NSYM),NASHT
      WRITE(6,101)'ADDED VALENCE ORB',(NVAL(I),I=1,NSYM),NVALT
      WRITE(6,101)'VIRTUAL ORBITALS',(NVIR(I),I=1,NSYM),NVIRT
      WRITE(6,*)
      WRITE(6,101)'SUM:CORREL ORBITALS',(NCSH(I),I=1,NSYM),NCSHT
      WRITE(6,*)
      WRITE(6,101)'FROZEN ORBITALS',(NFRO(I),I=1,NSYM),NFROT
      WRITE(6,101)'DELETED ORBITALS',(NDEL(I),I=1,NSYM),NDELT
      WRITE(6,*)
      WRITE(6,101)'SUM:ORBITALS IN CI',(NORB(I),I=1,NSYM),NORBT
      CALL XFLUSH(6)
      WRITE(6,*)
      WRITE(6,101)'PRE-FROZEN ORBITALS',(NFMO(I),I=1,NSYM),NFMOT
      WRITE(6,101)'PRE-DELETED ORBITALS',(NDMO(I),I=1,NSYM),NDMOT
      WRITE(6,101)'SUM:   TOTAL BASIS',(NBAS(I),I=1,NSYM),NBAST
101   FORMAT(6X,A,T47,9I5)
      WRITE(6,*)
      CALL XFLUSH(6)
      IF(LN1.EQ.0) THEN
         WRITE(6,*)'      ONE CLOSED SHELL REFERENCE STATE'
      CALL XFLUSH(6)
      ELSE
        WRITE(6,'(6X,I4,A)') NREF,' REFERENCE STATES'
        NREFWR=MIN(NREF,1000/LN1)
        LN2=MIN(32,LN1)
        WRITE(6,'(6X,A,T47)') 'Occupation of the reference states'
        IF(NREFWR.LT.NREF) THEN
          WRITE(6,'(6X,A,I3,A)')'( Only the ',NREFWR,
     &                                         ' first are listed)'
        END IF
        Write(6,'(6X,A,T25,32I2)')'Active orbital nr.',(I,I=1,LN2)
        jEnd=0
        Do iRef=1,NREFWR
          jStart=jEnd+1
          jEnd=jEnd+LN1
          Write(6,'(6X,A,I3,T25,32I2)')'Ref nr',IREF,
     &                              (IOCR(j),j=jStart,jStart-1+LN2)
        End Do
        CALL XFLUSH(6)
      END IF
      WRITE(6,*)
      CALL XFLUSH(6)
      IF(INTNUM.NE.0) WRITE(6,*)'      FIRST ORDER INTERACTING SPACE.'
      CALL XFLUSH(6)
      IX1=IRC(1)
      IX2=IRC(2)-IRC(1)
      ISC(1)=IX1
      ISC(2)=ISC(1)+IX2*NVIRT
      IY1=ISC(1)
      IY2=ISC(2)-ISC(1)
      IF(IFIRST.EQ.0) THEN
        ILIM=4
        IX3=IRC(3)-IRC(2)
        IX4=IRC(4)-IRC(3)
        ISC(3)=ISC(2)+IX3*NVT2
        ISC(4)=ISC(3)+IX4*NVT
        IY3=ISC(3)-ISC(2)
        IY4=ISC(4)-ISC(3)
        IF(IPRINT.GE.10) THEN
          WRITE(6,*)
      CALL XFLUSH(6)
          WRITE(6,*)'      INTERNAL WALKS:'
      CALL XFLUSH(6)
          WRITE(6,215)IX1,IX2,IX3,IX4
      CALL XFLUSH(6)
215       FORMAT(/,6X,'                 VALENCE',I7,
     *           /,6X,' DOUBLET COUPLED SINGLES',I7,
     *           /,6X,' TRIPLET COUPLED DOUBLES',I7,
     *           /,6X,' SINGLET COUPLED DOUBLES',I7)
          WRITE(6,*)
      CALL XFLUSH(6)
          WRITE(6,*)'      FORMAL CONFIGURATIONS:'
      CALL XFLUSH(6)
          WRITE(6,215)IY1,IY2,IY3,IY4
      CALL XFLUSH(6)
          WRITE(6,'(6X,A,I7)')'                  TOTAL:',ISC(ILIM)
      CALL XFLUSH(6)
        END IF
      ELSE
        ILIM=2
        IF(IPRINT.GE.10) THEN
          WRITE(6,*)
      CALL XFLUSH(6)
          WRITE(6,*)'      INTERNAL WALKS:'
      CALL XFLUSH(6)
          WRITE(6,216)IX1,IX2
      CALL XFLUSH(6)
216       FORMAT(/,6X,'                 VALENCE',I7,
     *           /,6X,' DOUBLET COUPLED SINGLES',I7)
          WRITE(6,*)
      CALL XFLUSH(6)
          WRITE(6,*)'      FORMAL CONFIGURATIONS:'
      CALL XFLUSH(6)
          WRITE(6,216)IY1,IY2
      CALL XFLUSH(6)
          WRITE(6,'(6X,A,I7)')'                  TOTAL:',ISC(ILIM)
      CALL XFLUSH(6)
        END IF
      END IF
      NIWLK=IRC(ILIM)
      NCVAL=IRC(1)
C ----------------------------------------------------------------
      IF (NVIRT.GT.255) THEN
        Write(6,*)
        Write(6,*)' Sorry -- The MRCI code uses internal integer codes'
        Write(6,*)' where the index of virtual orbitals is kept in'
        Write(6,*)' 8-bit fields. This cannot easily be increased'
        Write(6,*)' and limits the number of virtual orbitals to '
        Write(6,*)' 255. Your input asks for more virtuals than this.'
        Write(6,*)' The program cannot run.'
        Call Quit(_RC_INPUT_ERROR_)
      END IF
C ----------------------------------------------------------------
C ALLOCATION OF DATA PERMANENTLY IN CORE
*
*PAM04      LCSPCK=1
*
C ICSPCK - ARRAY OF BIT-PACKED GUGA CASE NUMBERS OF INTERNAL WALKS.
C CONSISTS OF NCSPCK INTEGERS.
*
*PAM04      CALL iDAFILE(LUSYMB,2,iHWork(iPointer(LCSPCK)),NCSPCK,IADD10)
      CALL GETMEM('CSPCK','ALLO','INTE',LCSPCK,NCSPCK)
      CALL iDAFILE(LUSYMB,2,IWORK(LCSPCK),NCSPCK,IADD10)
*
C INTSYM - ARRAY OF BIT-PACKED SYMMETRY LABELS OF INTERNAL WALKS.
C CONSISTS OF NINTSY INTEGERS.
*
*PAM04      LINTSY=LCSPCK+(NCSPCK+(RTOI-1))/RTOI
*PAM04      Changed to following line:
      CALL GETMEM('INTSY','ALLO','INTE',LINTSY,NINTSY)
*PAM04      CALL iDAFILE(LUSYMB,2,iHWork(iPointer(LINTSY)),NINTSY,IADD10)
      CALL iDAFILE(LUSYMB,2,IWork(LINTSY),NINTSY,IADD10)
*
C INDX - START POSITION IN CI ARRAY OF EACH INTERNAL-WALK-BLOCK
*
*PAM04      LINDX=LINTSY+(NINTSY+(RTOI-1))/RTOI
*PAM04      Changed to following line:
      CALL GETMEM('INDX','ALLO','INTE',LINDX,NIWLK)
*
C ISAB - ORDERING NR OF EACH VIRTUAL PAIR WITHIN ITS COMB-SYMM
*
*PAM04      LISAB=LINDX+NIWLK
*PAM04      Changed to following line:
      CALL GETMEM('ISAB','ALLO','INTE',LISAB,NVIRT**2)
*
C JREFX - FOR EACH VALENCE CSF, EITHER 0 OR ITS REFERENCE NR.
*
*PAM04      LJREFX=LISAB+(1+NVIRT**2)/RTOI
*PAM04      Changed to following line:
      CALL GETMEM('JREFX','ALLO','INTE',LJREFX,NCVAL)
      IADD10=IAD10(2)
      CALL iDAFILE(LUSYMB,2,iWork(LJREFX),NCVAL,IADD10)
*
C PROJECTION SELECTION VECTORS
*
*PAM04      LCISEL=LJREFX+(1+NCVAL)/RTOI
*PAM04      Changed to following line:
      CALL GETMEM('CISEL','ALLO','REAL',LCISEL,NSEL*NREF)
*
C START OF NON-PERMANENT AREA:
*
*PAM04      LPERMA=LCISEL+NSEL*NREF
*PAM04      CALL INDMAT(HWork(LCSPCK),HWork(LINTSY),HWork(LINDX),
*PAM04     *            HWork(LISAB),HWork(LJREFX),HWork(LCISEL))
      CALL INDMAT(IWork(LCSPCK),IWork(LINTSY),IWork(LINDX),
     *            IWork(LISAB),IWork(LJREFX),Work(LCISEL))
      IF(NREF.GT.MXREF) THEN
        WRITE(6,*)'READIN Error: Too many references.'
        WRITE(6,'(1X,A,2I6)')' actual, allowed:',NREF,MXREF
        CALL QUIT(_RC_INPUT_ERROR_)
      END IF

* Total available memory (at start of program) is MEMTOT
* Available now is MEMWRK
* Already (permanently) allocated is MEMPRM

      CALL GETMEM('HowMuch','MAX','REAL',LDUM,MemWrk)
      MEMPRM=MEMTOT-MEMWRK

      CALL ALLOC_MRCI
      RETURN
991   Continue
      WRITE(6,*)'READIN Error: Premature end of file while reading.'
      Call Quit(_RC_IO_ERROR_READ_)
992   Continue
      WRITE(6,*)'READIN Error: I/O error during internal read.'
      WRITE(6,*)'The line that could not be read is:'
      WRITE(6,*) Line
      Call Quit(_RC_IO_ERROR_READ_)
      END
