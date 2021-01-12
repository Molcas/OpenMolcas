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
*               1986, Margareta R. A. Blomberg                         *
************************************************************************
      Subroutine ReadIn_CPF(H,iH)
      Implicit Real*8 (A-H,O-Z)
C Read input and allocate memory
#include "SysDef.fh"
#include "cpfmcpf.fh"
#include "files_cpf.fh"
#include "niocr.fh"
      Dimension IOCR(nIOCR)
      Dimension H(*), iH(*)
      Logical LWSP
      Common /SPIN/ LWSP
      Parameter ( nCmd=18 )
      Parameter ( mxTit=10 )
      Character*4 Command,Cmd(nCmd)
      Character*72  Line,Title(mxTit)
      Character*88  ModLine
      Data Cmd/'TITL','MAXP','LEVS','THRP','PRIN',
     *         'FROZ','DELE','MAXI','ECON','ETRS',
     *         'REST','MCPF','CPF ','SDCI','ACPF',
     *         'LOW ','EXTR','END '/
*
*---- convert a pointer in H to a pointer for iH
      ipointer(i)=(i-1)*RtoI+1
*
*---  Initialize arrays and variables ---------------------------------*
      KBUFF1=2*9600
      D0=0.0D0
      D1=1.0D0
      D2=2.0D0
      LWSP=.FALSE.
      SQ2=SQRT(D2)
      ETHRE=1.0D-06
      CTRSH=5.0D-02
      IPRINT=5
      MAXIT=20
      IREST=0
CPAM97      IRHP=0
      ICPF=0
      ISDCI=0
      INCPF=0
      ICONV=0
      MAXITP=6
      WLEV=0.3D0
      ETOT=0.0D0
      DO I=1,8
        NPFRO(I)=0
        NFRO(I)=0
        NDEL(I)=0
        NPDEL(I)=0
        NISH(I)=0
        NASH(I)=0
        NVAL(I)=0
        NVIR(I)=0
        NORB(I)=0
        NBAS(I)=0
      END DO
      NPFROT=0
      NFROT=0
      NDELT=0
      NPDELT=0
      NISHT=0
      NASHT=0
      NVALT=0
      NVIRT=0
      NORBT=0
      NBAST=0
      DO I=1,MXORB+1
        IROW(I)=I*(I-1)/2
      END DO
      DO I=1,99
        LW(I)=0
      END DO
      nTit=0
*
*---  read the header of TRAONE ---------------------------------------*
C Note: NORB(i)=NBAS(i)-NPFRO(i)-NPDEL(i)
      NAMSIZ=LENIN8*MXORB
      IDISK=0
      CALL WR_MOTRA_Info(Lu_TraOne,2,iDisk,
     &                   ITOC17,64, POTNUC,NSYM,
     *                   NBAS,NORB,NPFRO,NPDEL,8,NAME,NAMSIZ)
*
*---  Read input from standard input ----------------------------------*
      Call RdNLst(5,'CPF')
10    Read(5,'(A)',End=991) Line
      Command=Line(1:4)
      Call UpCase(Command)
      If ( Command(1:1).eq.'*' ) Goto 10
      jCmd=0
      Do iCmd=1,nCmd
         If ( Command.eq.Cmd(iCmd) ) jCmd=iCmd
      End Do
20    Goto ( 100, 200, 300, 400, 500, 600, 700 ,800, 900,1000,
     &      1100,1200,1300,1400,1500,1600,1700,1800           ) jCmd
      WRITE(6,*)'READIN Error: Command not recognized.'
      WRITE(6,*)'The command is:'//''''//Command//''''
      CALL QUIT_OnUserError()
*
*---  process TITL command --------------------------------------------*
100   Continue
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
*---  process MAXP command --------------------------------------------*
200   Continue
      Read(5,'(A)',End=991) Line
      If ( Line(1:1).eq.'*' ) Goto 200
      Read(Line,*,Err=992) MaxItP
      Goto 10
*
*---  process LEVS command --------------------------------------------*
300   Continue
      Read(5,'(A)',End=991) Line
      If ( Line(1:1).eq.'*' ) Goto 300
      Read(Line,*,Err=992) WLev
      Goto 10
*
*---  process THRP command --------------------------------------------*
400   Continue
      Read(5,'(A)',End=991) Line
      If ( Line(1:1).eq.'*' ) Goto 400
      Read(Line,*,Err=992) CTrsh
      Goto 10
*
*---  process PRIN command --------------------------------------------*
500   Continue
      Read(5,'(A)',End=991) Line
      If ( Line(1:1).eq.'*' ) Goto 500
      Read(Line,*,Err=992) iPrint
      Goto 10
*
*---  process FROZ command --------------------------------------------*
600   Continue
      Read(5,'(A)',End=991) Line
      If ( Line(1:1).eq.'*' ) Goto 600
      ModLine=Line//' 0 0 0 0 0 0 0 0'
      Read(ModLine,*,Err=992) (nFro(i),i=1,8)
      Goto 10
*
*---  process DELE command --------------------------------------------*
700   Continue
      Read(5,'(A)',End=991) Line
      If ( Line(1:1).eq.'*' ) Goto 700
      ModLine=Line//' 0 0 0 0 0 0 0 0'
      Read(ModLine,*,Err=992) (NDEL(i),i=1,8)
      Goto 10
*
*---  process MAXI command --------------------------------------------*
800   Continue
      Read(5,'(A)',End=991) Line
      If ( Line(1:1).eq.'*' ) Goto 800
      Read(Line,*,Err=992) MaxIt
      MaxIt=Min(MaxIt,75)
      Goto 10
*
*---  process ECON command --------------------------------------------*
900   Continue
      Read(5,'(A)',End=991) Line
      If ( Line(1:1).eq.'*' ) Goto 900
      Read(Line,*,Err=992) EThre
      Goto 10
*
*---  process ETRS command --------------------------------------------*
1000  Continue
      Read(5,'(A)',End=991) Line
      If ( Line(1:1).eq.'*' ) Goto 1000
CPAM97      Read(Line,*,Err=992) ETrsh
      Read(Line,*)
      WRITE(6,*)' WARNING: The obsolete ETRS command is ignored.'
      Goto 10
*
*---  process REST command --------------------------------------------*
1100  Continue
      iRest=1
      Goto 10
*
*---  process MCPF command --------------------------------------------*
1200  Continue
      iCPF=0
      iSDCI=0
      iNCPF=0
      Goto 10
*
*---  process CPF  command --------------------------------------------*
1300  Continue
      iCPF=1
      iSDCI=0
      iNCPF=0
      Goto 10
*
*---  process SDCI command --------------------------------------------*
1400  Continue
      iSDCI=1
      iCPF=0
      iNCPF=0
      Goto 10
*
*---  process ACPF command --------------------------------------------*
1500  Continue
      iNCPF=1
      iCPF=0
      iSDCI=0
      Goto 10
*
*---  process LOW  command --------------------------------------------*
1600  Continue
      LWSP=.true.
      Goto 10
*
*---  process EXTR command --------------------------------------------*
1700  WRITE(6,*) 'The EXTRACT option is redundant and is ignored!'
      Goto 10
*
*---  The end of the input is reached, print the title ----------------*
1800  Continue
      if(ntit.eq.0) then
        ntit=1
        title(1)=' ( No title was given )'
      end if
      WRITE(6,*)
      WRITE(6,'(6X,120A1)') ('*',i=1,120)
      WRITE(6,'(6X,120A1)') '*',(' ',i=1,118),'*'
      WRITE(6,'(6X,57A1,A6,57A1)')
     &        '*',(' ',i=1,56),'Title:',(' ',i=1,56),'*'
      Do i=1,nTit
         Call Center(Title(i))
         WRITE(6,'(6X,24A1,A72,24A1)')
     &        '*',(' ',j=1,23),Title(i),(' ',j=1,23),'*'
      End Do
      WRITE(6,'(6X,120A1)') '*',(' ',i=1,118),'*'
      WRITE(6,'(6X,120A1)') ('*',i=1,120)
      WRITE(6,*)
*
*---  print the coordinates of the system -----------------------------*
      Call PrCoor
*
*---  print the method used -------------------------------------------*
       WRITE(6,*)
      If ( iSDCI.eq.1 ) then
         WRITE(6,'(6X,A)') 'This is an  S D C I  calculation'
      Else If ( iCPF.eq.1 ) then
         WRITE(6,'(6X,A)') 'This is a  C P F  calculation'
      Else If( INCPF.eq.1 ) then
         WRITE(6,'(6X,A)') 'This is an  A C P F calculation'
      Else
         WRITE(6,'(6X,A)') 'This is an  M C P F  calculation'
      End If
      If ( LWSP ) WRITE(6,'(6X,A)') 'This is a LOW SPIN calculation'
      CALL XFLUSH(6)
*
*---  read the header of CIGUGA ---------------------------------------*
      IADD10=0
      CALL iDAFILE(Lu_CIGuga,2,IAD10,9,IADD10)
      iOpt=2
      nMUL=64
      nJJS=18
      nIRC=4
      Call WR_GUGA(Lu_CIGuga,iOpt,IADD10,
     &             NFREF,S,N,LN,NSYM,IR1,IRJ,IFIRST,INTNUM,
     &             LSYM,NREF,LN1,NRLN1,MUL,nMUL,NASH,NISH,8,
     &             IRC,nIRC,JJS,nJJS,NVAL,IOCR,nIOCR)
      If ( LN.ge.MXORB ) THEN
        WRITE(6,*)'READIN Error: Too many orbitals.'
        WRITE(6,'(1X,A,2I5)')'LN,MXORB:',LN,MXORB
        CALL QUIT_OnUserError()
      END IF
      LW(1)=1
      CALL iDAFILE(Lu_CIGuga,2,iH(iPointer(LW(1))),IR1,IADD10)
      LW(2)=LW(1)+(IR1+(RTOI-1))/RTOI
      CALL iDAFILE(Lu_CIGuga,2,iH(iPointer(LW(2))),IRJ,IADD10)
*
*---  update orbital specifications -----------------------------------*
      IV0=0
      IV1=1
      IV2=2
      IV3=3
      DO 811 I=1,NSYM
        NVIR(I)=NORB(I)-NFRO(I)-NISH(I)-NASH(I)
     &                 -NVAL(I)-NDEL(I)
        NPFROT=NPFROT+NPFRO(I)
        NFROT=NFROT+NFRO(I)
        NISHT=NISHT+NISH(I)
        NASHT=NASHT+NASH(I)
        NVALT=NVALT+NVAL(I)
        NVIRT=NVIRT+NVIR(I)
        NDELT=NDELT+NDEL(I)
        NPDELT=NPDELT+NPDEL(I)
        NORBT=NORBT+NORB(I)
        NBAST=NBAST+NBAS(I)
811   CONTINUE
      IN=0
      IR=0
      IVA=0
      IU=NISHT+NVALT
      IT=NVALT
      IV=LN
      NVIRT=0
      DO 7 I=1,NSYM
        NFROI=NFRO(I)
        NISHI=NISH(I)
        NASHI=NASH(I)
        NVALI=NVAL(I)
        NDELI=NDEL(I)
CPAM97      NVIRDI=NORB(I)-NFROI-NASHI-NISHI-NVALI
CPAM97      NVIR(I)=NVIRDI-NDEL(I)
        NVIRI=NVIR(I)
        NVIRT=NVIRT+NVIRI
        DO J=1,NFROI
          IN=IN+1
          IR=IR-1
          ICH(IN)=IR
        END DO
        DO J=1,NISHI
          IN=IN+1
          IT=IT+1
          ICH(IN)=IT
          NSM(IT)=I
        END DO
        DO J=1,NASHI
          IN=IN+1
          IU=IU+1
          ICH(IN)=IU
          NSM(IU)=I
        END DO
        DO J=1,NVALI
          IN=IN+1
          IVA=IVA+1
          ICH(IN)=IVA
          NSM(IVA)=I
        END DO
        DO J=1,NVIRI
          IN=IN+1
          IV=IV+1
          ICH(IN)=IV
          NSM(IV)=I
        END DO
        DO J=1,NDELI
          IN=IN+1
          ICH(IN)=0
        END DO
 7    CONTINUE
      NVT=IROW(NVIRT+1)
      NVT2=IROW(NVIRT)
*
*---  report input specifications -------------------------------------*
      WRITE(6,*)
      WRITE(6,'(6X,A)') 'ONE-ELECTRON BASIS:'
      WRITE(6,'(6X,A)') '----------------------------'
      WRITE(6,*)
      WRITE(6,'(6X,A,T47,4X,4X,8I4)') 'Symmetry species',
     &     (iSym,iSym=1,nSym)
      WRITE(6,*)
      WRITE(6,'(6X,A,T47,I4,4x,8I4)') 'Orbitals pre-frozen in MOTRA',
     &     nPFroT,(nPFro(iSym),iSym=1,nSym)
      WRITE(6,'(6X,A,T47,I4,4x,8I4)') 'Orbitals used by this program',
     &     norbT,(nOrb(iSym),iSym=1,nSym)
      WRITE(6,'(6X,A,T47,I4,4x,8I4)') 'Pre-deleted in MOTRA',
     &     nPDelT,(nPDel(iSym),iSym=1,nSym)
      WRITE(6,'(6X,A,T47,I4,4x,8I4)') 'Sum: No. of basis functions',
     &     nBasT,(nBas(iSym),iSym=1,nSym)
      CALL XFLUSH(6)
      WRITE(6,*)
      WRITE(6,'(6X,A)') 'ORBITAL SPECIFICATION:'
      WRITE(6,'(6X,A)') '-------------------------------'
      WRITE(6,*)
      WRITE(6,'(6X,A,T47,I4,4x,8I4)') 'Orbitals frozen here',
     &     nFroT,(nFro(iSym),iSym=1,nSym)
      WRITE(6,'(6X,A,T47,I4,4x,8I4)') 'Inactive orbitals',
     &     nIShT,(nISh(iSym),iSym=1,nSym)
      WRITE(6,'(6X,A,T47,I4,4x,8I4)') 'Active orbitals',
     &     nAShT,(nASh(iSym),iSym=1,nSym)
      WRITE(6,'(6X,A,T47,I4,4x,8I4)') 'Additional valence orbitals',
     &     nValT,(nVal(iSym),iSym=1,nSym)
      WRITE(6,'(6X,A,T47,I4,4x,8I4)') 'Virtual orbitals',
     &     nVirT,(nVir(iSym),iSym=1,nSym)
      WRITE(6,'(6X,A,T47,I4,4x,8I4)') 'Orbitals deleted here',
     &     nDelT,(nDel(iSym),iSym=1,nSym)
      WRITE(6,'(6X,A,T47,I4,4x,8I4)') 'Sum: Total no. of orbitals',
     &     norbT,(nOrb(iSym),iSym=1,nSym)
      WRITE(6,*)
      CALL XFLUSH(6)
      WRITE(6,*)
      WRITE(6,'(6X,A)') 'WAVE FUNCTION SPECIFICATION:'
      WRITE(6,'(6X,A)') '----------------------------'
      WRITE(6,*)
      WRITE(6,'(6X,A,T47,I4)') 'Number of electrons in CI',N
      WRITE(6,'(6X,A,T47,I4)') 'Internal orbitals in CI',LN
      WRITE(6,'(6X,A,T47,I4)') 'External orbitals in CI',NVIRT
      WRITE(6,'(6X,A,T47,I4)') 'Number of irreps',NSYM
      WRITE(6,'(6X,A,T47,F4.1)') 'Spin quantum number',S
      WRITE(6,'(6X,A,T47,I4)') 'State symmetry',LSYM
      CALL XFLUSH(6)
      WRITE(6,*)
      WRITE(6,'(6X,A)') 'REFERENCE STATE:'
      WRITE(6,'(6X,A)') '------------------------------'
      WRITE(6,*)
      WRITE(6,'(6X,A,T47,I4)') 'Number of reference states',NREF
      LN2=MIN(16,LN1)
      If ( LN1.eq.0 ) then
         WRITE(6,'(6X,A,T47)') 'One closed shell reference state'
      Else
         WRITE(6,'(6X,A,T47)') 'Occupation of active orbitals in'
     &                       //' the reference state:'
         Write(6,'(6X,A,T25,16I4)')'Active orbital nr.',(I,I=1,LN2)
         jEnd=0
         Do iRef=1,nRef
            jStart=jEnd+1
            jEnd=jEnd+LN1
            Write(6,'(6X,A,I3,T25,16I4)')'Ref nr',IREF,
     &                                      (IOCR(j),j=jStart,jEnd)
         End Do
      End If
      WRITE(6,*)
      CALL XFLUSH(6)
      WRITE(6,'(6X,A)') 'OPTIONS:'
      WRITE(6,'(6X,A)') '--------'
      WRITE(6,*)
      WRITE(6,'(6X,A,T47,I4)') 'Print parameter',iPrint
      WRITE(6,'(6X,A,T47)') 'Pulay diagonalization'
      If ( INTNUM.ne.0 ) WRITE(6,'(6X,A,T47)')
     &                        'First order interacting space'
CPAM97      If ( IRHP.ne.0 ) WRITE(6,'(6X,A,T47)') 'Root homing'
      IX1=IRC(1)
      IX2=IRC(2)-IRC(1)
      ISC(1)=IX1
      ISC(2)=ISC(1)+IX2*NVIRT
      IY1=ISC(1)
      IY2=ISC(2)-ISC(1)
      WRITE(6,214)
      CALL XFLUSH(6)
214   FORMAT(//,6X,'INTERNAL CONFIGURATIONS')
      IF(IFIRST.NE.0)GO TO 205
      IX3=IRC(3)-IRC(2)
      IX4=IRC(4)-IRC(3)
      ISC(3)=ISC(2)+IX3*NVT2
      ISC(4)=ISC(3)+IX4*NVT
      IY3=ISC(3)-ISC(2)
      IY4=ISC(4)-ISC(3)
      WRITE(6,215)IX1,IX2,IX3,IX4
      CALL XFLUSH(6)
215   FORMAT(/,6X,'NUMBER OF VALENCE STATES',I16,
     */,6X,'NUMBER OF DOUBLET COUPLED SINGLES',I7,
     */,6X,'NUMBER OF TRIPLET COUPLED DOUBLES',I7,
     */,6X,'NUMBER OF SINGLET COUPLED DOUBLES',I7)
      WRITE(6,213)
      CALL XFLUSH(6)
213   FORMAT(//,6X,'FULL-SPACE CONFIGURATIONS (FORMAL)')
      WRITE(6,215)IY1,IY2,IY3,IY4
      CALL XFLUSH(6)
      GO TO 206
205   WRITE(6,216)IX1,IX2
      CALL XFLUSH(6)
216   FORMAT(/,6X,'NUMBER OF VALENCE STATES',I16,
     */,6X,'NUMBER OF DOUBLET COUPLED SINGLES',I7)
      WRITE(6,213)
      CALL XFLUSH(6)
      WRITE(6,216)IY1,IY2
      CALL XFLUSH(6)
206   ILIM=4
      IF(IFIRST.NE.0)ILIM=2
C ERROR CONDITIONS:
C     IF(LN.NE.NISHT+NASHT+NVALT) THEN
C       WRITE(6,*)' ERROR: Orbital specifications do not match'
C       WRITE(6,*)' input to GUGA. The number of internal'
C       WRITE(6,*)' orbitals must equal the number of inactive,'
C       WRITE(6,*)' active, and additional valence orbitals.'
C       CALL QUIT(20)
C     END IF
C ALLOCATION FOR INDEX VECTORS
C THESE VECTORS ARE PERMANENTLY IN CORE
C -- INDEX
CPAM97      LW(3)=LW(2)+IRJ
      LW(3)=LW(2)+(IRJ+(RTOI-1))/RTOI
C -- ISAB
      LW(4)=LW(3)+IRC(ILIM)
      NVIR2=NVIRT*NVIRT
C -- JREFX
      LW(5)=LW(4)+NVIR2
      IADD10=IAD10(2)
      CALL iDAFILE(Lu_CIGuga,2,iH(iPointer(LW(5))),ISC(1),IADD10)
C -- ADDRESSES NOT USED
      LW(6)=LW(5)+ISC(1)
      LW(7)=LW(6)
      LW(8)=LW(7)
      LW(9)=LW(8)
      LW(10)=LW(9)
C -- LIMIT FOR PERMANENT VECTORS
      LPERMA=LW(10)
      CALL dINDMAT(H)
      CALL ALLOC_CPF(ISMAX,LPERMA)
      RETURN
*
991   Continue
      WRITE(6,*)'READIN Error: Premature end of file while reading.'
      Call Quit_OnUserError()
992   Continue
      WRITE(6,*)'READIN Error: I/O error during internal read.'
      WRITE(6,*)'The line that could not be read is:'
      WRITE(6,*) Line
      Call Quit_OnUserError()
*
*     This is to allow type punning without an explicit interface
      CONTAINS
      SUBROUTINE dINDMAT(H)
      USE ISO_C_BINDING
      REAL*8, TARGET :: H(*)
      INTEGER, POINTER :: iH2(:),iH3(:),iH4(:),iH5(:)
      CALL C_F_POINTER(C_LOC(H(LW(2))),iH2,[1])
      CALL C_F_POINTER(C_LOC(H(LW(3))),iH3,[1])
      CALL C_F_POINTER(C_LOC(H(LW(4))),iH4,[1])
      CALL C_F_POINTER(C_LOC(H(LW(5))),iH5,[1])
      CALL INDMAT_CPF(iH2,iH3,iH4,ISMAX,iH5)
      NULLIFY(iH2,iH3,iH4,iH5)
      END SUBROUTINE dINDMAT
*
      END
