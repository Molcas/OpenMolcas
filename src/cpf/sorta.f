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
      SUBROUTINE SORTA_CPF(BUFOUT,INDOUT,ICAD,IBUFL,TIBUF,ISAB,BUFBI,
     *INDBI,BIAC,BICA,NINTGR)
      IMPLICIT REAL*8 (A-H,O-Z)
      EXTERNAL COUNT_CPF
#include "SysDef.fh"
#include "cpfmcpf.fh"
#include "files_cpf.fh"
      DIMENSION BUFOUT(*),INDOUT(*)
      DIMENSION ICAD(*),IBUFL(*),TIBUF(NTIBUF),ISAB(*)
      DIMENSION BUFBI(*),INDBI(*),BIAC(*),BICA(*)
C     SORTS INTEGRALS (AB/CI)
C     FOR FIXED B,I ALL A,C
C     FIRST CHAIN FOR IJKL
      DIMENSION NORB0(9)
*
      CALL QENTER('SORTA_CPF')
      CALL COUNT_CPF(NINTGR,NSYM,NORB,MUL)
      IF (IPRINT.GE.2) THEN
        WRITE(6,*)' NUMBER OF TWO-ELECTRON INTEGRALS:',NINTGR
      END IF
      IAD50=0
      CALL iDAFILE(Lu_TraInt,2,iTraToc,nTraToc,IAD50)
      KKBUF0=(RTOI*(KBUFF1+2)-2)/(RTOI+1)
      KKBUF1=RTOI*KKBUF0+KKBUF0+1
      KKBUF2=KKBUF1+1
      NOV=LN*NVIRT+1
      IF(IFIRST.NE.0)NOV=1
      IDISK=0
      KBUF0=RTOI*KBUF
      KBUF1=KBUF0+KBUF+1
      KBUF2=KBUF1+1
      IDIV=RTOI
      ID=0
      DO 5 IREC=1,NOV
        IBUFL(IREC)=0
        ICAD(IREC)=ID
        INDOUT(ID+KBUF2)=-1
        ID=ID+KBUF2
5     CONTINUE
      NORB0(1)=0
      DO 4 I=1,NSYM
        NORB0(I+1)=NORB0(I)+NORB(I)
4     CONTINUE
C
C     TWO-ELECTRON INTEGRALS
C
      DO 313 NSP=1,NSYM
        NOP=NORB(NSP)
        DO 312 NSQ=1,NSP
          NSPQ=MUL(NSP,NSQ)
          NOQ=NORB(NSQ)
          DO 311 NSR=1,NSP
            NSPQR=MUL(NSPQ,NSR)
            NOR=NORB(NSR)
            NSSM=NSR
            IF(NSR.EQ.NSP)NSSM=NSQ
            DO 310 NSS=1,NSSM
              IF(NSS.NE.NSPQR)GO TO 310
              NOS=NORB(NSS)
              NORBP=NOP*NOQ*NOR*NOS
              IF(NORBP.EQ.0)GO TO 310
              CALL dDAFILE(Lu_TraInt,2,TIBUF,NTIBUF,IAD50)
              IOUT=0
              DO 309 NV=1,NOR
                NXM=NOS
                IF(NSR.EQ.NSS)NXM=NV
                DO 308 NX=1,NXM
                  NTM=1
                  IF(NSP.EQ.NSR)NTM=NV
                  DO 307 NT=NTM,NOP
                    NUMIN=1
                    IF(NSP.EQ.NSR.AND.NT.EQ.NV)NUMIN=NX
                    NUMAX=NOQ
                    IF(NSP.EQ.NSQ)NUMAX=NT
                    DO 306 NU=NUMIN,NUMAX
                      IOUT=IOUT+1
                      IF(IOUT.GT.NTIBUF) THEN
                        CALL dDAFILE(Lu_TraInt,2,TIBUF,NTIBUF,IAD50)
                        IOUT=1
                      END IF
                      M1=ICH(NORB0(NSP)+NT)
                      M2=ICH(NORB0(NSQ)+NU)
                      M3=ICH(NORB0(NSR)+NV)
                      M4=ICH(NORB0(NSS)+NX)
                      IF(M1.LE.0.OR.M2.LE.0)GO TO 306
                      IF(M3.LE.0.OR.M4.LE.0)GO TO 306
C                     ORDER THESE INDICES CANONICALLY
                      N1=M1
                      N2=M2
                      IF(M1.GT.M2)GO TO 11
                      N1=M2
                      N2=M1
11                    N3=M3
                      N4=M4
                      IF(M3.GT.M4)GO TO 12
                      N3=M4
                      N4=M3
12                    NI=N1
                      NJ=N2
                      NK=N3
                      NL=N4
                      IF(NI.GT.NK)GO TO 502
                      IF(NI.EQ.NK)GO TO 14
                      NI=N3
                      NJ=N4
                      NK=N1
                      NL=N2
                      GO TO 502
14                    IF(NJ.GT.NL)GO TO 502
                      NL=N2
                      NJ=N4
502                   FINI=TIBUF(IOUT)
                      IF(ABS(FINI).LT.1.D-09)GO TO 306
                      IF(NI.LE.LN)GO TO 109
                      IF(NK.LE.LN)GO TO 306
                      IF(NJ.LE.LN)GO TO 42
                      IF(NL.GT.LN)GO TO 306
                      IF(IFIRST.NE.0)GO TO 306
C
C                     ABCI
C
                      NA=NI-LN
                      NB=NJ-LN
                      NC=NK-LN
                      NI=NL
108                   ITURN=0
107                   NIB=(NI-1)*NVIRT+NB+1
                      IBUFL(NIB)=IBUFL(NIB)+1
                      ICQ=ICAD(NIB)
                      ICP=ICQ/IDIV+IBUFL(NIB)
                      BUFOUT(ICP)=FINI
                      ICPP=ICQ+KBUF0+IBUFL(NIB)
                      INDOUT(ICPP)=(NA-1)*NVIRT+NC
                      IF(IBUFL(NIB).LT.KBUF)GO TO 106
                      INDOUT(ICQ+KBUF1)=KBUF
                      JDISK=IDISK
                      CALL iDAFILE(Lu_TiABIJ,1,INDOUT(ICQ+1),KBUF2,
     &                             IDISK)
                      INDOUT(ICQ+KBUF2)=JDISK
                      IBUFL(NIB)=0
106                   IF(ITURN.EQ.1.OR.NA.EQ.NB)GO TO 306
                      ITURN=1
                      NAT=NA
                      NA=NB
                      NB=NAT
                      GO TO 107
42                    IF(NJ.LE.0.OR.NL.LE.LN)GO TO 306
C
C                     CIAB
C
                      IF(IFIRST.NE.0)GO TO 306
                      NA=NK-LN
                      NB=NL-LN
                      NC=NI-LN
                      NI=NJ
                      GO TO 108
C
C                     IJKL
C
109                   IIJ=IROW(NI)+NJ
                      KL=IROW(NK)+NL
                      IJKL=IIJ*(IIJ-1)/2+KL
                      IJ=1
                      IBUFL(IJ)=IBUFL(IJ)+1
                      ICQ=ICAD(IJ)
                      ICP=ICQ/IDIV+IBUFL(IJ)
                      BUFOUT(ICP)=FINI
                      ICPP=ICQ+KBUF0+IBUFL(IJ)
                      INDOUT(ICPP)=IJKL
                      IF(IBUFL(IJ).LT.KBUF)GO TO 306
                      INDOUT(ICQ+KBUF1)=KBUF
                      JDISK=IDISK
                      CALL iDAFILE(Lu_TiABIJ,1,INDOUT(ICQ+1),KBUF2,
     &                             IDISK)
                      INDOUT(ICQ+KBUF2)=JDISK
                      IBUFL(IJ)=0
306                 CONTINUE
307               CONTINUE
308             CONTINUE
309           CONTINUE
310         CONTINUE
311       CONTINUE
312     CONTINUE
313   CONTINUE
C     EMPTY LAST BUFFERS
CFUE Start of insertion
      If ( NOV.gt.mAdr ) then
        WRITE(6,*)'SORTA_CPF Error: NOV > MADR (See code).'
        CALL Abend
      End If
CFUE End of insertion
      DO 150 I=1,NOV
        ICQ=ICAD(I)
        INDOUT(ICQ+KBUF1)=IBUFL(I)
        JDISK=IDISK
        CALL iDAFILE(Lu_TiABIJ,1,INDOUT(ICQ+1),KBUF2,IDISK)
        LASTAD(I)=JDISK
150   CONTINUE
C
C     IJKL
C
      IDISK=0
      IBUFIJ=0
      INDBI(KKBUF2)=-1
      IADR=LASTAD(1)
201   CALL iDAFILE(Lu_TiABIJ,2,INDOUT,KBUF2,IADR)
      LENGTH=INDOUT(KBUF1)
      IADR=INDOUT(KBUF2)
      IF(LENGTH.EQ.0)GO TO 209
      DO 202 I=1,LENGTH
        IBUFIJ=IBUFIJ+1
        BUFBI(IBUFIJ)=BUFOUT(I)
        INDBI(RTOI*KKBUF0+IBUFIJ)=INDOUT(KBUF0+I)
        IF(IBUFIJ.LT.KKBUF0)GO TO 202
        INDBI(KKBUF1)=KKBUF0
        JDISK=IDISK
        CALL iDAFILE(Lu_TiABCI,1,INDBI,KKBUF2,IDISK)
        INDBI(KKBUF2)=JDISK
        IBUFIJ=0
202   CONTINUE
209   IF(IADR.NE.-1) GO TO 201
C     EMPTY LAST BUFFER
      INDBI(KKBUF1)=IBUFIJ
      JDISK=IDISK
      CALL iDAFILE(Lu_TiABCI,1,INDBI,KKBUF2,IDISK)
      LASTAD(1)=JDISK
C
C     ABCI
C
      ICHK=0
      IAD15=IDISK
      IADABCI=IAD15
      INSOUT=0
      NOVST=1
      IADD10=IAD10(4)
      CALL dDAFILE(Lu_CIGuga,2,COP,nCOP,IADD10)
      CALL iDAFILE(Lu_CIGuga,2,iCOP1,nCOP+1,IADD10)
      LEN=ICOP1(nCOP+1)
      IN=2
      NSAVE=ICOP1(IN)
100   NI=NSAVE
      IOUT=0
110   IN=IN+1
      IF(IN.LE.LEN)GO TO 15
      CALL dDAFILE(Lu_CIGuga,2,COP,nCOP,IADD10)
      CALL iDAFILE(Lu_CIGuga,2,iCOP1,nCOP+1,IADD10)
      LEN=ICOP1(nCOP+1)
      IF(LEN.LE.0)GO TO 6
      IN=1
15    IF(ICHK.NE.0)GO TO 460
      IF(ICOP1(IN).EQ.0)GO TO 10
      IOUT=IOUT+1
      GO TO 110
10    ICHK=1
      GO TO 110
460   ICHK=0
      NSAVE=ICOP1(IN)
6     CONTINUE
      NIB=(NI-1)*NVIRT+NOVST
      DO 20 NB=1,NVIRT
        NSIB=MUL(NSM(LN+NB),NSM(NI))
        INS=NNS(NSIB)
        IF(INS.EQ.0)GO TO 18
        DO 21 I=1,INS
          BIAC(I)=0.0D0
          BICA(I)=0.0D0
21      CONTINUE
18      NIB=NIB+1
        IADR=LASTAD(NIB)
203     CALL iDAFILE(Lu_TiABIJ,2,INDOUT,KBUF2,IADR)
        LENGTH=INDOUT(KBUF1)
        IADR=INDOUT(KBUF2)
        IF(LENGTH.EQ.0)GO TO 210
        DO 204 KK=1,LENGTH
          INND=INDOUT(KBUF0+KK)
          NA=(INND-1)/NVIRT+1
          NC=INND-(NA-1)*NVIRT
          NAC=(NA-1)*NVIRT+NC
          IACS=ISAB(NAC)
          BIAC(IACS)=BIAC(IACS)+BUFOUT(KK)
          IF(NA.GT.NC)BICA(IACS)=BICA(IACS)-BUFOUT(KK)
          IF(NA.LT.NC)BICA(IACS)=BICA(IACS)+BUFOUT(KK)
204     CONTINUE
210     IF(IADR.NE.-1) GO TO 203
        ILOOP=0
72      DO 75 I=1,INS
          INSOUT=INSOUT+1
          IF(ILOOP.EQ.0)BUFBI(INSOUT)=BIAC(I)
          IF(ILOOP.EQ.1)BUFBI(INSOUT)=BICA(I)
          IF(INSOUT.LT.KBUFF1)GO TO 75
          CALL dDAFILE(Lu_TiABCI,1,BUFBI,KBUFF1,IAD15)
          INSOUT=0
75      CONTINUE
        ILOOP=ILOOP+1
        IF(ILOOP.EQ.1)GO TO 72
20    CONTINUE
      IF(LEN.GE.0)GO TO 100
C     EMPTY LAST BUFFER
      IF ( INSOUT.EQ.0 ) THEN
        CALL QEXIT('SORTA_CPF')
        RETURN
      END IF
      CALL dDAFILE(Lu_TiABCI,1,BUFBI,KBUFF1,IAD15)
      CALL QEXIT('SORTA_CPF')
      RETURN
      END
