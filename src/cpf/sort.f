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
      SUBROUTINE SORT_CPF(BUFOUT,INDOUT,ICAD,IBUFL,FC,FIJ,FJI,TIBUF)
      IMPLICIT REAL*8 (A-H,O-Z)
#include "SysDef.fh"
#include "cpfmcpf.fh"
#include "files_cpf.fh"
      DIMENSION BUFOUT(*),INDOUT(*)
      DIMENSION ICAD(*),FC(*),IBUFL(*),FIJ(*),FJI(*),TIBUF(*)
      DIMENSION IVEC(20),IPOF(65)
      DIMENSION NORB0(9)
*
      CALL QENTER('SORT_CPF')
      IAD50=0
      CALL iDAFILE(Lu_TraInt,2,iTraToc,nTraToc,IAD50)
      NVT=IROW(NVIRT+1)
      DO 50 I=1,20
      IVEC(I)=0
50    CONTINUE
      IN=1
      DO 3 I=1,NSYM
      CALL IPO(IPOF(IN),NVIR,MUL,NSYM,I,-1)
      IN=IN+NSYM
3     CONTINUE
C ORDER OF RECORD-CHAINS IS
C 1. NOT2 CHAINS (AB/IJ)
C 2. NOT2 CHAINS (AI/BJ)
C 3. NOT2 CHAINS (AI/JK)
C RECORD STRUCTURE IS
C 1. LBUF INTEGRALS
C 2. LBUF INDICES
C 3. NUMBER OF INTEGRALS IN THIS RECORD
C 4. ADDRESS OF LAST RECORD
      NOT2=IROW(LN+1)
      NOV=3*NOT2
      NOTT=2*NOT2
      NOVST=LN*NVIRT+1+NVT
      IDISK=0
      LBUF0=RTOI*LBUF
      LBUF1=LBUF0+LBUF+1
      LBUF2=LBUF1+1
      IDIV=RTOI
      ID=0
      DO 5 IREC=1,NOV
      IBUFL(IREC)=0
      ICAD(IREC)=ID
      INDOUT(ID+LBUF2)=-1
      ID=ID+LBUF2
5     CONTINUE
      NORB0(1)=0
      DO 2 I=1,NSYM
      NORB0(I+1)=NORB0(I)+NORB(I)
2     CONTINUE
C
C     ONE ELECTRON INTEGRALS
C
      NORBTT=0
      DO 7654 ISYM=1,nsym
        NORBTT=NORBTT+(NORB(ISYM)*(NORB(ISYM)+1))/2
 7654 CONTINUE
      EMY=POTNUC
      NOB2=IROW(NORBT+1)
      IADD17=ITOC17(2)
      CALL dDAFILE(Lu_TraOne,2,FIJ,NORBTT,IADD17)
      CALL DCOPY_(NOB2,0.0D0,0,FC,1)
      IBUF=0
      KORBI=0
      DO 200 ISYM=1,NSYM
        DO 198 JORBI=KORBI+1,KORBI+NORB(ISYM)
          DO 198 IORBI=KORBI+1,JORBI
            IBUF=IBUF+1
            ONEHAM=FIJ(IBUF)
            NI=ICH(IORBI)
            NJ=ICH(JORBI)
            IF(NI.EQ.0.OR.NJ.EQ.0)GO TO 198
            IF(NI.LT.NJ) THEN
              NTMP=NI
              NI=NJ
              NJ=NTMP
            END IF
            IF(NJ.GT.0) THEN
              IJT=IROW(NI)+NJ
              FC(IJT)=FC(IJT)+ONEHAM
            ELSE IF(NI.EQ.NJ) THEN
              EMY=EMY+2.0D0*ONEHAM
            END IF
198     CONTINUE
        KORBI=KORBI+NORB(ISYM)
200   CONTINUE
      CALL DCOPY_(NOB2,0.0D0,0,FIJ,1)
      CALL DCOPY_(NOB2,0.0D0,0,FJI,1)
      IF( IPRINT.GE.20 ) THEN
         CALL TRIPRT('FC IN SORT BEFORE TWOEL',' ',FC,NORBT)
         WRITE(6,'(A,F20.8)') ' EMY:',EMY
      CALL XFLUSH(6)
      END IF
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
      IF(M1.EQ.0.OR.M2.EQ.0)GO TO 306
      IF(M3.EQ.0.OR.M4.EQ.0)GO TO 306
C     ORDER THESE INDICES CANONICALLY
      N1=M1
      N2=M2
      IF(M1.GT.M2)GO TO 11
      N1=M2
      N2=M1
11    N3=M3
      N4=M4
      IF(M3.GT.M4)GO TO 12
      N3=M4
      N4=M3
12    NI=N1
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
14    IF(NJ.GT.NL)GO TO 502
      NL=N2
      NJ=N4
502   FINI=TIBUF(IOUT)
      IF(ABS(FINI).LT.1.D-09)GO TO 306
      IF(NI.LE.0 .OR. NJ.LE.0)GO TO 41
      IF(NK.LE.0 .OR. NL.LE.0)GO TO 41
      DFINI=ABS(FINI)
      IEXP=INT(-LOG10(DFINI+1.0D-20)+5)
      IF(IEXP.GT.20) IEXP=20
      IF(IEXP.LT.1) IEXP=1
      IVEC(IEXP)=IVEC(IEXP)+1
      IF(NI.NE.NJ.OR.NK.NE.NL)GO TO 42
      IJ=IROW(NI)+NK
      FIJ(IJ)=FINI
C     SKIP (AA/II) INTEGRALS
      GO TO 306
42    IF(NI.NE.NK.OR.NJ.NE.NL)GO TO 43
      IJ=IROW(NI)+NJ
      FJI(IJ)=FINI
43    IF(NI.LE.LN)GO TO 306
      IF(NJ.GT.LN)GO TO 102
      IF(NK.GT.LN)GO TO 103
C     AIJK
      JK=NOTT+IROW(NK)+NL
      IBUFL(JK)=IBUFL(JK)+1
      ICQ=ICAD(JK)
      ICP=ICQ/IDIV+IBUFL(JK)
      BUFOUT(ICP)=FINI
      ICPP=ICQ+LBUF0+IBUFL(JK)
      INDOUT(ICPP)=IROW(NI)+NJ
      IF(IBUFL(JK).LT.LBUF)GO TO 306
      INDOUT(ICQ+LBUF1)=LBUF
      JDISK=IDISK
      CALL iDAFILE(Lu_TiABIJ,1,INDOUT(ICQ+1),LBUF2,IDISK)
      INDOUT(ICQ+LBUF2)=JDISK
      IBUFL(JK)=0
      GO TO 306
103   IF(NL.GT.LN)GO TO 306
C     AIBJ
      IIJ=NOT2+IROW(NJ)+NL
      IF(NL.GT.NJ)IIJ=NOT2+IROW(NL)+NJ
      IBUFL(IIJ)=IBUFL(IIJ)+1
      ICQ=ICAD(IIJ)
      ICP=ICQ/IDIV+IBUFL(IIJ)
      BUFOUT(ICP)=FINI
      NSA=NSM(NI)
      NAV=NI-LN-NSYS(NSA)
      NSB=NSM(NK)
      NBV=NK-LN-NSYS(NSB)
      NSIJT=(MUL(NSM(NJ),NSM(NL))-1)*NSYM
      IF(NL.GT.NJ)GO TO 105
      INAV=IPOF(NSIJT+NSA)+(NBV-1)*NVIR(NSA)+NAV
      GO TO 104
105   INAV=IPOF(NSIJT+NSB)+(NAV-1)*NVIR(NSB)+NBV
104   ICPP=ICQ+LBUF0+IBUFL(IIJ)
      INDOUT(ICPP)=INAV
      IF(IBUFL(IIJ).LT.LBUF)GO TO 108
      INDOUT(ICQ+LBUF1)=LBUF
      JDISK=IDISK
      CALL iDAFILE(Lu_TiABIJ,1,INDOUT(ICQ+1),LBUF2,IDISK)
      INDOUT(ICQ+LBUF2)=JDISK
      IBUFL(IIJ)=0
108   IF(NJ.NE.NL)GO TO 306
      IF(NI.EQ.NK)GO TO 306
      JNAV=IROW(NI)+NK
      FC(JNAV)=FC(JNAV)-FINI
      GO TO 306
102   IF(NK.GT.LN)GO TO 306
C     ABIJ
      IIJ=IROW(NK)+NL
      IBUFL(IIJ)=IBUFL(IIJ)+1
      ICQ=ICAD(IIJ)
      ICP=ICQ/IDIV+IBUFL(IIJ)
      BUFOUT(ICP)=FINI
      NSA=NSM(NI)
      NAV=NI-LN-NSYS(NSA)
      NSB=NSM(NJ)
      NBV=NJ-LN-NSYS(NSB)
      NSIJT=(MUL(NSM(NK),NSM(NL))-1)*NSYM
      INAV=IPOF(NSIJT+NSA)+(NBV-1)*NVIR(NSA)+NAV
      ICPP=ICQ+LBUF0+IBUFL(IIJ)
      INDOUT(ICPP)=INAV
      IF(IBUFL(IIJ).LT.LBUF)GO TO 106
      INDOUT(ICQ+LBUF1)=LBUF
      JDISK=IDISK
      CALL iDAFILE(Lu_TiABIJ,1,INDOUT(ICQ+1),LBUF2,IDISK)
      INDOUT(ICQ+LBUF2)=JDISK
      IBUFL(IIJ)=0
106   IF(NK.NE.NL)GO TO 306
      IF(NI.EQ.NJ)GO TO 306
      JNAV=IROW(NI)+NJ
      FC(JNAV)=FC(JNAV)+D2*FINI
      GO TO 306
C     CHECK FOR FOCK-MATRIX CONTRIBUTION
41    IF(NI.NE.NJ)GO TO 51
      II=1
      CALL IFOCK(FC,NI,NK,NL,FINI,II)
      IF(NK.NE.NL)GO TO 52
      IF(NI.GT.0.OR.NK.GT.0)GO TO 57
      EMY=EMY+D2*FINI
      IF(NI.NE.NK)EMY=EMY+D2*FINI
57    IF(NI.EQ.NK)GO TO 52
      CALL IFOCK(FC,NK,NI,NJ,FINI,II)
      GO TO 52
51    IF(NK.NE.NL)GO TO 52
      II=1
      CALL IFOCK(FC,NK,NI,NJ,FINI,II)
52    II=0
      IF(NI.NE.NK)GO TO 53
      CALL IFOCK(FC,NI,NJ,NL,FINI,II)
      IF(NJ.NE.NL)GO TO 306
      IF(NI.GT.0.OR.NJ.GT.0)GO TO 58
      EMY=EMY-FINI
      IF(NI.NE.NJ)EMY=EMY-FINI
58    IF(NI.EQ.NJ)GO TO 306
      CALL IFOCK(FC,NJ,NI,NK,FINI,II)
      GO TO 306
53    IF(NI.NE.NL)GO TO 54
      CALL IFOCK(FC,NI,NJ,NK,FINI,II)
      GO TO 306
54    IF(NJ.NE.NK)GO TO 55
      CALL IFOCK(FC,NJ,NI,NL,FINI,II)
      GO TO 306
55    IF(NJ.NE.NL)GO TO 306
      CALL IFOCK(FC,NJ,NI,NK,FINI,II)
306   CONTINUE
307   CONTINUE
308   CONTINUE
309   CONTINUE
310   CONTINUE
311   CONTINUE
312   CONTINUE
313   CONTINUE
C     EMPTY LAST BUFFERS
      If ( (NOVST+NOV).gt.mAdr ) then
         WRITE(6,*)'SORT Error: NOVST+NOV>MADR (See code).'
         CALL Abend
      End If
      DO 150 I=1,NOV
      ICQ=ICAD(I)
      INDOUT(ICQ+LBUF1)=IBUFL(I)
      JDISK=IDISK
      CALL iDAFILE(Lu_TiABIJ,1,INDOUT(ICQ+1),LBUF2,IDISK)
      LASTAD(NOVST+I)=JDISK
150   CONTINUE
      DO 95 J=1,NORBT
      IND=IROW(J+1)
      FC(IND)=FC(IND)+EMY/N
95    CONTINUE
      IADD25=0
      CALL dDAFILE(Lu_25,1,FC,NOB2,IADD25)
      IAD25S=IADD25
      WRITE(6,154)
      CALL XFLUSH(6)
154   FORMAT(//6X,'STATISTICS FOR INTEGRALS, FIRST ENTRY 10**3-10**4',
     */)
      WRITE(6,155)(IVEC(I),I=1,20)
      CALL XFLUSH(6)
155   FORMAT(6X,5I10)
      CALL QEXIT('SORT_CPF')
      RETURN
      END
