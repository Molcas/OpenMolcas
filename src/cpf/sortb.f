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
      SUBROUTINE SORTB_CPF(BUFOUT,INDOUT,ICAD,IBUFL,TIBUF,ACBDS,ACBDT,
     *ISAB,BUFACBD)
      IMPLICIT REAL*8 (A-H,O-Z)
#include "SysDef.fh"
#include "cpfmcpf.fh"
#include "files_cpf.fh"
      DIMENSION BUFOUT(*),INDOUT(*)
      DIMENSION ICAD(*),IBUFL(*),TIBUF(NTIBUF),ACBDS(*),ACBDT(*)
      DIMENSION ISAB(*),BUFACBD(*)
      DIMENSION NORB0(9)
      PARAMETER (IPOW8=2**8)
C SORTS INTEGRALS (AB/CD) FOR FIXED A,C ALL B,D
*
      KBUFF1=2*9600
      NVT=IROW(NVIRT+1)
      NOV=(NVT-1)/IPASS+1
      NOVST=LN*NVIRT+1
      IAD16=0
      JBUF0=RTOI*JBUF
      JBUF1=JBUF0+JBUF+1
      JBUF2=JBUF1+1
      IDIV=RTOI
      NORB0(1)=0
      DO 4 I=1,NSYM
      NORB0(I+1)=NORB0(I)+NORB(I)
4     CONTINUE
      INSOUT=0
      IACMAX=0
      DO 50 ISTEP=1,IPASS
      IAD50=0
      CALL iDAFILE(Lu_TraInt,2,iTraToc,nTraToc,IAD50)
      IDISK=0
      IACMIN=IACMAX+1
      IACMAX=IACMAX+NOV
      IF(IACMAX.GT.NVT)IACMAX=NVT
      IF(IACMIN.GT.IACMAX)GO TO 50
      ID=0
      DO 5 IREC=1,NOV
      IBUFL(IREC)=0
      ICAD(IREC)=ID
      INDOUT(ID+JBUF2)=-1
      ID=ID+JBUF2
5     CONTINUE
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
      IF(M1.LE.LN.OR.M2.LE.LN)GO TO 306
      IF(M3.LE.LN.OR.M4.LE.LN)GO TO 306
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
      NA=NI-LN
      NB=NJ-LN
      NC=NK-LN
      ND=NL-LN
      ITURN=0
      IF(NA.EQ.NB.AND.NC.EQ.ND)GO TO 306
107   IAC=IROW(NA)+NC
      IF(IAC.LT.IACMIN)GO TO 106
      IF(IAC.GT.IACMAX)GO TO 106
      IF(NA.EQ.NC.AND.NB.EQ.ND)FINI=FINI/D2
      NAC=IAC-IACMIN+1
      IBUFL(NAC)=IBUFL(NAC)+1
      ICQ=ICAD(NAC)
      ICP=ICQ/IDIV+IBUFL(NAC)
      BUFOUT(ICP)=FINI
      ICPP=ICQ+JBUF0+IBUFL(NAC)
CPAM97      INDOUT(ICPP)=IOR(NB,ISHFT(ND,8))
      INDOUT(ICPP)=NB+ND*IPOW8
      IF(IBUFL(NAC).LT.JBUF)GO TO 106
      INDOUT(ICQ+JBUF1)=JBUF
      JDISK=IDISK
      CALL iDAFILE(Lu_TiABIJ,1,INDOUT(ICQ+1),JBUF2,IDISK)
      INDOUT(ICQ+JBUF2)=JDISK
      IBUFL(NAC)=0
106   IF(ITURN.EQ.1)GO TO 306
      IF(NA.EQ.NC.AND.NB.EQ.ND)GO TO 306
      IF(NA.EQ.NB.OR.NC.EQ.ND)GO TO 306
      ITURN=1
      NC=NL-LN
      ND=NK-LN
      GO TO 107
306   CONTINUE
307   CONTINUE
308   CONTINUE
309   CONTINUE
310   CONTINUE
311   CONTINUE
312   CONTINUE
313   CONTINUE
C     EMPTY LAST BUFFERS
      NOVM=IACMAX-IACMIN+1
      If ( (NOVST+IACMIN-1+NOVM).gt.mAdr ) then
        WRITE(6,*)'SORTB_CPF Error: NOVST+IACMIN-1+NOVM > MADR'
        WRITE(6,*)'  (See code).'
        CALL Abend
      End If
      DO 150 I=1,NOVM
      ICQ=ICAD(I)
      INDOUT(ICQ+JBUF1)=IBUFL(I)
      JDISK=IDISK
      CALL iDAFILE(Lu_TiABIJ,1,INDOUT(ICQ+1),JBUF2,IDISK)
      LASTAD(NOVST+IACMIN-1+I)=JDISK
150   CONTINUE
      DO 40 ISYM=1,NSYM
      IST1=IRC(3)+JJS(ISYM+9)+1
      IFIN1=IRC(3)+JJS(ISYM+10)
      INPS=IFIN1-IST1+1
      IST2=IRC(2)+JJS(ISYM)+1
      IFIN2=IRC(2)+JJS(ISYM+1)
      INPT=IFIN2-IST2+1
      ITAIL=INPS+INPT
      IF(ITAIL.EQ.0)GO TO 40
      IN1=-NVIRT
      DO 55 NA=1,NVIRT
      IN1=IN1+NVIRT
      DO 60 NC=1,NA
      IAC=IROW(NA)+NC
      IF(IAC.LT.IACMIN)GO TO 60
      IF(IAC.GT.IACMAX)GO TO 60
      IF(NA.EQ.1)GO TO 60
      NSAC=MUL(NSM(LN+NA),NSM(LN+NC))
      NSACL=MUL(NSAC,LSYM)
      IF(NSACL.NE.ISYM)GO TO 60
      NDMAX=NSYS(NSM(LN+NC)+1)
      IF(NDMAX.GT.NA)NDMAX=NA
      INS=ISAB(IN1+NDMAX)
      DO 65 I=1,INS
      ACBDS(I)=D0
      ACBDT(I)=D0
65    CONTINUE
      IADR=LASTAD(NOVST+IAC)
201   CALL iDAFILE(Lu_TiABIJ,2,INDOUT,JBUF2,IADR)
      LENGTH=INDOUT(JBUF1)
      IADR=INDOUT(JBUF2)
      IF(LENGTH.EQ.0)GO TO 209
      DO 202 KK=1,LENGTH
      INND=INDOUT(JBUF0+KK)
*      NB=MOD(INND,IPOW8)
*      ND=MOD(INND/IPOW8,IPOW8)
      NB=IBITS(INND,0,8)
      ND=IBITS(INND,8,8)
      NBD=(NB-1)*NVIRT+ND
      IBDS=ISAB(NBD)
      ACBDS(IBDS)=ACBDS(IBDS)+BUFOUT(KK)
      IF(NB.GT.ND)ACBDT(IBDS)=ACBDT(IBDS)+BUFOUT(KK)
      IF(NB.LT.ND)ACBDT(IBDS)=ACBDT(IBDS)-BUFOUT(KK)
202   CONTINUE
209   IF(IADR.NE.-1) GO TO 201
      ILOOP=0
72    DO 75 I=1,INS
      INSOUT=INSOUT+1
      IF(ILOOP.EQ.0)BUFACBD(INSOUT)=ACBDS(I)
      IF(ILOOP.EQ.1)BUFACBD(INSOUT)=ACBDT(I)
      IF(INSOUT.LT.KBUFF1)GO TO 75
      CALL dDAFILE(Lu_TiABCD,1,BUFACBD,KBUFF1,IAD16)
      INSOUT=0
75    CONTINUE
      ILOOP=ILOOP+1
      IF(ILOOP.EQ.1)GO TO 72
60    CONTINUE
55    CONTINUE
40    CONTINUE
50    CONTINUE
C     EMPTY LAST BUFFER
      IF(INSOUT.EQ.0) THEN
         RETURN
      END IF
      CALL dDAFILE(Lu_TiABCD,1,BUFACBD,KBUFF1,IAD16)
      RETURN
      END
