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
*      SUBROUTINE SORTB(BUFOUT,INDOUT,ACBDS,ACBDT,ISAB,BFACBD,NINTGR)
      SUBROUTINE SORTB(BUFS,INDS,ACBDS,ACBDT,ISAB,BFACBD,NINTGR)
      IMPLICIT REAL*8 (A-H,O-Z)
#include "SysDef.fh"
#include "warnings.fh"
#include "mrci.fh"
c      DIMENSION BUFOUT(NBSIZ2,NCHN2)
*PAM04      DIMENSION BUFOUT(*)
      DIMENSION BUFS(NBITM2,NCHN2)
c      DIMENSION INDOUT(RTOI*NBSIZ2,NCHN2)
*PAM04      DIMENSION INDOUT(*)
      DIMENSION INDS(NBITM2+2,NCHN2)
      DIMENSION ACBDS(*),ACBDT(*), BFACBD(*)
c      DIMENSION ISAB(NVIRT,NVIRT)
      DIMENSION ISAB(*)
C SORTS INTEGRALS (AB/CD)
C FOR FIXED A,C ALL B,D
      DIMENSION NORB0(9)
      CALL QENTER('SORTB')

      NVT=IROW(NVIRT+1)
      NOVST=LN*NVIRT+1
      IAD16=0

*PAM04C Buffer layout:
*PAM04      IBOFF2=RTOI*NBITM2
*PAM04      IBBC2=IBOFF2+NBITM2+1
*PAM04      IBDA2=IBBC2+1

      NORB0(1)=0
      DO I=1,NSYM
        NORB0(I+1)=NORB0(I)+NORB(I)
      END DO

      INSOUT=0
      IACMAX=0
      DO 50 ISTEP=1,IPASS
      IAD50=0
      CALL iDAFILE(LUTRA,2,iTraToc,nTraToc,IAD50)
      IDISK=0
      IACMIN=IACMAX+1
      IACMAX=IACMAX+NCHN2
      IF(IACMAX.GT.NVT)IACMAX=NVT
      IF(IACMIN.GT.IACMAX)GO TO 50

C Initialize Buffer Counts and BackChain Links.
      DO IREC=1,NCHN2
*PAM04        INDOUT(IBBC2+(IREC-1)*RTOI*NBSIZ2)=0
*PAM04        INDOUT(IBDA2+(IREC-1)*RTOI*NBSIZ2)=-1
       INDS(NBITM2+1,IREC)=0
       INDS(NBITM2+2,IREC)=-1
      END DO

C Loop over symmetry blocks of all-virtual integrals.
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

      CALL dDAFILE(LUTRA,2,TIBUF,NTIBUF,IAD50)

C Loop over index quadruples in this symm block
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

C MO integral value is made accessable at TIBUF(IOUT)
      IOUT=IOUT+1
      IF(IOUT.GT.NTIBUF) THEN
         CALL dDAFILE(LUTRA,2,TIBUF,NTIBUF,IAD50)
         IOUT=1
      END IF

C M1..M4 seqential number of orbitals.
      M1=ICH(NORB0(NSP)+NT)
      M2=ICH(NORB0(NSQ)+NU)
      M3=ICH(NORB0(NSR)+NV)
      M4=ICH(NORB0(NSS)+NX)
      IF(M1.LE.LN.OR.M2.LE.LN)GO TO 306
      IF(M3.LE.LN.OR.M4.LE.LN)GO TO 306

C Permute orbital indices to canonical order
C and put integral value in FINI
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

C Compute virtual indices.
      NA=NI-LN
      NB=NJ-LN
      NC=NK-LN
      ND=NL-LN
      ITURN=0
      IF(NA.EQ.NB.AND.NC.EQ.ND)GO TO 306
107   IAC=IROW(NA)+NC
      IF(IAC.LT.IACMIN)GO TO 106
      IF(IAC.GT.IACMAX)GO TO 106
      IF(NA.EQ.NC.AND.NB.EQ.ND)FINI=FINI/2
      NAC=IAC-IACMIN+1
*PAM04      IPOS=INDOUT(IBBC2+(NAC-1)*RTOI*NBSIZ2)+1
*PAM04      INDOUT(IBBC2+(NAC-1)*RTOI*NBSIZ2)=IPOS
      IPOS=INDS(NBITM2+1,NAC)+1
      INDS(NBITM2+1,NAC)=IPOS
*PAM04      BUFOUT(IPOS+(NAC-1)*NBSIZ2)=FINI
*PAM04      INDOUT(IBOFF2+IPOS+(NAC-1)*RTOI*NBSIZ2)=NB+2**8*ND
      INDS(IPOS,NAC)=NB+2**8*ND
      BUFS(IPOS,NAC)=FINI
      IF(IPOS.LT.NBITM2)GO TO 106
C Save this buffer if filled up.
      JDISK=IDISK
*PAM04      CALL dDAFILE(Lu_60,1,INDOUT(1+(NAC-1)*RTOI*NBSIZ2),NBSIZ2,IDISK)
      CALL iDAFILE(Lu_60,1,INDS(1,NAC),NBITM2+2,IDISK)
      CALL dDAFILE(Lu_60,1,BUFS(1,NAC),NBITM2,IDISK)
*PAM04      INDOUT(IBBC2+(NAC-1)*RTOI*NBSIZ2)=0
*PAM04      INDOUT(IBDA2+(NAC-1)*RTOI*NBSIZ2)=JDISK
      INDS(NBITM2+1,NAC)=0
      INDS(NBITM2+2,NAC)=JDISK

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
C EMPTY LAST BUFFERS
      NOVM=IACMAX-IACMIN+1
      If ( (NOVST+IACMIN-1+NOVM).gt.mChain ) then
        WRITE(6,*)'SORTB Error: NOVST+IACMIN-1+NOVM > MCHAIN'
        WRITE(6,*)'NOVST =',NOVST
        WRITE(6,*)'IACMIN=',IACMIN
        WRITE(6,*)'NOVM  =',NOVM
        WRITE(6,*)'MCHAIN=',MCHAIN
        WRITE(6,*)'  (See code).'
        CALL QUIT(_RC_GENERAL_ERROR_)
      End If
      DO 150 I=1,NOVM
      JDISK=IDISK
*PAM04      CALL dDAFILE(Lu_60,1,INDOUT(1+(I-1)*RTOI*NBSIZ2),NBSIZ2,IDISK)
      CALL iDAFILE(Lu_60,1,INDS(1,I),NBITM2+2,IDISK)
      CALL dDAFILE(Lu_60,1,BUFS(1,I),NBITM2,IDISK)
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
            ISAC=ISAB(NA+(NC-1)*NVIRT)
            NSC=NSM(LN+NC)
            NDMAX=NVIRP(NSC)+NVIR(NSC)
            IF(NDMAX.GT.NA)NDMAX=NA
            INS=ISAB(NA+(NDMAX-1)*NVIRT)
            CALL FZERO(ACBDS,INS)
            CALL FZERO(ACBDT,INS)
            IADR=LASTAD(NOVST+IAC)
201         CONTINUE
*PAM04            CALL dDAFILE(Lu_60,2,INDOUT,NBSIZ2,IADR)
            CALL iDAFILE(Lu_60,2,INDS,NBITM2+2,IADR)
            CALL dDAFILE(Lu_60,2,BUFS,NBITM2,IADR)
*PAM04            LENGTH=INDOUT(IBBC2)
*PAM04            IADR=INDOUT(IBDA2)
            LENGTH=INDS(NBITM2+1,1)
            IADR=INDS(NBITM2+2,1)
            IF(LENGTH.EQ.0)GO TO 209
            DO 202 KK=1,LENGTH
*PAM04              INND=INDOUT(IBOFF2+KK)
              INND=INDS(KK,1)
CPAM96              NB=IAND(INND,255)
CPAM96              ND=IAND(ISHFT(INND,-8),255)
*              NB=MOD(INND,2**8)
*              ND=MOD(INND/2**8,2**8)
              NB=IBITS(INND, 0,8)
              ND=IBITS(INND, 8,8)

              IBDS=ISAB(NB+(ND-1)*NVIRT)
*PAM04              ACBDS(IBDS)=ACBDS(IBDS)+BUFOUT(KK)
*PAM04              IF(NB.GT.ND)ACBDT(IBDS)=ACBDT(IBDS)+BUFOUT(KK)
*PAM04              IF(NB.LT.ND)ACBDT(IBDS)=ACBDT(IBDS)-BUFOUT(KK)
              ACBDS(IBDS)=ACBDS(IBDS)+BUFS(KK,1)
              IF(NB.GT.ND)ACBDT(IBDS)=ACBDT(IBDS)+BUFS(KK,1)
              IF(NB.LT.ND)ACBDT(IBDS)=ACBDT(IBDS)-BUFS(KK,1)
202         CONTINUE
209         IF(IADR.NE.-1) GO TO 201
            ILOOP=0
72          INSB=INS
73          INB=KBUFF1-INSOUT
            INUMB=INSB
            IF(INSB.GT.INB)INUMB=INB
            IST=INS-INSB+1
            IF(ILOOP.EQ.0)
     *      CALL DCOPY_(INUMB,ACBDS(IST),1,BFACBD(INSOUT+1),1)
            IF(ILOOP.EQ.1)
     *      CALL DCOPY_(INUMB,ACBDT(IST),1,BFACBD(INSOUT+1),1)
            INSOUT=INSOUT+INUMB
            IF(INSOUT.LT.KBUFF1)GO TO 75
            CALL dDAFILE(Lu_80,1,BFACBD,KBUFF1,IAD16)
            INSOUT=0
75          INSB=INSB-INUMB
            IF(INSB.GT.0)GO TO 73
            ILOOP=ILOOP+1
            IF(ILOOP.EQ.1)GO TO 72
60        CONTINUE
55      CONTINUE
40    CONTINUE
50    CONTINUE
C EMPTY LAST BUFFER
      IF(INSOUT.NE.0) CALL dDAFILE(Lu_80,1,BFACBD,KBUFF1,IAD16)
      CALL QEXIT('SORTB')
      RETURN
c Avoid unused argument warnings
      IF (.FALSE.) CALL Unused_integer(NINTGR)
      END
