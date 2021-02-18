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
      SUBROUTINE SORTA(BUFS,INDS,ISAB,BUFBI,BIAC,BICA,NINTGR)
      IMPLICIT REAL*8 (A-H,O-Z)
      External COUNT
#include "SysDef.fh"
#include "warnings.fh"
#include "mrci.fh"
      DIMENSION BUFS(NBITM1,NCHN1)
      DIMENSION INDS(NBITM1+2,NCHN1)
      DIMENSION BUFBI(KBUFF1)
      DIMENSION BIAC(ISMAX),BICA(ISMAX)
      DIMENSION ISAB(*)
C SORTS INTEGRALS (AB/CI)
C FOR FIXED B,I ALL A,C
C FIRST CHAIN FOR IJKL
      DIMENSION NORB0(9)
      CALL COUNT(NINTGR,NSYM,NORB,MUL)
      IF(IPRINT.GE.6) WRITE(6,1234)NINTGR
      CALL XFLUSH(6)
 1234 FORMAT(//6X,'NUMBER OF TWO-ELECTRON INTEGRALS',I10)

      IAD50=0
      CALL iDAFILE(LUTRA,2,iTraToc,nTraToc,IAD50)
*      IBOFF1=RTOI*NBITM1
*      IBBC1=IBOFF1+NBITM1+1
*      IBDA1=IBBC1+1

      IDISK=0
      ICHK=0
      DO 5 IREC=1,NCHN1
*        INDOUT(IBBC1+(IREC-1)*RTOI*NBSIZ1)=0
*        INDOUT(IBDA1+(IREC-1)*RTOI*NBSIZ1)=-1
       INDS(NBITM1+1,IREC)=0
       INDS(NBITM1+2,IREC)=-1
5     CONTINUE
      NORB0(1)=0
      DO 4 I=1,NSYM
        NORB0(I+1)=NORB0(I)+NORB(I)
4     CONTINUE
C
C TWO-ELECTRON INTEGRALS
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
          CALL dDAFILE(LUTRA,2,TIBUF,NTIBUF,IAD50)
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
                CALL dDAFILE(LUTRA,2,TIBUF,NTIBUF,IAD50)
                IOUT=1
              END IF
              FINI=TIBUF(IOUT)
              NI=ICH(NORB0(NSP)+NT)
              IF(NI.LE.0) GOTO 306
              NJ=ICH(NORB0(NSQ)+NU)
              IF(NJ.LE.0) GOTO 306
              NK=ICH(NORB0(NSR)+NV)
              IF(NK.LE.0) GOTO 306
              NL=ICH(NORB0(NSS)+NX)
              IF(NL.LE.0) GOTO 306
C ORDER THESE INDICES CANONICALLY
              IF(NI.LT.NJ) THEN
                M=NI
                NI=NJ
                NJ=M
              END IF
              IF(NK.LT.NL) THEN
                M=NK
                NK=NL
                NL=M
              END IF
              IF(NI.LT.NK) THEN
                M=NK
                NK=NI
                NI=M
                M=NL
                NL=NJ
                NJ=M
              ELSE IF((NI.EQ.NK).AND.(NJ.LT.NL)) THEN
                M=NL
                NL=NJ
                NJ=M
              END IF
              IF(NI.LE.LN)GO TO 109
              IF(NK.LE.LN)GO TO 306
              IF(IFIRST.NE.0)GO TO 306
              IF(NJ.LE.LN)GO TO 42
              IF(NL.GT.LN)GO TO 306
C ABCI
      NA=NI-LN
      NB=NJ-LN
      NC=NK-LN
      NI=NL
      GO TO 108
42    CONTINUE
      IF(NL.LE.LN) GOTO 306
C CIAB
      NA=NK-LN
      NB=NL-LN
      NC=NI-LN
      NI=NJ
108   NIB=(NI-1)*NVIRT+NB+1
*      IPOS=INDOUT(IBBC1+(NIB-1)*RTOI*NBSIZ1)+1
*      INDOUT(IBBC1+(NIB-1)*RTOI*NBSIZ1)=IPOS
      IPOS=INDS(NBITM1+1,NIB)+1
      INDS(NBITM1+1,NIB)=IPOS
*      BUFOUT(IPOS+(NIB-1)*NBSIZ1)=FINI
*      INDOUT(IBOFF1+IPOS+(NIB-1)*RTOI*NBSIZ1)=(NA-1)*NVIRT+NC
      BUFS(IPOS,NIB)=FINI
      INDS(IPOS,NIB)=(NA-1)*NVIRT+NC
      IF(IPOS.GE.NBITM1) THEN
        JDISK=IDISK
*        CALL dDAFILE(Lu_60,1,INDOUT(1+(NIB-1)*RTOI*NBSIZ1),NBSIZ1,IDISK)
        CALL iDAFILE(Lu_60,1,INDS(1,NIB),NBITM1+2,IDISK)
        CALL dDAFILE(Lu_60,1,BUFS(1,NIB),NBITM1,IDISK)
*        INDOUT(IBDA1+(NIB-1)*RTOI*NBSIZ1)=JDISK
*        INDOUT(IBBC1+(NIB-1)*RTOI*NBSIZ1)=0
        INDS(NBITM1+1,NIB)=0
        INDS(NBITM1+2,NIB)=JDISK
      END IF
      IF(NA.EQ.NB)GO TO 306
      NAT=NA
      NA=NB
      NB=NAT
      NIB=(NI-1)*NVIRT+NB+1
*      IPOS=INDOUT(IBBC1+(NIB-1)*RTOI*NBSIZ1)+1
*      INDOUT(IBBC1+(NIB-1)*RTOI*NBSIZ1)=IPOS
      IPOS=INDS(NBITM1+1,NIB)+1
      INDS(NBITM1+1,NIB)=IPOS
*      BUFOUT(IPOS+(NIB-1)*NBSIZ1)=FINI
*      INDOUT(IBOFF1+IPOS+(NIB-1)*RTOI*NBSIZ1)=(NA-1)*NVIRT+NC
      BUFS(IPOS,NIB)=FINI
      INDS(IPOS,NIB)=(NA-1)*NVIRT+NC
      IF(IPOS.GE.NBITM1) THEN
        JDISK=IDISK
*        CALL dDAFILE(Lu_60,1,INDOUT(1+(NIB-1)*RTOI*NBSIZ1),NBSIZ1,IDISK)
        CALL iDAFILE(Lu_60,1,INDS(1,NIB),NBITM1+2,IDISK)
        CALL dDAFILE(Lu_60,1,BUFS(1,NIB),NBITM1,IDISK)
*        INDOUT(IBDA1+(NIB-1)*RTOI*NBSIZ1)=JDISK
*        INDOUT(IBBC1+(NIB-1)*RTOI*NBSIZ1)=0
        INDS(NBITM1+1,NIB)=0
        INDS(NBITM1+2,NIB)=JDISK
      END IF
      GOTO 306
C IJKL
109   IIJ=IROW(NI)+NJ
      KL=IROW(NK)+NL
      IJKL=IIJ*(IIJ-1)/2+KL
      IJ=1
*      IPOS=INDOUT(IBBC1+(IJ-1)*RTOI*NBSIZ1)+1
*      INDOUT(IBBC1+(IJ-1)*RTOI*NBSIZ1)=IPOS
      IPOS=INDS(NBITM1+1,IJ)+1
      INDS(NBITM1+1,IJ)=IPOS
*      BUFOUT(IPOS+(IJ-1)*NBSIZ1)=FINI
*      INDOUT(IBOFF1+IPOS+(IJ-1)*RTOI*NBSIZ1)=IJKL
      BUFS(IPOS,IJ)=FINI
      INDS(IPOS,IJ)=IJKL
      IF(IPOS.EQ.NBITM1) THEN
        JDISK=IDISK
*        CALL dDAFILE(Lu_60,1,INDOUT(1+(IJ-1)*RTOI*NBSIZ1),NBSIZ1,IDISK)
        CALL iDAFILE(Lu_60,1,INDS(1,IJ),NBITM1+2,IDISK)
        CALL dDAFILE(Lu_60,1,BUFS(1,IJ),NBITM1,IDISK)
*        INDOUT(IBDA1+(IJ-1)*RTOI*NBSIZ1)=JDISK
*        INDOUT(IBBC1+(IJ-1)*RTOI*NBSIZ1)=0
        INDS(NBITM1+1,IJ)=0
        INDS(NBITM1+2,IJ)=JDISK
      END IF
306          CONTINUE
307         CONTINUE
308        CONTINUE
309       CONTINUE
310      CONTINUE
311     CONTINUE
312    CONTINUE
313   CONTINUE
C EMPTY LAST BUFFERS
      If ( NChn1.gt.mChain ) then
        WRITE(6,*)'SORTA Error: NCHN1 > MCHAIN (See code).'
        CALL QUIT(_RC_GENERAL_ERROR_)
      End If
      DO 150 I=1,NCHN1
        JDISK=IDISK
*        CALL dDAFILE(Lu_60,1,INDOUT(1+(I-1)*RTOI*NBSIZ1),NBSIZ1,IDISK)
        CALL iDAFILE(Lu_60,1,INDS(1,I),NBITM1+2,IDISK)
        CALL dDAFILE(Lu_60,1,BUFS(1,I),NBITM1,IDISK)
        LASTAD(I)=JDISK
150   CONTINUE
C IJKL
      IDISK=0
      IBUFIJ=0
      ISRTAD=-1
      IADR=LASTAD(1)
201   CONTINUE
*      CALL dDAFILE(Lu_60,2,INDOUT,NBSIZ1,IADR)
      CALL iDAFILE(Lu_60,2,INDS,NBITM1+2,IADR)
      CALL dDAFILE(Lu_60,2,BUFS,NBITM1,IADR)
*      LENGTH=INDOUT(IBBC1)
*      IADR=INDOUT(IBDA1)
      LENGTH=INDS(NBITM1+1,1)
      IADR=INDS(NBITM1+2,1)
      DO 202 I=1,LENGTH
        IBUFIJ=IBUFIJ+1
*        VALSRT(IBUFIJ)=BUFOUT(I)
*        INDSRT(IBUFIJ)=INDOUT(IBOFF1+I)
        VALSRT(IBUFIJ)=BUFS(I,1)
        INDSRT(IBUFIJ)=INDS(I,1)
        IF(IBUFIJ.LT.NSRTMX)GO TO 202
        NSRTCN=NSRTMX
        JDISK=IDISK
*
        INDSRT(NSRTMX+1)=NSRTCN
        INDSRT(NSRTMX+2)=ISRTAD
        CALL dDAFILE(Lu_70,1,VALSRT,NSRTMX,IDISK)
        CALL iDAFILE(Lu_70,1,INDSRT,NSRTMX+2,IDISK)
*
        ISRTAD=JDISK
        IBUFIJ=0
202   CONTINUE
      IF(IADR.NE.-1) GO TO 201
C EMPTY LAST BUFFER
      NSRTCN=IBUFIJ
      JDISK=IDISK
*
      INDSRT(NSRTMX+1)=NSRTCN
      INDSRT(NSRTMX+2)=ISRTAD
      CALL dDAFILE(Lu_70,1,VALSRT,NSRTMX,IDISK)
      CALL iDAFILE(Lu_70,1,INDSRT,NSRTMX+2,IDISK)
*
      LASTAD(1)=JDISK
C ABCI
      IAD15=IDISK
      IADABCI=IAD15
      INSOUT=0
      IADD10=IAD10(4)
      CALL dDAFILE(LUSYMB,2,COP,nCOP,IADD10)
      CALL iDAFILE(LUSYMB,2,iCOP1,nCOP+1,IADD10)
      LEN=ICOP1(nCOP+1)
      IN=2
      NSAVE=ICOP1(IN)
100   NI=NSAVE
      IOUT=0
110   IN=IN+1
      IF(IN.GT.LEN) THEN
        CALL dDAFILE(LUSYMB,2,COP,nCOP,IADD10)
        CALL iDAFILE(LUSYMB,2,iCOP1,nCOP+1,IADD10)
        LEN=ICOP1(nCOP+1)
        IF(LEN.LE.0)GO TO 6
        IN=1
      END IF
      IF(ICHK.NE.0)GO TO 460
      IF(ICOP1(IN).EQ.0) THEN
        ICHK=1
      ELSE
        IOUT=IOUT+1
      END IF
      GO TO 110
460   ICHK=0
      NSAVE=ICOP1(IN)
6     CONTINUE
      NIB=1+(NI-1)*NVIRT
C LOOP OVER VIRTUAL ORBITAL INDEX B:
      DO 20 NB=1,NVIRT
        NSIB=MUL(NSM(LN+NB),NSM(NI))
        INS=NVPAIR(NSIB)
        IF(INS.NE.0) THEN
          CALL FZERO(BIAC,INS)
          CALL FZERO(BICA,INS)
        END IF
        NIB=NIB+1
C READ & PROCESS INTEGRAL BUFFERS ON UNIT 14:
        IADR=LASTAD(NIB)
203     CONTINUE
*        CALL dDAFILE(Lu_60,2,INDOUT,NBSIZ1,IADR)
        CALL iDAFILE(Lu_60,2,INDS,NBITM1+2,IADR)
        CALL dDAFILE(Lu_60,2,BUFS,NBITM1,IADR)
*        LENGTH=INDOUT(IBBC1)
*        IADR=INDOUT(IBDA1)
        LENGTH=INDS(NBITM1+1,1)
        IADR=INDS(NBITM1+2,1)
        DO 204 KK=1,LENGTH
*          INND=INDOUT(IBOFF1+KK)
          INND=INDS(KK,1)
          NA=(INND-1)/NVIRT+1
          NC=INND-(NA-1)*NVIRT
          IACS=ISAB(NA+(NC-1)*NVIRT)
*          BIAC(IACS)=BIAC(IACS)+BUFOUT(KK)
          BIAC(IACS)=BIAC(IACS)+BUFS(KK,1)
*          IF(NA.GT.NC)BICA(IACS)=BICA(IACS)-BUFOUT(KK)
*          IF(NA.LT.NC)BICA(IACS)=BICA(IACS)+BUFOUT(KK)
          IF(NA.GT.NC)BICA(IACS)=BICA(IACS)-BUFS(KK,1)
          IF(NA.LT.NC)BICA(IACS)=BICA(IACS)+BUFS(KK,1)
204     CONTINUE
        IF(IADR.NE.-1) GO TO 203
        DO 72 ILOOP=0,1
        INSB=INS
73      INB=KBUFF1-INSOUT
        INUMB=INSB
        IF(INSB.GT.INB)INUMB=INB
        IST=INS-INSB+1
        IF(ILOOP.EQ.0)CALL DCOPY_(INUMB,BIAC(IST),1,BUFBI(INSOUT+1),1)
        IF(ILOOP.EQ.1)CALL DCOPY_(INUMB,BICA(IST),1,BUFBI(INSOUT+1),1)
        INSOUT=INSOUT+INUMB
        IF(INSOUT.GT.KBUFF1) THEN
          WRITE(6,*) 'SortA: INSOUT.GT.KBUFF1'
          WRITE(6,*) 'INSOUT=',INSOUT
          WRITE(6,*) 'KBUFF1=',KBUFF1
          CALL ABEND
        END IF
        IF(INSOUT.EQ.KBUFF1) THEN
          CALL dDAFILE(Lu_70,1,BUFBI,KBUFF1,IAD15)
          INSOUT=0
        END IF
        INSB=INSB-INUMB
        IF(INSB.GT.0)GO TO 73
72      CONTINUE
20    CONTINUE
      IF(LEN.GE.0)GO TO 100
C EMPTY LAST BUFFER IF NOT EMPTY
      IF(INSOUT.GT.0) CALL dDAFILE(Lu_70,1,BUFBI,KBUFF1,IAD15)
      RETURN
      END
