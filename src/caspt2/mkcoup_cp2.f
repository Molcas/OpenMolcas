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
* Copyright (C) 1994, Per Ake Malmqvist                                *
************************************************************************
*--------------------------------------------*
* 1994  PER-AAKE MALMQUIST                   *
* DEPARTMENT OF THEORETICAL CHEMISTRY        *
* UNIVERSITY OF LUND                         *
* SWEDEN                                     *
*--------------------------------------------*
      SUBROUTINE MKCOUP_CP2(IVR,IMAW,ISGMNT,VSGMNT,NOW,IOW,
     &                  NOCP,IOCP,ILNDW,ICASE,ICOUP,
     &                  NVTAB_TMP,VTAB_TMP,NVTAB_FINAL,
     &                  ISCR,VALUE)

      IMPLICIT REAL*8 (A-H,O-Z)

#include "rasdim.fh"
#include "caspt2.fh"
#include "output.fh"
#include "pt2_guga.fh"
      COMMON /SEGTAB/ IC1(26),IC2(26),ITVPT(26),IBVPT(26),ISVC(26),
     &                NIVR,LIVR,NSGMNT,LSGMNT

C INPUT PARAMETERS:
      DIMENSION IVR(NVERT,2),IMAW(NVERT,0:3)
      DIMENSION ISGMNT(NVERT,26), VSGMNT(NVERT,26)
      DIMENSION NOW(2,NSYM,NMIDV), NOCP(MXEO,NSYM,NMIDV)
C OUTPUT PARAMETERS:
      DIMENSION IOW(2,NSYM,NMIDV), IOCP(MXEO,NSYM,NMIDV)
      DIMENSION ILNDW(NWALK),ICASE(NICASE)
      DIMENSION VTAB_TMP(NVTAB_TMP)
      INTEGER ICOUP
      DIMENSION ICOUP(3,NICOUP)
C SCRATCH PARAMETERS:
      DIMENSION ISCR(7,0:NLEV), VALUE(0:NLEV)
      PARAMETER (IVLFT=1,ITYPE=2,IAWSL=3,IAWSR=4,ILS=5,ICS=6)
      PARAMETER (ISEG=7)

C NOW IS ZEROED AND WILL BE USED AS AN ARRAY OF COUNTERS, BUT WILL
C    BE RESTORED FINALLY.
      DO IHALF=1,2
        DO MV=1,NMIDV
          DO IS=1,NSYM
            NOW(IHALF,IS,MV)=0
          END DO
        END DO
      END DO
C SIMILAR FOR THE COUPLING COEFFICIENT TABLE:
      DO INDEO=1,MXEO
        DO MV=1,NMIDV
          DO IS=1,NSYM
            NOCP(INDEO,IS,MV)=0
          END DO
        END DO
      END DO

C COUPLING COEFFICIENT VALUE TABLE:
      NVTAB_FINAL=2
      VTAB_TMP(1)=1.0D00
      VTAB_TMP(2)=-1.0D00

      NCHECK=0

      DO IHALF=1,2
        IF(IHALF.EQ.1) THEN
          IVTSTA=1
          IVTEND=1
          LEV1=NLEV
          LEV2=MIDLEV
          ITYPMX=0
        ELSE
          IVTSTA=MIDV1
          IVTEND=MIDV2
          LEV1=MIDLEV
          LEV2=0
          ITYPMX=2
        END IF
        DO IVTOP=IVTSTA,IVTEND
        DO ITYP=0,ITYPMX
          IVRTOP=IVTOP
          IF(ITYP.GT.0)IVRTOP=IVR(IVTOP,ITYP)
          IF(IVRTOP.EQ.0) GOTO 400
          LEV=LEV1
          ISCR(IVLFT,LEV)=IVTOP
          ISCR(ITYPE,LEV)=ITYP
          ISCR(IAWSL,LEV)=0
          ISCR(IAWSR,LEV)=0
          ISCR(ILS,LEV)=1
          ISCR(ISEG,LEV)=0
          VALUE(LEV)=1.0D00
 100      IF(LEV.GT.LEV1) GOTO 400
          ITYPT=ISCR(ITYPE,LEV)
          IVLT=ISCR(IVLFT,LEV)
          DO ISGT=ISCR(ISEG,LEV)+1,26
            IVLB=ISGMNT(IVLT,ISGT)
            IF(IVLB.EQ.0) GOTO 110
            IF(ITYPT.EQ.ITVPT(ISGT)) GOTO 200
 110        CONTINUE
          END DO
          ISCR(ISEG,LEV)=0
          LEV=LEV+1
          GOTO 100

 200      ISCR(ISEG,LEV)=ISGT
          ICL=IC1(ISGT)
          ISYM=1
          IF((ICL.EQ.1).OR.(ICL.EQ.2)) ISYM=ISM(LEV)
          IVRT=IVLT
          IF((ITYPT.EQ.1).OR.(ITYPT.EQ.2)) IVRT=IVR(IVLT,ITYPT)
          ICR=IC2(ISGT)
          ISCR(ICS,LEV)=ICL
          LEV=LEV-1
          ISCR(IAWSL,LEV)=ISCR(IAWSL,LEV+1)+IMAW(IVLT,ICL)
          ISCR(IAWSR,LEV)=ISCR(IAWSR,LEV+1)+IMAW(IVRT,ICR)
          VALUE(LEV)=VALUE(LEV+1)*VSGMNT(IVLT,ISGT)
          ISCR(ILS,LEV)=MUL(ISYM,ISCR(ILS,LEV+1))
          ISCR(IVLFT,LEV)=IVLB
          ISCR(ITYPE,LEV)=IBVPT(ISGT)
          ISCR(ISEG,LEV)=0
          IF (LEV.GT.LEV2) GOTO 100

          MV=ISCR(IVLFT,MIDLEV)+1-MIDV1
          LFTSYM=ISCR(ILS,LEV2)
          IT=ISCR(ITYPE,MIDLEV)
          IF(IT.EQ.0) IT=3
          IF(ISCR(ITYPE,LEV2).EQ.0) IT=0

          IF(IT.EQ.0) THEN
            ILND=1+NOW(IHALF,LFTSYM,MV)
            IAWS=ISCR(IAWSL,LEV2)
            ILNDW(IAWS)=ILND
            NOW(IHALF,LFTSYM,MV)=ILND
            IPOS=IOW(IHALF,LFTSYM,MV)+(ILND-1)*NIPWLK
            DO LL=LEV2+1,LEV1,15
              IC=0
              DO L=MIN(LL+14,LEV1),LL,-1
                IC=4*IC+ISCR(ICS,L)
              END DO
              IPOS=IPOS+1
              ICASE(IPOS)=IC
            END DO
          ELSE
            IP=0
            IQ=0
            DO L=LEV2+1,LEV1
              ISG=ISCR(ISEG,L)
              IF((ISG.GE.5).AND.(ISG.LE.8))IP=L
              IF((ISG.GE.19).AND.(ISG.LE.22))IQ=L
            END DO
            IF(IP.EQ.0) IP=IQ
            INDEO=NLEV*(IT-1)+IP
            IF(IT.EQ.3) INDEO=2*NLEV+(IP*(IP-1))/2+IQ
            ICOP=1+NOCP(INDEO,LFTSYM,MV)
            NOCP(INDEO,LFTSYM,MV)=ICOP
            ICOP=IOCP(INDEO,LFTSYM,MV)+ICOP
            NCHECK=NCHECK+1
            IF (ICOP.GT.NICOUP) THEN
              WRITE(6,*)' ERROR: NICOUP=',NICOUP
              WRITE(6,*)' NR OF COUPS PRODUCED:',NCHECK
              WRITE(6,*)'           TYPE NR IT:',IT
              WRITE(6,*)'            IP,IQ    :',IP,IQ
              WRITE(6,*)'            INDEO    :',INDEO
              WRITE(6,*)'        MIDVERTEX MV :',MV
              WRITE(6,*)' LEFT SYMMETRY LFTSYM:',LFTSYM
              WRITE(6,*)' COUP OFFSET IOCP    :',IOCP(INDEO,LFTSYM,MV)
              WRITE(6,*)' COUP SERIAL NR ICOP :',ICOP
              WRITE(6,*)' D:O, WITHOUT OFFSET :',
     &                    ICOP-IOCP(INDEO,LFTSYM,MV)
              WRITE(6,*)' CURRENT NOCP NUMBER :',NOCP(INDEO,LFTSYM,MV)
              CALL ABEND()
            END IF
C
            C=VALUE(LEV2)
            DO I=1,NVTAB_FINAL
              IVTAB=I
              IF(ABS(C-VTAB_TMP(I)).LT.1.0D-10) GOTO 212
            END DO
            NVTAB_FINAL=NVTAB_FINAL+1
            IF(NVTAB_FINAL.GT.NVTAB_TMP) THEN
              WRITE(6,*)'MKCOUP: NVTAB_FINAL=',NVTAB_FINAL
              WRITE(6,*)'NVTAB_FINAL should not be allowed to grow'
              WRITE(6,*)'beyond NVTAB_TMP which was set provisionally'
              WRITE(6,*)'in subroutine GINIT in file ginit.f.'
              WRITE(6,*)'Now NVTAB_TMP=',NVTAB_TMP
              WRITE(6,*)'This may indicate a problem with your input.'
              WRITE(6,*)'If you do want to do this big calculation, try'
              WRITE(6,*)'increasing NVTAB_TMP in GINIT and recompile.'
              CALL ABEND()
            END IF
            VTAB_TMP(NVTAB_FINAL)=C
            IVTAB=NVTAB_FINAL
 212        ICOUP(1,ICOP)=ISCR(IAWSL,LEV2)
            ICOUP(2,ICOP)=ISCR(IAWSR,LEV2)
            ICOUP(3,ICOP)=IVTAB
            IF (ICOP.GT.NICOUP) THEN
              WRITE(6,*)'MKCOUP: ICOP>NICOUP!'
              CALL ABEND()
            END IF
          END IF

          LEV=LEV+1
          GOTO 100

 400      CONTINUE
          END DO
        END DO
      END DO
C RENUMBER THE COUPLING COEFFICIENT INDICES BY LUND SCHEME:
      DO ICOP=1,NICOUP
        I1=ICOUP(1,ICOP)
        I2=ICOUP(2,ICOP)
        ICOUP(1,ICOP)=ILNDW(I1)
        ICOUP(2,ICOP)=ILNDW(I2)
      END DO

#ifdef _DEBUGPRINT_
        ICOP1=0
        ICOP2=0
        WRITE(6,*)' NR OF DIFFERENT VALUES OF COUP:',NVTAB_FINAL
        DO ICOP=1,NICOUP
          I3=ICOUP(3,ICOP)
          IF(I3.EQ.1) ICOP1=ICOP1+1
          IF(I3.EQ.2) ICOP2=ICOP2+1
        END DO
        WRITE(6,*)
        WRITE(6,*)' NR OF COUPS WITH VALUE  1.0:',ICOP1
        WRITE(6,*)' NR OF COUPS WITH VALUE -1.0:',ICOP2
        WRITE(6,*)
        WRITE(6,*)' COUPLING COEFFICIENTS:'
        WRITE(6,*)'    IP    IQ    MV LFTSYM NOCP'
        WRITE(6,*)' 1. OPEN LOOPS TYPE 1.'
        DO IP=1,NLEV
          DO MV=1,NMIDV
            DO LFTSYM=1,NSYM
              N=NOCP(IP,LFTSYM,MV)
              WRITE(6,'(6X,6(1X,I5),F10.7)') IP,MV,LFTSYM,N
            END DO
          END DO
        END DO
        WRITE(6,*)' 2. OPEN LOOPS TYPE 2.'
        DO IP=1,NLEV
          INDEO=NLEV+IP
          DO MV=1,NMIDV
            DO LFTSYM=1,NSYM
              N=NOCP(INDEO,LFTSYM,MV)
              WRITE(6,'(6X,6(1X,I5),F10.7)') IP,MV,LFTSYM,N
            END DO
          END DO
        END DO
        WRITE(6,*)' 3. CLOSED LOOPS.'
        DO IP=1,NLEV
         DO IQ=1,IP
          INDEO=2*NLEV+(IP*(IP-1))/2+IQ
          DO MV=1,NMIDV
            DO LFTSYM=1,NSYM
              N=NOCP(INDEO,LFTSYM,MV)
              WRITE(6,'(7(1X,I5),F10.7)') IP,IQ,MV,LFTSYM,N
            END DO
          END DO
         END DO
        END DO
        WRITE(6,*)
        WRITE(6,*)' COUPLING COEFFICIENTS:'
        WRITE(6,*)'    IP    IQ    MV LFTSYM ICOP ICOUP1&2   COUP'
        WRITE(6,*)' 1. OPEN LOOPS TYPE 1.'
        DO IP=1,NLEV
          DO MV=1,NMIDV
            DO LFTSYM=1,NSYM
              N=NOCP(IP,LFTSYM,MV)
              ICOP=IOCP(IP,LFTSYM,MV)
              DO I=1,N
                ICOP=ICOP+1
                ICP1=ICOUP(1,ICOP)
                ICP2=ICOUP(2,ICOP)
                CP=VTAB_TMP(ICOUP(3,ICOP))
                WRITE(6,'(6X,6(1X,I5),F10.7)') IP,MV,LFTSYM,
     &                                      ICOP,ICP1,ICP2,CP
            END DO
          END DO
         END DO
        END DO
        WRITE(6,*)' 2. OPEN LOOPS TYPE 2.'
        DO IP=1,NLEV
          INDEO=NLEV+IP
          DO MV=1,NMIDV
            DO LFTSYM=1,NSYM
              N=NOCP(INDEO,LFTSYM,MV)
              ICOP=IOCP(INDEO,LFTSYM,MV)
              DO I=1,N
                ICOP=ICOP+1
                ICP1=ICOUP(1,ICOP)
                ICP2=ICOUP(2,ICOP)
                CP=VTAB_TMP(ICOUP(3,ICOP))
                WRITE(6,'(6X,6(1X,I5),F10.7)') IP,MV,LFTSYM,
     &                                      ICOP,ICP1,ICP2,CP
              END DO
            END DO
          END DO
        END DO
        WRITE(6,*)' 3. CLOSED LOOPS.'
        DO IP=1,NLEV
         DO IQ=1,IP
          INDEO=2*NLEV+(IP*(IP-1))/2+IQ
          DO MV=1,NMIDV
            DO LFTSYM=1,NSYM
              N=NOCP(INDEO,LFTSYM,MV)
              ICOP=IOCP(INDEO,LFTSYM,MV)
              DO I=1,N
                ICOP=ICOP+1
                ICP1=ICOUP(1,ICOP)
                ICP2=ICOUP(2,ICOP)
                CP=VTAB_TMP(ICOUP(3,ICOP))
                WRITE(6,'(7(1X,I5),F10.7)') IP,IQ,MV,LFTSYM,
     &                                      ICOP,ICP1,ICP2,CP
              END DO
            END DO
          END DO
         END DO
        END DO
      WRITE(6,*)
      WRITE(6,*)' CONVENTIONAL NR OF COUPLING COEFFS, BY PAIR:'
      NRC=0
      DO IP=2,MIDLEV
        DO IQ=1,IP-1
          NRCPQ=0
          INDEO=2*NLEV+(IP*(IP-1))/2+IQ
          DO LFTSYM=1,NSYM
            ISYM=LFTSYM
            DO MV=1,NMIDV
              NRCPQ=NRCPQ+NOCP(INDEO,LFTSYM,MV)*NOW(1,ISYM,MV)
            END DO
          END DO
          WRITE(6,'(1X,2I5,5X,I9)') IP,IQ,NRCPQ
          NRC=NRC+NRCPQ
        END DO
      END DO
      DO IP=MIDLEV+1,NLEV
        DO IQ=1,MIDLEV
          NRCPQ=0
          DO LFTSYM=1,NSYM
            ISYM=LFTSYM
            DO MV=1,NMIDV
              NRCPQ=NRCPQ+NOCP(IP,LFTSYM,MV)*NOCP(IQ,ISYM,MV)
              NRCPQ=NRCPQ+NOCP(NLEV+IP,LFTSYM,MV)*NOCP(NLEV+IQ,ISYM,MV)
            END DO
          END DO
          WRITE(6,'(1X,2I5,5X,I9)') IP,IQ,NRCPQ
          NRC=NRC+NRCPQ
        END DO
      END DO
      DO IP=MIDLEV+2,NLEV
        DO IQ=MIDLEV+1,IP-1
          NRCPQ=0
          INDEO=2*NLEV+(IP*(IP-1))/2+IQ
          DO LFTSYM=1,NSYM
            ISYM=LFTSYM
            DO MV=1,NMIDV
              NRCPQ=NRCPQ+NOCP(INDEO,LFTSYM,MV)*NOW(2,ISYM,MV)
            END DO
          END DO
          WRITE(6,'(1X,2I5,5X,I9)') IP,IQ,NRCPQ
          NRC=NRC+NRCPQ
        END DO
      END DO
      WRITE(6,*)
      WRITE(6,*)' TOTAL CONVENTIONAL COUPLING COEFFS:',NRC
#endif

      RETURN
      END
