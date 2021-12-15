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
      SUBROUTINE MKCOUP(nLev,iSm,nVert,MidLev,nMidV,MVSta,MVEnd,
     &                  MxEO,nICoup,nWalk,nICase,nVTab,
     &                  IVR,IMAW,ISGMNT,VSGMNT,NOW,IOW,
     &                  NOCP,IOCP,ILNDW,ICASE,ICOUP,VTAB,
     &                  ISGPTH,VALUE)
      IMPLICIT REAL*8 (A-H,O-Z)
#include "prgm.fh"
      CHARACTER*16 ROUTINE
      PARAMETER (ROUTINE='MKCOUP')
#include "symmul.fh"
      COMMON /SEGTAB/ IC1(26),IC2(26),ITVPT(26),IBVPT(26),ISVC(26),
     *                NIVR,LIVR,NSGMNT,LSGMNT
C Purpose: Compute and return the table ICOUP(1..3,ICOP).
C The number of coupling coeffs is obtained from NOCP, the offset to
C the ICOP numbering is given by IOCP. The numbers ICOUP(1..3,ICOP) are
C then the ket and bra half-walks, enumerated by the Lund scheme,
C and the index into the VTAB table (the pool of possible values of
C coupling coefficients).

C Any loop is regarded as a segment path from top to midlevel, or
C from midlevel to bottom.
C The segment path is described by the table ISGPTH. It is
C essentially a list of which one of segments nr 1..26 that are
C used at each level. The segments are of three types:
C Type 0: Upwalk segment or top loop segment.
C Type 1: Mid segment.
C Type 2: Bottom segment.
C Type 3: Downwalk segment.
C ISGPTH(IVLFT,LEV)=Left upper vertex.
C ISGPTH(ITYPE,LEV)=Type of segment, (0..3).
C ISGPTH(IAWSL,LEV)=Left arc weight sum (from top, or from bottom).
C ISGPTH(IAWSR,LEV)=Similar, right.
C ISGPTH(ILS  ,LEV)=Left symmetry label (accum from top or bottom).
C ISGPTH(ICS  ,LEV)=Left coupling case number (0..3).
C ISGPTH(ISEG ,LEV)=Segment type (1..26).
C These indices are used to denote the columns of table ISGPTH.
      PARAMETER (IVLFT=1,ITYPE=2,IAWSL=3,IAWSR=4,ILS=5,ICS=6)
      PARAMETER (ISEG=7)
C INPUT CALL PARAMETERS:
      DIMENSION ISM(NLEV)
      DIMENSION IVR(NVERT,2),IMAW(NVERT,0:3)
      DIMENSION ISGMNT(NVERT,26), VSGMNT(NVERT,26)
      DIMENSION NOW(2,NSYM,NMIDV), NOCP(MXEO,NSYM,NMIDV)
C OUTPUT CALL PARAMETERS:
      DIMENSION IOW(2,NSYM,NMIDV), IOCP(MXEO,NSYM,NMIDV)
      DIMENSION ILNDW(NWALK),ICASE(NICASE)
      DIMENSION VTAB(NVTAB)
      DIMENSION ICOUP(3,NICOUP)
C SCRATCH CALL PARAMETERS:
      DIMENSION ISGPTH(7,0:NLEV), VALUE(0:NLEV)




      nIpWlk=1+(MidLev-1)/15
      nIpWlk=max(nIpWlk,1+(nLev-MidLev-1)/15)
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
          DO  IS=1,NSYM
            NOCP(INDEO,IS,MV)=0
          END DO
        END DO
      END DO
C COUPLING COEFFICIENT VALUE TABLE:
      NVTAB=2
      VTAB(1)=1.0D00
      VTAB(2)=-1.0D00
      NCHECK=0
      DO IHALF=1,2
        IF(IHALF.EQ.1) THEN
          IVTSTA=1
          IVTEND=1
          LEV1=NLEV
          LEV2=MIDLEV
          ITYPMX=0
        ELSE
          IVTSTA=MVSTA
          IVTEND=MVEND
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
          ISGPTH(IVLFT,LEV)=IVTOP
          ISGPTH(ITYPE,LEV)=ITYP
          ISGPTH(IAWSL,LEV)=0
          ISGPTH(IAWSR,LEV)=0
          ISGPTH(ILS,LEV)=1
          ISGPTH(ISEG,LEV)=0
          VALUE(LEV)=1.0D00
100       CONTINUE
           IF(LEV.GT.LEV1) GOTO 400
           ITYPT=ISGPTH(ITYPE,LEV)
           IVLT=ISGPTH(IVLFT,LEV)
           DO ISGT=ISGPTH(ISEG,LEV)+1,26
             IVLB=ISGMNT(IVLT,ISGT)
             IF(IVLB.NE.0 .AND. ITYPT.EQ.ITVPT(ISGT)) GOTO 200
           END DO
           ISGPTH(ISEG,LEV)=0
           LEV=LEV+1
          GOTO 100

200       CONTINUE
          ISGPTH(ISEG,LEV)=ISGT
          ICL=IC1(ISGT)
          ISYM=1
          IF((ICL.EQ.1).OR.(ICL.EQ.2)) ISYM=ISM(LEV)
          IVRT=IVLT
          IF((ITYPT.EQ.1).OR.(ITYPT.EQ.2)) IVRT=IVR(IVLT,ITYPT)
          ICR=IC2(ISGT)
          ISGPTH(ICS,LEV)=ICL
          LEV=LEV-1
          ISGPTH(IAWSL,LEV)=ISGPTH(IAWSL,LEV+1)+IMAW(IVLT,ICL)
          ISGPTH(IAWSR,LEV)=ISGPTH(IAWSR,LEV+1)+IMAW(IVRT,ICR)
          VALUE(LEV)=VALUE(LEV+1)*VSGMNT(IVLT,ISGT)
          ISGPTH(ILS,LEV)=MUL(ISYM,ISGPTH(ILS,LEV+1))
          ISGPTH(IVLFT,LEV)=IVLB
          ISGPTH(ITYPE,LEV)=IBVPT(ISGT)
          ISGPTH(ISEG,LEV)=0
          IF (LEV.GT.LEV2) GOTO 100

          MV=ISGPTH(IVLFT,MIDLEV)+1-MVSTA
          LFTSYM=ISGPTH(ILS,LEV2)
          IT=ISGPTH(ITYPE,MIDLEV)
          IF(IT.EQ.0) IT=3
          IF(ISGPTH(ITYPE,LEV2).EQ.0) IT=0
          IF(IT.EQ.0) THEN
            ILND=1+NOW(IHALF,LFTSYM,MV)
            IAWS=ISGPTH(IAWSL,LEV2)
            ILNDW(IAWS)=ILND
            NOW(IHALF,LFTSYM,MV)=ILND
            IPOS=IOW(IHALF,LFTSYM,MV)+(ILND-1)*NIPWLK
            DO LL=LEV2+1,LEV1,15
              IC=0
              DO L=MIN(LL+14,LEV1),LL,-1
                IC=4*IC+ISGPTH(ICS,L)
              END DO
              IPOS=IPOS+1
              ICASE(IPOS)=IC
            END DO
          ELSE
            IP=0
            IQ=0
            DO L=LEV2+1,LEV1
              ISG=ISGPTH(ISEG,L)
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
     *                    ICOP-IOCP(INDEO,LFTSYM,MV)
              WRITE(6,*)' CURRENT NOCP NUMBER :',NOCP(INDEO,LFTSYM,MV)
              CALL ABEND()
            END IF
C
            C=VALUE(LEV2)

C Determine value code, IVTAB:
            DO I=1,NVTAB
              IVTAB=I
              IF(ABS(C-VTAB(I)).LT.1.0D-10) GOTO 212
            END DO
            NVTAB=NVTAB+1
            IF(NVTAB.GT.5000)THEN
              WRITE(6,*)' MKCOUP: NVTAB > 5000!'
              WRITE(6,*)' Need recompilation.'
              CALL ABEND()
            END IF
            VTAB(NVTAB)=C
            IVTAB=NVTAB
212         CONTINUE

            ICOUP(1,ICOP)=ISGPTH(IAWSL,LEV2)
            ICOUP(2,ICOP)=ISGPTH(IAWSR,LEV2)
            ICOUP(3,ICOP)=IVTAB
            IF (ICOP.GT.NICOUP) THEN
              WRITE(6,*)' ICOP > NICOUP!'
              WRITE(6,*)' (This should never happen!)'
              CALL ABEND()
            END IF
          END IF
          LEV=LEV+1
          GOTO 100
* End of long loops
400      CONTINUE
         END DO
        END DO
* End of a long loop
      END DO
C RENUMBER THE COUPLING COEFFICIENT INDICES BY LUND SCHEME:
      DO ICOP=1,NICOUP
        I1=ICOUP(1,ICOP)
        I2=ICOUP(2,ICOP)
        ICOUP(1,ICOP)=ILNDW(I1)
        ICOUP(2,ICOP)=ILNDW(I2)
      END DO

      IF(IPGLOB.GE.DEBUG) THEN
        ICOP1=0
        ICOP2=0
        WRITE(6,*)' NR OF DIFFERENT VALUES OF COUP:',NVTAB
        DO ICOP=1,NICOUP
          I3=ICOUP(3,ICOP)
          IF(I3.EQ.1) ICOP1=ICOP1+1
          IF(I3.EQ.2) ICOP2=ICOP2+1
        END DO
        WRITE(6,*)
        WRITE(6,*)' NR OF COUPS WITH VALUE  1.0:',ICOP1
        WRITE(6,*)' NR OF COUPS WITH VALUE -1.0:',ICOP2
      END IF
      IF(IPGLOB.GE.INSANE) THEN
        WRITE(6,*)
        WRITE(6,*)' COUPLING COEFFICIENTS:'
        WRITE(6,*)'    IP    IQ    MV LFTSYM NOCP'
        WRITE(6,*)' 1. OPEN LOOPS TYPE 1.'
        DO IP=1,NLEV
          DO MV=1,NMIDV
            DO LFTSYM=1,NSYM
              N=NOCP(IP,LFTSYM,MV)
              WRITE(6,'(6X,6(1X,I5),1X,F10.7)') IP,MV,LFTSYM,N
            END DO
          END DO
        END DO
        WRITE(6,*)' 2. OPEN LOOPS TYPE 2.'
        DO IP=1,NLEV
          INDEO=NLEV+IP
          DO MV=1,NMIDV
            DO LFTSYM=1,NSYM
              N=NOCP(INDEO,LFTSYM,MV)
              WRITE(6,'(6X,6(1X,I5),1X,F10.7)') IP,MV,LFTSYM,N
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
              WRITE(6,'(7(1X,I5),1X,F10.7)') IP,IQ,MV,LFTSYM,N
            END DO
          END DO
         END DO
        END DO
      END IF
      IF(IPGLOB.GE.INSANE) THEN
        WRITE(6,*)
        WRITE(6,*)' COUPLING COEFFICIENTS:'
        WRITE(6,*)'    IP    IQ    MV LFTSYM ICOP   ICOUP1&2   COUP'
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
                CP=VTAB(ICOUP(3,ICOP))
                WRITE(6,'(6X,6(1X,I5),1X,F10.7)')
     &          IP,MV,LFTSYM,ICOP,ICP1,ICP2,CP
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
                CP=VTAB(ICOUP(3,ICOP))
                WRITE(6,'(6X,6(1X,I5),1X,F10.7)')
     &            IP,MV,LFTSYM,ICOP,ICP1,ICP2,CP
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
                CP=VTAB(ICOUP(3,ICOP))
                WRITE(6,'(7(1X,I5),1X,F10.7)')
     &            IP,IQ,MV,LFTSYM,ICOP,ICP1,ICP2,CP
              END DO
            END DO
          END DO
         END DO
        END DO
      END IF
      IF(IPGLOB.GE.INSANE) THEN
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
      END IF

      RETURN
      END
