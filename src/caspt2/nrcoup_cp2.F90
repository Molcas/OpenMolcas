!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 1994, Per Ake Malmqvist                                *
!***********************************************************************
!--------------------------------------------*
! 1994  PER-AAKE MALMQUIST                   *
! DEPARTMENT OF THEORETICAL CHEMISTRY        *
! UNIVERSITY OF LUND                         *
! SWEDEN                                     *
!--------------------------------------------*
      SUBROUTINE NRCOUP_CP2(  &
     &                  NVERT,NMIDV,MXEO,ISM,  &
     &                  IDRT,ISGMNT,NOW,IOW,NOCP,IOCP,                  &
     &                  NOCSF,IOCSF,NCSF,   &
     &                  NRL,MVL,MVR,NICOUP)
      use UNIXInfo, only: ProgName
      use gugx, only: SGS, EXS, CIS
      use stdalloc, only: mma_allocate
      IMPLICIT None

#include "rasdim.fh"
#include "caspt2.fh"
#include "pt2_guga.fh"
#include "segtab.fh"

! INPUT PARAMETERS:
      Integer, Intent(In) :: nVert, nMidV, MxEO
      Integer MVL(NMIDV,2),MVR(NMIDV,2)
      Integer IDRT(NVERT,5),ISGMNT(NVERT,26)
      Integer ISM(*)
      Integer, PARAMETER :: LTAB=1
! OUTPUT PARAMETERS:
      Integer, Intent(Out):: nICoup
      Integer NOW(2,NSYM,NMIDV),NOCP(MXEO,NSYM,NMIDV)
      Integer IOW(2,NSYM,NMIDV),IOCP(MXEO,NSYM,NMIDV)
      Integer NOCSF(NSYM,NMIDV,NSYM),IOCSF(NSYM,NMIDV,NSYM)
      Integer NCSF(NSYM)
! SCRATCH PARAMETERS:
      Integer NRL(NSYM,NVERT,0:MXEO)

      Integer IBSYM, ICL, INDEO, INDEOB, INDEOT, IP, IPQ, IQ, ISGT,     &
     &        ISYDS1, ISYM, ISYUS1, ITSYM, IV, IVLB, IVLT, LEV, LFTSYM, &
     &        MV, MV1, MV2, MV3, MV4, MV5, MXDWN, MXUP, N, NSGMX,       &
     &        NT1TMP, NT2TMP, NT3TMP, NT4TMP, NUPS1, NSGTMP, NLW, NUW,  &
     &        NWALK, ISYDWN, ISYTOT, ISYUP, NDWNS1
#ifdef _DEBUGPRINT_
      Integer IS, IST, NCP, NLW, NUW
#endif
      Logical :: IF_RASSI=.FALSE.
      Integer :: nLev, MVSTA, MVEND, NIPWLK

! Dereference SGS, CIS for some other data
      nLev  =SGS%nLev
      MVSTA =SGS%MVSta
      MVEND =SGS%MVEnd
      NIPWLK=CIS%nIpWlk

      Select Case (ProgName(1:6))
         Case ('rassi')
            IF_RASSI = .TRUE.
         Case Default
            IF_RASSI = .FALSE.
      End Select

      DO INDEO=0,MXEO
        DO IV=1,MVEnd
          DO LFTSYM=1,NSYM
            NRL(LFTSYM,IV,INDEO)=0
          END DO
        END DO
      END DO
      NRL(1,1,0)=1

      DO IVLT=1,MVSta-1
        LEV=IDRT(IVLT,LTAB)
        DO ISGT=1,26
          IVLB=ISGMNT(IVLT,ISGT)
          IF (IVLB.EQ.0) Cycle
          ICL=IC1(ISGT)
          ISYM=1
          IF ((ICL.EQ.1).OR.(ICL.EQ.2)) ISYM=ISM(LEV)
          DO ITSYM=1,NSYM
            IBSYM=MUL(ITSYM,ISYM)

            IF (ISGT.LE.4) THEN
! THIS IS AN UPPER WALK.
              NRL(IBSYM,IVLB,0)=NRL(IBSYM,IVLB,0)+NRL(ITSYM,IVLT,0)
              Cycle
            END IF
            IF(ISGT.LE.8) THEN
! THIS IS AN TOP SEGMENT.
              INDEO=LEV+(IBVPT(ISGT)-1)*NLEV
              NRL(IBSYM,IVLB,INDEO)=NRL(IBSYM,IVLB,INDEO)+NRL(ITSYM,IVLT,0)
              Cycle
            END IF
            IF (ISGT.LE.18) THEN
! THIS IS A MID-SEGMENT.
              DO IP=LEV+1,NLEV
                INDEOT = IP+(ITVPT(ISGT)-1)*NLEV
                INDEOB = IP+(IBVPT(ISGT)-1)*NLEV
                NRL(IBSYM,IVLB,INDEOB) = NRL(IBSYM,IVLB,INDEOB) + NRL(ITSYM,IVLT,INDEOT)
              END DO
              Cycle
            END IF
            IF (ISGT.LE.22) THEN
! THIS IS A BOTTOM SEGMENT.
              DO IP=LEV+1,NLEV
                INDEOT=IP+(ITVPT(ISGT)-1)*NLEV
                IPQ=(IP*(IP-1))/2 + LEV
                INDEOB=IPQ+2*NLEV
                NRL(IBSYM,IVLB,INDEOB) = NRL(IBSYM,IVLB,INDEOB) + NRL(ITSYM,IVLT,INDEOT)
              END DO
              Cycle
            END IF
! THIS IS A LOWER WALK.
            DO INDEO=2*NLEV+1,MXEO
              NRL(IBSYM,IVLB,INDEO) = NRL(IBSYM,IVLB,INDEO) + NRL(ITSYM,IVLT,INDEO)
            END DO
          END DO
        END DO
      END DO

      MXUP=0
      DO MV=1,NMIDV
        IVLT=MV+MVSta-1
        DO LFTSYM=1,NSYM
          IF (IF_RASSI) NOW(1,LFTSYM,MV)=NRL(LFTSYM,IVLT,0)
          MXUP=MAX(MXUP,NOW(1,LFTSYM,MV))
          DO INDEO=1,MXEO
            NOCP(INDEO,LFTSYM,MV)=NRL(LFTSYM,IVLT,INDEO)
          END DO
        END DO
      END DO

      DO INDEO=0,MXEO
        DO IV=MVSta,NVERT
          DO LFTSYM=1,NSYM
            NRL(LFTSYM,IV,INDEO)=0
          END DO
        END DO
      END DO

      NRL(1,NVERT,0)=1
      DO 201 IVLT=NVERT-1,MVSta,-1
        LEV=IDRT(IVLT,LTAB)
        DO 202 ISGT=1,26
          IVLB=ISGMNT(IVLT,ISGT)
          IF(IVLB.EQ.0) GOTO 202
          ICL=IC1(ISGT)
          ISYM=1
          IF((ICL.EQ.1).OR.(ICL.EQ.2)) ISYM=ISM(LEV)
          DO 200 ITSYM=1,NSYM
            IBSYM=MUL(ITSYM,ISYM)
      IF(ISGT.GE.23) THEN
! THIS IS A LOWER WALK.
        NRL(ITSYM,IVLT,0)=NRL(ITSYM,IVLT,0)+NRL(IBSYM,IVLB,0)
        GOTO 200
      END IF
      IF(ISGT.GE.19) THEN
! THIS IS AN BOTTOM SEGMENT.
        INDEO=LEV+(ITVPT(ISGT)-1)*NLEV
        NRL(ITSYM,IVLT,INDEO)=NRL(ITSYM,IVLT,INDEO)+NRL(IBSYM,IVLB,0)
        GOTO 200
      END IF
      IF(ISGT.GE.9) THEN
! THIS IS A MID-SEGMENT.
        DO 210 IQ=1,LEV-1
          INDEOT=IQ+(ITVPT(ISGT)-1)*NLEV
          INDEOB=IQ+(IBVPT(ISGT)-1)*NLEV
          NRL(ITSYM,IVLT,INDEOT) = NRL(ITSYM,IVLT,INDEOT) + NRL(IBSYM,IVLB,INDEOB)
 210    CONTINUE
        GOTO 200
      END IF
      IF(ISGT.GE.5) THEN
! THIS IS AN TOP SEGMENT.
        DO 220 IQ=1,LEV-1
          INDEOB=IQ+(IBVPT(ISGT)-1)*NLEV
          IPQ=(LEV*(LEV-1))/2+IQ
          INDEOT=IPQ+2*NLEV
          NRL(ITSYM,IVLT,INDEOT) = NRL(ITSYM,IVLT,INDEOT) + NRL(IBSYM,IVLB,INDEOB)
 220    CONTINUE
        GOTO 200
      END IF
! THIS IS AN UPPER WALK.
      DO 230 IPQ=1,(LEV*(LEV-1))/2
        INDEO=2*NLEV+IPQ
        NRL(ITSYM,IVLT,INDEO) = NRL(ITSYM,IVLT,INDEO) + NRL(IBSYM,IVLB,INDEO)
 230  CONTINUE

 200  CONTINUE
 202  CONTINUE
 201  CONTINUE

      MXDWN=0
      DO MV=1,NMIDV
        IVLT=MV+MVSta-1
        DO LFTSYM=1,NSYM
          IF (IF_RASSI) NOW(2,LFTSYM,MV)=NRL(LFTSYM,IVLT,0)
          MXDWN=MAX(MXDWN,NOW(2,LFTSYM,MV))
          DO INDEO=1,MXEO
            N=NRL(LFTSYM,IVLT,INDEO)
            IF(N.NE.0) NOCP(INDEO,LFTSYM,MV)=N
          END DO
        END DO
      END DO

      IF (IF_RASSI) Then
      NUW=0
      DO MV=1,NMIDV
        DO ISYM=1,NSYM
          IOW(1,ISYM,MV)=NUW*NIPWLK
          NUW=NUW+NOW(1,ISYM,MV)
        END DO
      END DO
      NWALK=NUW
      DO MV=1,NMIDV
        DO ISYM=1,NSYM
          IOW(2,ISYM,MV)=NWALK*NIPWLK
          NWALK=NWALK+NOW(2,ISYM,MV)
        END DO
      END DO
      NLW=NWALK-NUW
      End If

      NICOUP=0
      DO INDEO=1,MXEO
        DO MV=1,NMIDV
          DO LFTSYM=1,NSYM
            IOCP(INDEO,LFTSYM,MV)=NICOUP
            NICOUP=NICOUP+NOCP(INDEO,LFTSYM,MV)
          END DO
        END DO
      END DO

      If (IF_RASSI) Then
      DO ISYTOT=1,NSYM
        NCSF(ISYTOT)=0
        DO MV=1,NMIDV
          DO ISYUP=1,NSYM
            ISYDWN=MUL(ISYTOT,ISYUP)
            N=NOW(1,ISYUP,MV)*NOW(2,ISYDWN,MV)
            NOCSF(ISYUP,MV,ISYTOT)=N
            IOCSF(ISYUP,MV,ISYTOT)=NCSF(ISYTOT)
            NCSF(ISYTOT)=NCSF(ISYTOT)+N
          END DO
        END DO
      END DO
      End If

!AR   INSERT FOR US IN SIGMA ROUTINE
!
      NSGMX=1
      NSGTMP=MAX(MXUP,MXDWN)
      DO 550 MV3=1,NMIDV
        MV1=MVL(MV3,2)
        MV2=MVL(MV3,1)
        MV4=MVR(MV3,1)
        MV5=MVR(MV3,2)
      DO 551 ISYUS1=1,NSYM
        NUPS1=NOW(1,ISYUS1,MV3)
      DO 552 ISYDS1=1,NSYM
        If (IF_RASSI) NDWNS1=NOW(2,ISYDS1,MV3)
        NSGMX=MAX(NSGMX,NOCSF(ISYUS1,MV3,ISYDS1))
        IF (MV1.NE.0)THEN
          NT4TMP=NUPS1*NOW(2,ISYDS1,MV1)
          NSGTMP=MAX(NSGTMP,NT4TMP)
        ENDIF
        IF (MV2.NE.0)THEN
          NT3TMP=NUPS1*NOW(2,ISYDS1,MV2)
          NSGTMP=MAX(NSGTMP,NT3TMP)
        ENDIF
        IF (MV4.NE.0)THEN
          NT1TMP=NUPS1*NOW(2,ISYDS1,MV4)
          NSGTMP=MAX(NSGTMP,NT1TMP)
        ENDIF
        IF (MV5.NE.0)THEN
          NT2TMP=NUPS1*NOW(2,ISYDS1,MV5)
          NSGTMP=MAX(NSGTMP,NT2TMP)
        ENDIF
  552 CONTINUE
  551 CONTINUE
  550 CONTINUE
      CALL mma_allocate(EXS%SGTMP,NSGTMP,Label='EXS%SGTMP')
!
#ifdef _DEBUGPRINT_
      NUW=NOW(1,NSYM,NMIDV)
      NLW=NOW(2,NSYM,NMIDV)-NOW(1,NSYM,NMIDV)
      WRITE(6,600)MXUP,MXDWN,                                           &
     &            NSGMX,NSGMX,NSGTMP
  600 FORMAT(/,' MAXIMUM NUMBER OF WALKS',                              &
     &       /,' UPPER ',I6,' LOWER ',I6,                               &
     &       /,' LENGTH OF LARGEST WORK ARRAYS IN SUBROUTINE SIGMA',    &
     &       /,' TEMPORARY SGM1 ',I7,                                   &
     &       /,' TEMPORARY SGM2 ',I7,                                   &
     &       /,' NSGTMP         ',I7)
!
!AR   END OF INSERT
        WRITE(6,*)
        WRITE(6,*)' TOTAL NR OF WALKS: UPPER ',NUW
        WRITE(6,*)'                    LOWER ',NLW
        WRITE(6,*)'                     SUM  ',CIS%NWALK
        WRITE(6,*)' TOTAL NR OF COUPL COEFFS ',NICOUP
        INDEO=2*NLEV+1
        WRITE(6,*)'         OF TYPE 1&2 ONLY:',IOCP(INDEO,1,1)
        WRITE(6,*)
        WRITE(6,*)' NR OF CONFIGURATIONS/SYMM:'
        WRITE(6,'(8(1X,I8))')(CIS%NCSF(IS),IS=1,NSYM)
        WRITE(6,*)

        WRITE(6,*)
        WRITE(6,*)' NR OF WALKS AND CONFIGURATIONS IN NRCOUP'
        WRITE(6,*)' BY MIDVERTEX AND SYMMETRY.'
        DO 310 MV=1,NMIDV
          WRITE(6,*)
          WRITE(6,1234) MV,(NOW(1,IS,MV),IS=1,NSYM)
          WRITE(6,1235)    (NOW(2,IS,MV),IS=1,NSYM)
          DO 305 IST=1,NSYM
            WRITE(6,1236)IST,(NOCSF(IS,MV,IST),IS=1,NSYM)
 305      CONTINUE
 1234 FORMAT('  MV=',I2,'    UPPER WALKS:',8I6)
 1235 FORMAT('           LOWER WALKS:',8I6)
 1236 FORMAT(' IST=',I2,'  CONFIGURATIONS:',8I6)
 310    CONTINUE
        WRITE(6,*)
        WRITE(6,*)' NR OF COUPLING COEFFICIENTS:'
        WRITE(6,*)' 1. OPEN LOOPS TYPE 1:'
        DO 320 IP=1,NLEV
          DO 321 MV=1,NMIDV
            DO 322 IS=1,NSYM
              NCP=NOCP(IP,IS,MV)
              IF(NCP.EQ.0) GOTO 322
              WRITE(6,2345) IP,MV,IS,NCP
 2345 FORMAT(' P=',I2,'  MV=',I2,' SYMM ',I1,' NOCP=',I4)
 322        CONTINUE
 321      CONTINUE
 320    CONTINUE
        WRITE(6,*)
        WRITE(6,*)' 2. OPEN LOOPS TYPE 2:'
        DO 330 IP=1,NLEV
          DO 331 MV=1,NMIDV
            DO 332 IS=1,NSYM
              NCP=NOCP(NLEV+IP,IS,MV)
              IF(NCP.EQ.0) GOTO 332
              WRITE(6,2345) IP,MV,IS,NCP
 332       CONTINUE
 331     CONTINUE
 330    CONTINUE
        WRITE(6,*)
        WRITE(6,*)' 3. CLOSED LOOPS:'
        DO 340 IP=2,NLEV
          DO 341 IQ=1,IP-1
            INDEO=2*NLEV+(IP*(IP-1))/2+IQ
            DO 342 MV=1,NMIDV
              DO 343 IS=1,NSYM
                NCP=NOCP(INDEO,IS,MV)
                IF(NCP.EQ.0) GOTO 343
                WRITE(6,2346) IP,IQ,MV,IS,NCP
 2346 FORMAT(' P=',I2,'  Q=',I2,'  MV=',I2,' SYMM ',I1,' NOCP=',I4)
 343          CONTINUE
 342        CONTINUE
 341      CONTINUE
 340    CONTINUE
#endif

      END SUBROUTINE NRCOUP_CP2
