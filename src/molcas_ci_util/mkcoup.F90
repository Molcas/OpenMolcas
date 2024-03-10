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
      SUBROUTINE MKCOUP(SGS,CIS,EXS)

      use Symmetry_Info, only: Mul
      use stdalloc, only: mma_allocate, mma_deallocate
      use struct, only: SGStruct, CIStruct, EXStruct
      IMPLICIT None

#include "segtab.fh"
! Purpose: Compute and return the table ICOUP(1..3,ICOP).
! The number of coupling coeffs is obtained from NOCP, the offset to
! the ICOP numbering is given by IOCP. The numbers ICOUP(1..3,ICOP) are
! then the ket and bra half-walks, enumerated by the Lund scheme,
! and the index into the VTAB table (the pool of possible values of
! coupling coefficients).

! Any loop is regarded as a segment path from top to midlevel, or
! from midlevel to bottom.
! The segment path is described by the table ISGPTH. It is
! essentially a list of which one of segments nr 1..26 that are
! used at each level. The segments are of three types:
! Type 0: Upwalk segment or top loop segment.
! Type 1: Mid segment.
! Type 2: Bottom segment.
! Type 3: Downwalk segment.
! ISGPTH(IVLFT,LEV)=Left upper vertex.
! ISGPTH(ITYPE,LEV)=Type of segment, (0..3).
! ISGPTH(IAWSL,LEV)=Left arc weight sum (from top, or from bottom).
! ISGPTH(IAWSR,LEV)=Similar, right.
! ISGPTH(ILS  ,LEV)=Left symmetry label (accum from top or bottom).
! ISGPTH(ICS  ,LEV)=Left coupling case number (0..3).
! ISGPTH(ISEG ,LEV)=Segment type (1..26).
! These indices are used to denote the columns of table ISGPTH.

! INPUT PARAMETERS:
      Type(SGStruct) SGS
      Type(CIStruct) CIS
      Type(EXStruct) EXS
! OUTPUT PARAMETERS:
      Integer, Parameter:: nVTab=5000
! SCRATCH PARAMETERS:
      Integer, PARAMETER :: IVLFT=1,ITYPE=2,IAWSL=3,IAWSR=4,ILS=5,ICS=6
      Integer, PARAMETER :: ISEG=7
      Integer, Allocatable:: ILNDW(:), ISGPTH(:,:)
      Real*8, Allocatable:: Value(:)
      Real*8, Allocatable:: VTab(:)

! local stuff
      Integer :: nVTab_Final, i, i1, i2, IAWS, IC, ICL, ICOP, ICR,      &
     &           IHALF, iLnd, IndEO, iP, iPos, iQ, iS, iSg, iSgt,       &
     &           iSym, iT, iTyp, iTypMx, iTypT, iVlb, iVlt, iVrt,       &
     &           iVrTop, iVTab, iVTEnd, iVTop, iVTSta, L, Lev, Lev1,    &
     &           Lev2, LftSym, LL, MV, nCheck
      Real*8 :: C
#ifdef _DEBUGPRINT_
      Real*8 :: CP
      Integer :: I3, ICOP1, ICOP2, ICP1, ICP2, N, NRC, NRCPQ
#endif

      Call mma_allocate(ILNDW,CIS%nWalk,Label='ILNDW')
      Call mma_allocate(ISGPTH,[1,7],[0,SGS%nLev],Label='ISGPTH')
      Call mma_allocate(Value,[0,SGS%nLev],Label='Value')
      Call mma_allocate(VTab,nVTab,Label='VTab')

      Call mma_allocate(EXS%ICoup,3,EXS%nICoup,Label='EXS%ICoup')
      If (.NOT.Allocated(CIS%ICase)) Call mma_allocate(CIS%ICase,       &
     &                    CIS%nWalk*CIS%nIpWlk,Label='CIS%ICase')

      Associate (ICoup=>EXS%ICoup, nICoup=>EXS%nICoup,                  &
     &           ISm=>SGS%ISm, IVR=>CIS%IVR, iMAW=>SGS%MAW,             &
     &           nLev=>SGS%nLev, ISGMNT=>CIS%ISGM,                      &
     &           VSGMNT=>CIS%VSGM, NOW=>CIS%NOW, IOW=>CIS%IOW,          &
     &           nVert=>SGS%nVert, iCase=>CIS%ICase,                    &
     &           NOCP=>EXS%NOCP, IOCP=>EXS%IOCP, nSym=>SGS%nSym,        &
     &           MxEO=>EXS%MxEO, nMidV=>CIS%nMidV, MidLev=>SGS%MidLev,  &
     &           MVSta=>SGS%MVSta, MVEnd=>SGS%MVEnd,                    &
     &           nWalk=>CIS%nWalk, nIpWlk=>CIS%nIpWlk)

!     nIpWlk=1+(MidLev-1)/15
!     nIpWlk=max(nIpWlk,1+(nLev-MidLev-1)/15)
! NOW IS ZEROED AND WILL BE USED AS AN ARRAY OF COUNTERS, BUT WILL
!    BE RESTORED FINALLY.
      DO IHALF=1,2
        DO MV=1,NMIDV
          DO IS=1,NSYM
            NOW(IHALF,IS,MV)=0
          END DO
        END DO
      END DO
! SIMILAR FOR THE COUPLING COEFFICIENT TABLE:
      DO INDEO=1,MXEO
        DO MV=1,NMIDV
          DO IS=1,NSYM
            NOCP(INDEO,IS,MV)=0
          END DO
        END DO
      END DO

! COUPLING COEFFICIENT VALUE TABLE:
      NVTAB_FINAL=2
      VTab(1)=1.0D00
      VTab(2)=-1.0D00

      NCHECK=0

      DO IHALF=1,2
        IF(IHALF.EQ.1) THEN
          IVTSTA=1
          IVTEND=1
          LEV1=NLEV
          LEV2=MIDLEV
          ITYPMX=0
        ELSE
          IVTSTA=MVSta
          IVTEND=MVEnd
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
 100      IF(LEV.GT.LEV1) GOTO 400
          ITYPT=ISGPTH(ITYPE,LEV)
          IVLT=ISGPTH(IVLFT,LEV)
          DO ISGT=ISGPTH(ISEG,LEV)+1,26
            IVLB=ISGMNT(IVLT,ISGT)
            IF(IVLB.EQ.0) GOTO 110
            IF(ITYPT.EQ.ITVPT(ISGT)) GOTO 200
 110        CONTINUE
          END DO
          ISGPTH(ISEG,LEV)=0
          LEV=LEV+1
          GOTO 100

 200      ISGPTH(ISEG,LEV)=ISGT
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

          MV=ISGPTH(IVLFT,MIDLEV)+1-MVSta
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
              WRITE(6,*)' D:O, WITHOUT OFFSET :',                       &
     &                    ICOP-IOCP(INDEO,LFTSYM,MV)
              WRITE(6,*)' CURRENT NOCP NUMBER :',NOCP(INDEO,LFTSYM,MV)
              CALL ABEND()
            END IF
!
            C=VALUE(LEV2)
            DO I=1,NVTAB_FINAL
              IVTAB=I
              IF(ABS(C-VTab(I)).LT.1.0D-10) GOTO 212
            END DO
            NVTAB_FINAL=NVTAB_FINAL+1
            IF(NVTAB_FINAL.GT.nVTab) THEN
              WRITE(6,*)'MKCOUP: NVTAB_FINAL=',NVTAB_FINAL
              WRITE(6,*)'NVTAB_FINAL should not be allowed to grow'
              WRITE(6,*)'beyond nVTab which was set provisionally'
              WRITE(6,*)'in subroutine GINIT in file ginit.f.'
              WRITE(6,*)'Now nVTab=',nVTab
              WRITE(6,*)'This may indicate a problem with your input.'
              WRITE(6,*)'If you do want to do this big calculation, try'
              WRITE(6,*)'increasing nVTab in GINIT and recompile.'
              CALL ABEND()
            END IF
            VTab(NVTAB_FINAL)=C
            IVTAB=NVTAB_FINAL
 212        ICOUP(1,ICOP)=ISGPTH(IAWSL,LEV2)
            ICOUP(2,ICOP)=ISGPTH(IAWSR,LEV2)
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
! RENUMBER THE COUPLING COEFFICIENT INDICES BY LUND SCHEME:
      DO ICOP=1,NICOUP
        I1=ICOUP(1,ICOP)
        I2=ICOUP(2,ICOP)
        ICOUP(1,ICOP)=ILNDW(I1)
        ICOUP(2,ICOP)=ILNDW(I2)
      END DO

      Call mma_deallocate(Value)
      Call mma_deallocate(ISGPTH)
      Call mma_deallocate(ILNDW)

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
                CP=VTab(ICOUP(3,ICOP))
                WRITE(6,'(6X,6(1X,I5),F10.7)') IP,MV,LFTSYM,            &
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
                CP=VTab(ICOUP(3,ICOP))
                WRITE(6,'(6X,6(1X,I5),F10.7)') IP,MV,LFTSYM,            &
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
                CP=VTab(ICOUP(3,ICOP))
                WRITE(6,'(7(1X,I5),F10.7)') IP,IQ,MV,LFTSYM,            &
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

      End Associate

      Call mma_allocate(EXS%VTab,nVTab_Final,Label='EXS%VTab')
      EXS%VTab(1:nVTab_final) = VTab(1:nVTab_final)
      Call mma_deallocate(VTab)

      END SUBROUTINE MKCOUP
