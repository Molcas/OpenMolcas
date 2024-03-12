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

      SUBROUTINE SIGMA1(SGS,CIS,EXS,IP,IQ,CPQ,ISYCI,CI,SGM)
      use struct, only: SGStruct, CIStruct, EXStruct
      use Symmetry_Info, only: Mul
      IMPLICIT REAL*8 (A-H,O-Z)
      Type (SGStruct) SGS
      Type (CIStruct) CIS
      Type (EXStruct) EXS
      Integer, Intent(in):: IP, IQ, ISYCI
      Real*8, Intent(in):: CPQ
      Real*8 CI(*)
      Real*8 SGM(*)

      INTRINSIC MOD

      Associate (ISM => SGS%ISM, nLev => SGS%nLev, MidLev => SGS%MidLev,&
     &           nIpWlk => CIS%nIpWlk, MVL => EXS%MVL,                  &
     &           MVR => EXS%MVR, ICOUP => EXS%ICOUP,                    &
     &           NOCP => EXS%NOCP, IOCP => EXS%IOCP,                    &
     &           VTab => EXS%VTab, ICASE => CIS%ICASE,                  &
     &           NOW => CIS%NOW, IOW => CIS%IOW,                        &
     &           NOCSF => CIS%NOCSF, IOCSF => CIS%IOCSF,                &
     &           nSym => SGS%nSym, nMidV => CIS%nMidV,                  &
     &           SGTMP=>EXS%SGTMP)

!****************************************************************
!  GIVEN ACTIVE LEVEL INDICES IP AND IQ, AND INPUT CI ARRAYS
!  CI AND SGM, THIS ROUTINE ADDS TO SGM THE RESULT OF ACTING ON
!  CI WITH THE NUMBER CPQ TIMES THE EXCITATION OPERATOR E(IP,IQ).
!  THE ADDITIONAL ENTRIES IN THE PARAMETER LIST ARE TABLES THAT
!  WERE PREPARED BY GINIT AND ITS SUBROUTINES.
!****************************************************************

      IF(ABS(CPQ).LT.1.0D-12) RETURN

! SYMMETRY OF ORBITALS:
      ISYP=ISM(IP)
      ISYQ=ISM(IQ)
      ISYPQ=MUL(ISYP,ISYQ)
! SYMMETRY OF SIGMA ARRAY:
      ISYSGM=MUL(ISYPQ,ISYCI)

IF(IP.GT.IQ) THEN
! THEN THIS IS AN EXCITING OPERATOR.
  IF(IP.LE.MIDLEV) Then

! EXCITING CASE, IQ<IP<=MIDLEV.
      DO MVSGM=1,NMIDV
        DO ISYUSG=1,NSYM
          NS1=NOCSF(ISYUSG,MVSGM,ISYSGM)
          IF(NS1.EQ.0) Cycle
          ISYDSG=MUL(ISYUSG,ISYSGM)
          ISYDC=MUL(ISYPQ,ISYDSG)
          NDWNC=NOW(2,ISYDC,MVSGM)
          IF(NDWNC.EQ.0) Cycle
          ISGSTA=1+IOCSF(ISYUSG,MVSGM,ISYSGM)
          NUPSG=NOW(1,ISYUSG,MVSGM)
          IOC=IOCSF(ISYUSG,MVSGM,ISYCI)
          INDEO=2*NLEV+(IP*(IP-1))/2+IQ
          NCP=NOCP(INDEO,ISYDC,MVSGM)
          IF(NCP.GT.0) THEN
            LICP=1+IOCP(INDEO,ISYDC,MVSGM)
! CASE IS: LOWER HALF, EXCITE:
            CALL EXC1 (CPQ,NUPSG,CI(IOC+1),SGM(ISGSTA),NCP,ICOUP(1,LICP),VTAB)
          END IF
        End Do
      End Do

  ELSE IF(IQ.GT.MIDLEV) Then

! EXCITING CASE, MIDLEV<IQ<IP
      DO MVSGM=1,NMIDV
        DO ISYUSG=1,NSYM
          NS1=NOCSF(ISYUSG,MVSGM,ISYSGM)
          IF(NS1.EQ.0) Cycle
          ISYUC=MUL(ISYPQ,ISYUSG)
          NUPC=NOW(1,ISYUC,MVSGM)
          IF (NUPC.EQ.0) Cycle
          ISGSTA=1+IOCSF(ISYUSG,MVSGM,ISYSGM)
          NUPSG=NOW(1,ISYUSG,MVSGM)
          ISYDSG=MUL(ISYUSG,ISYSGM)
          NDWNSG=NOW(2,ISYDSG,MVSGM)
          IOC=IOCSF(ISYUC,MVSGM,ISYCI)
          INDEO=2*NLEV+(IP*(IP-1))/2+IQ
          NCP=NOCP(INDEO,ISYUC,MVSGM)
          IF(NCP.GT.0) THEN
            LICP=1+IOCP(INDEO,ISYUC,MVSGM)
! CASE IS: UPPER HALF, EXCITE:
            CALL EXC2 (CPQ,NDWNSG,NUPC,CI(IOC+1),NUPSG,SGM(ISGSTA),NCP,ICOUP(1,LICP),VTAB)
          END IF
        End Do
      End Do

  ELSE

! EXCITING CASE, IQ<=MIDLEV<IP
      DO MVSGM=1,NMIDV
        MV1=MVL(MVSGM,2)
        MV2=MVL(MVSGM,1)
        IF((MV1.EQ.0).AND.(MV2.EQ.0)) Cycle
        DO ISYUSG=1,NSYM
          NS1=NOCSF(ISYUSG,MVSGM,ISYSGM)
          IF(NS1.EQ.0) Cycle
          ISGSTA=1+IOCSF(ISYUSG,MVSGM,ISYSGM)
          NUPSG=NOW(1,ISYUSG,MVSGM)
          ISYDSG=MUL(ISYUSG,ISYSGM)
          ISYUC=MUL(ISYP,ISYUSG)
          ISYDC=MUL(ISYQ,ISYDSG)
          IF (MV2/=0) Then
            NUPC=NOW(1,ISYUC,MV2)
            IF(NUPC/=0) Then
              NDWNC=NOW(2,ISYDC,MV2)
              IF(NDWNC/=0) Then
                INDEO=IP
                NCP=NOCP(INDEO,ISYUC,MV2)
                IF(NCP/=0) Then
                  NTMP=NUPSG*NDWNC
                  SGTMP(1:NTMP)=0.0D0
                  LICP=1+IOCP(INDEO,ISYUC,MV2)
                  IOC=IOCSF(ISYUC,MV2,ISYCI)
! CASE IS: UPPER HALF, EXCITE:
                  CALL EXC2 (CPQ,NDWNC,NUPC,CI(IOC+1),NUPSG,SGTMP,NCP,ICOUP(1,LICP),VTAB)
                  INDEO=IQ
                  NCP=NOCP(INDEO,ISYDC,MV2)
                  IF(NCP/=0) Then
                    LICP=1+IOCP(INDEO,ISYDC,MV2)
! CASE IS: LOWER HALF, EXCITE:
                    X=1.0D00
                    CALL EXC1 (X,NUPSG,SGTMP,SGM(ISGSTA),NCP,ICOUP(1,LICP),VTAB)
                End If
                End If
              End If
            End If
          End If

          IF(MV1/=0) Then
            NUPC=NOW(1,ISYUC,MV1)
            IF(NUPC/=0) Then
              NDWNC=NOW(2,ISYDC,MV1)
              IF(NDWNC/=0) Then
                INDEO=NLEV+IP
                NCP=NOCP(INDEO,ISYUC,MV1)
                IF(NCP/=0) Then
                  NTMP=NUPSG*NDWNC
                  SGTMP(1:NTMP)=0.0D0
                  LICP=1+IOCP(INDEO,ISYUC,MV1)
                  IOC=IOCSF(ISYUC,MV1,ISYCI)
! CASE IS: UPPER HALF, EXCITE:
                  CALL EXC2 (CPQ,NDWNC,NUPC,CI(IOC+1),NUPSG,SGTMP,NCP,ICOUP(1,LICP),VTAB)
                  INDEO=NLEV+IQ
                  NCP=NOCP(INDEO,ISYDC,MV1)
                  IF(NCP/=0) Then
                    LICP=1+IOCP(INDEO,ISYDC,MV1)
! CASE IS: LOWER HALF, EXCITE:
                    X=1.0D00
                    CALL EXC1 (X,NUPSG,SGTMP,SGM(ISGSTA),NCP,ICOUP(1,LICP),VTAB)
                  End If
                End If
              End If
            End If
          End If
        END DO
      END DO

  END IF
ELSE IF(IQ.GT.IP) THEN
! THEN THIS IS A DEEXCITING OPERATOR.
  IF(IQ.LE.MIDLEV) Then

! DEEXCITING OPERATOR, IP<IQ<=MIDLEV.
      DO MVSGM=1,NMIDV
        DO ISYUSG=1,NSYM
          NS1=NOCSF(ISYUSG,MVSGM,ISYSGM)
          IF(NS1.EQ.0) Cycle
          ISYDSG=MUL(ISYUSG,ISYSGM)
          ISYDC=MUL(ISYPQ,ISYDSG)
          NDWNC=NOW(2,ISYDC,MVSGM)
          IF(NDWNC.EQ.0) Cycle
          ISGSTA=1+IOCSF(ISYUSG,MVSGM,ISYSGM)
          NUPSG=NOW(1,ISYUSG,MVSGM)
          IOC=IOCSF(ISYUSG,MVSGM,ISYCI)
          INDEO=2*NLEV+(IQ*(IQ-1))/2+IP
          NCP=NOCP(INDEO,ISYDSG,MVSGM)
          IF(NCP.GT.0) THEN
            LICP=1+IOCP(INDEO,ISYDSG,MVSGM)
! CASE IS: LOWER HALF, DEEXCITE:
            CALL DEX1 (CPQ,NUPSG,CI(IOC+1),SGM(ISGSTA),NCP,ICOUP(1,LICP),VTAB)
          END IF
        END DO
      END DO

  ELSE IF(IP.GT.MIDLEV) THEN

! DEEXCITING OPERATOR, MIDLEV<IP<IQ
      DO MVSGM=1,NMIDV
        DO ISYUSG=1,NSYM
          NS1=NOCSF(ISYUSG,MVSGM,ISYSGM)
          IF(NS1.EQ.0) Cycle
          ISYUC=MUL(ISYPQ,ISYUSG)
          NUPC=NOW(1,ISYUC,MVSGM)
          IF (NUPC.EQ.0) Cycle
          ISGSTA=1+IOCSF(ISYUSG,MVSGM,ISYSGM)
          NUPSG=NOW(1,ISYUSG,MVSGM)
          ISYDSG=MUL(ISYUSG,ISYSGM)
          NDWNSG=NOW(2,ISYDSG,MVSGM)
          IOC=IOCSF(ISYUC,MVSGM,ISYCI)
          INDEO=2*NLEV+(IQ*(IQ-1))/2+IP
          NCP=NOCP(INDEO,ISYUSG,MVSGM)
          IF(NCP.GT.0) THEN
            LICP=1+IOCP(INDEO,ISYUSG,MVSGM)
! CASE IS: UPPER HALF, DEEXCITE:
            CALL DEX2 (CPQ,NDWNSG,NUPC,CI(IOC+1),NUPSG,SGM(ISGSTA),NCP,ICOUP(1,LICP),VTAB)
          END IF
        END DO
      END DO

  ELse

! DEEXCITING CASE, IP<=MIDLEV<IQ.
! ALLOCATE TEMPORARY WORK AREA:
      DO MVSGM=1,NMIDV
        MV4=MVR(MVSGM,1)
        MV5=MVR(MVSGM,2)
        DO ISYUSG=1,NSYM
          NS1=NOCSF(ISYUSG,MVSGM,ISYSGM)
          IF(NS1.EQ.0) Cycle
          ISGSTA=1+IOCSF(ISYUSG,MVSGM,ISYSGM)
          NUPSG=NOW(1,ISYUSG,MVSGM)
          ISYDSG=MUL(ISYUSG,ISYSGM)
          ISYUC=MUL(ISYQ,ISYUSG)
          ISYDC=MUL(ISYP,ISYDSG)
          IF(MV4.EQ.0) GOTO 499
          NUPC=NOW(1,ISYUC,MV4)
          IF(NUPC.EQ.0) GOTO 499
          NDWNC=NOW(2,ISYDC,MV4)
          IF(NDWNC.EQ.0) GOTO 499
          INDEO=IQ
          NCP=NOCP(INDEO,ISYUSG,MVSGM)
          IF(NCP.EQ.0) GOTO 499
          NTMP=NUPSG*NDWNC
          SGTMP(1:NTMP)=0.0D0
          LICP=1+IOCP(INDEO,ISYUSG,MVSGM)
          IOC=IOCSF(ISYUC,MV4,ISYCI)
! CASE IS: UPPER HALF, DEEXCITE:
          CALL DEX2 (CPQ,NDWNC,NUPC,CI(IOC+1),NUPSG,SGTMP,NCP,ICOUP(1,LICP),VTAB)
          INDEO=IP
          NCP=NOCP(INDEO,ISYDSG,MVSGM)
          IF(NCP.GT.0) THEN
            LICP=1+IOCP(INDEO,ISYDSG,MVSGM)
! CASE IS: LOWER HALF, DEEXCITE:
            X=1.0D00
            CALL DEX1 (X,NUPSG,SGTMP,SGM(ISGSTA),NCP,ICOUP(1,LICP),VTAB)
          END IF
 499      CONTINUE
          IF(MV5.EQ.0) Cycle
          NUPC=NOW(1,ISYUC,MV5)
          IF(NUPC.EQ.0) Cycle
          NDWNC=NOW(2,ISYDC,MV5)
          IF(NDWNC.EQ.0) Cycle
          INDEO=NLEV+IQ
          NCP=NOCP(INDEO,ISYUSG,MVSGM)
          IF(NCP.EQ.0) Cycle
          NTMP=NUPSG*NDWNC
          SGTMP(1:NTMP)=0.0D0
          LICP=1+IOCP(INDEO,ISYUSG,MVSGM)
          IOC=IOCSF(ISYUC,MV5,ISYCI)
! CASE IS: UPPER HALF, DEEXCITE:
          CALL DEX2 (CPQ,NDWNC,NUPC,CI(IOC+1),NUPSG,SGTMP,NCP,ICOUP(1,LICP),VTAB)
          INDEO=NLEV+IP
          NCP=NOCP(INDEO,ISYDSG,MVSGM)
          IF(NCP.GT.0) THEN
            LICP=1+IOCP(INDEO,ISYDSG,MVSGM)
! CASE IS: LOWER HALF, DEEXCITE:
            X=1.0D00
            CALL DEX1 (X,NUPSG,SGTMP,SGM(ISGSTA),NCP,ICOUP(1,LICP),VTAB)
          END IF
        END DO
      END DO

  END IF
ELSE
! THEN THIS IS A SPECIAL CASE.
  IF(IP.GT.MIDLEV) THEN

! SPECIAL CASE: WEIGHT OPERATOR, IP=IQ.
! IP=IQ>MIDLEV
      DO MVSGM=1,NMIDV
        DO ISYUSG=1,NSYM
          NS1=NOCSF(ISYUSG,MVSGM,ISYSGM)
          IF(NS1.EQ.0) CYCLE
          ISGSTA=1+IOCSF(ISYUSG,MVSGM,ISYSGM)
          NUPSG=NOW(1,ISYUSG,MVSGM)
          ISYDSG=MUL(ISYUSG,ISYSGM)
          NDWNSG=NOW(2,ISYDSG,MVSGM)
          IOUW=IOW(1,ISYUSG,MVSGM)
          IPSHFT=2*(IP-1-MIDLEV)
          LUW=1+IOUW-NIPWLK+IPSHFT/30
          IPSHFT=MOD(IPSHFT,30)
          IPPOW=2**IPSHFT
          DO I=1,NUPSG
            IC=ICase(LUW+I*NIPWLK)
            ICS=MOD(IC/IPPOW,4)
            IF(ICS.EQ.0) cycle
            X=CPQ*DBLE((1+ICS)/2)
            ISTA=ISGSTA-1+I
            CALL DAXPY_(NDWNSG,X,CI(ISTA),NUPSG,SGM(ISTA),NUPSG)
          END DO
        END DO
      END DO

  ELSE

! SPECIAL CASE: WEIGHT OPERATOR, IP=IQ.
! IP=IQ < MIDLEV.
      DO MVSGM=1,NMIDV
        DO ISYUSG=1,NSYM
          NS1=NOCSF(ISYUSG,MVSGM,ISYSGM)
          IF(NS1.EQ.0) Cycle
          ISGSTA=1+IOCSF(ISYUSG,MVSGM,ISYSGM)
          NUPSG=NOW(1,ISYUSG,MVSGM)
          ISYDSG=MUL(ISYUSG,ISYSGM)
          NDWNSG=NOW(2,ISYDSG,MVSGM)
          IOLW=IOW(2,ISYDSG,MVSGM)
          IPSHFT=2*(IP-1)
          LLW=1+IOLW-NIPWLK+IPSHFT/30
          IPSHFT=MOD(IPSHFT,30)
          IPPOW=2**IPSHFT
          DO J=1,NDWNSG
            JC=ICase(LLW+J*NIPWLK)
            ICS=MOD(JC/IPPOW,4)
            IF(ICS.EQ.0) Cycle
            X=CPQ*DBLE((1+ICS)/2)
            JSTA=ISGSTA+NUPSG*(J-1)
            CALL DAXPY_(NUPSG,X,CI(JSTA),1,SGM(JSTA),1)
          END DO
        END DO
      END DO

  END IF
END IF

      End Associate

contains
      SUBROUTINE EXC1 (CPQ,NUP,A,B,NCP,ICOUP,VTAB)
      IMPLICIT NONE
      INTEGER NUP,NCP,ICOUP(3,NCP)
      REAL*8 CPQ,A(NUP,*),B(NUP,*),VTAB(*)
      INTEGER ICP,JLFT,JRGT,I
      REAL*8 X
! CASE: ADD EPQ*A TO B, WHERE Q<P<=MIDLEV AND A AND B ARE SINGLE
! MATRIX BLOCKS.
! ONLY LOWER WALKS AFFECTED => COLUMN OPERATION.
! P>Q, EXCITING OPERATOR  => A IS  LEFTHAND WALK.
! USE VECTOR ROUTINE CALL IF LONG ENOUGH:
      IF(NUP.GT.20) THEN
        DO ICP=1,NCP
          JLFT=ICOUP(1,ICP)
          JRGT=ICOUP(2,ICP)
          X=CPQ*VTAB(ICOUP(3,ICP))
          CALL DAXPY_(NUP,X,A(1,JLFT),1,B(1,JRGT),1)
        END DO
      ELSE
        DO ICP=1,NCP
          JLFT=ICOUP(1,ICP)
          JRGT=ICOUP(2,ICP)
          X=CPQ*VTAB(ICOUP(3,ICP))
          DO I=1,NUP
            B(I,JRGT)=B(I,JRGT)+X*A(I,JLFT)
          END DO
        END DO
      END IF

      END SUBROUTINE EXC1

      SUBROUTINE EXC2 (CPQ,NDWN,NUPA,A,NUPB,B,NCP,ICOUP,VTAB)
      IMPLICIT NONE
      INTEGER NDWN,NUPA,NUPB,NCP,ICOUP(3,NCP)
      REAL*8 CPQ,A(NUPA,NDWN),B(NUPB,NDWN),VTAB(*)
      INTEGER ICP,ILFT,IRGT,J
      REAL*8 X
! CASE: ADD EPQ*A TO B, WHERE MIDLEV<Q<P AND A AND B ARE SINGLE
! MATRIX BLOCKS.
! ONLY UPPER WALKS AFFECTED =>    ROW OPERATION.
! Q<P, EXCITING OPERATOR  => A IS  LEFTHAND WALK.
      IF(NDWN.GT.15) THEN
        DO ICP=1,NCP
          ILFT=ICOUP(1,ICP)
          IRGT=ICOUP(2,ICP)
          X=CPQ*VTAB(ICOUP(3,ICP))
          CALL DAXPY_(NDWN,X,A(ILFT,1),NUPA,B(IRGT,1),NUPB)
        END DO
      ELSE
        DO ICP=1,NCP
          ILFT=ICOUP(1,ICP)
          IRGT=ICOUP(2,ICP)
          X=CPQ*VTAB(ICOUP(3,ICP))
          DO J=1,NDWN
            B(IRGT,J)=B(IRGT,J)+X*A(ILFT,J)
          END DO
        END DO
      END IF

      END SUBROUTINE EXC2

      SUBROUTINE DEX1(CPQ,NUP,A,B,NCP,ICOUP,VTAB)
      IMPLICIT NONE
      INTEGER NUP,NCP,ICOUP(3,NCP)
      REAL*8 CPQ,A(NUP,*),B(NUP,*),VTAB(*)
      INTEGER ICP,JLFT,JRGT,I
      REAL*8 X
! CASE: ADD EPQ*A TO B, WHERE P<Q<=MIDLEV AND A AND B ARE SINGLE
! MATRIX BLOCKS.
! ONLY LOWER WALKS AFFECTED => COLUMN OPERATION.
! P<Q, DEEXCITING OPERATOR  => A IS RIGHTHAND WALK.
! CHOSE DAXPY IF LARGE VECTOR LENGTH:
      IF(NUP.GT.20) THEN
        DO ICP=1,NCP
          JLFT=ICOUP(1,ICP)
          JRGT=ICOUP(2,ICP)
          X=CPQ*VTAB(ICOUP(3,ICP))
          CALL DAXPY_(NUP,X,A(1,JRGT),1,B(1,JLFT),1)
        END DO
      ELSE
        DO ICP=1,NCP
          JLFT=ICOUP(1,ICP)
          JRGT=ICOUP(2,ICP)
          X=CPQ*VTAB(ICOUP(3,ICP))
          DO I=1,NUP
            B(I,JLFT)=B(I,JLFT)+X*A(I,JRGT)
          END DO
        END DO
      END IF

      END SUBROUTINE DEX1

      SUBROUTINE DEX2 (CPQ,NDWN,NUPA,A,NUPB,B,NCP,ICOUP,VTAB)
      IMPLICIT NONE
      INTEGER NDWN,NUPA,NUPB,NCP,ICOUP(3,NCP)
      REAL*8 CPQ,A(NUPA,NDWN),B(NUPB,NDWN),VTAB(*)
      INTEGER ICP,ILFT,IRGT,J
      REAL*8 X
! CASE: ADD EPQ*A TO B, WHERE MIDLEV<P<Q AND A AND B ARE SINGLE
! MATRIX BLOCKS.
! ONLY UPPER WALKS AFFECTED =>    ROW OPERATION.
! P<Q, DEEXCITING OPERATOR  => A IS RIGHTHAND WALK.
! CHOSE DAXPY IF LARGE LENGTH.
      IF(NDWN.GT.10) THEN
        DO ICP=1,NCP
          ILFT=ICOUP(1,ICP)
          IRGT=ICOUP(2,ICP)
          X=CPQ*VTAB(ICOUP(3,ICP))
          CALL DAXPY_(NDWN,X,A(IRGT,1),NUPA,B(ILFT,1),NUPB)
        END DO
      ELSE
        DO ICP=1,NCP
          ILFT=ICOUP(1,ICP)
          IRGT=ICOUP(2,ICP)
          X=CPQ*VTAB(ICOUP(3,ICP))
          DO J=1,NDWN
            B(ILFT,J)=B(ILFT,J)+X*A(IRGT,J)
          END DO
        END DO
      END IF

      END SUBROUTINE DEX2

      END SUBROUTINE SIGMA1
