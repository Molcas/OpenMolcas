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
      SUBROUTINE MKCLIST_MCLR
     &          (NSYM,NLEV,NVERT,MIDLEV,MIDV1,MIDV2,NMIDV,NICASE,NIPWLK,
     &           ISM,IDOWN,NOW,IOW,ICASE,ISCR,IPRINT)
C
C     PURPOSE: CONSTRUCT THE COMPRESSED CASE-LIST, I.E.,
C              STORE THE STEP VECTOR FOR ALL POSSIBLE WALKS
C              IN THE ARRAY ICASE. GROUPS OF 15 CASES ARE PACKED
C              INTO ONE INTEGER WORD.
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      DIMENSION ISM(NLEV),IDOWN(NVERT,0:3)
      DIMENSION NOW(2,NSYM,NMIDV),IOW(2,NSYM,NMIDV)
      DIMENSION ICASE(NICASE)
      DIMENSION ISCR(3,0:NLEV)
      PARAMETER(IVERT=1,ISYM=2,ISTEP=3)
C
C
C     CLEAR ARRAY NOW. IT WILL BE RESTORED FINALLY
C
      iZero=0
      DO 10 IHALF=1,2
        DO 11 MV=1,NMIDV
          DO 12 IS=1,NSYM
            NOW(IHALF,IS,MV)=0
12        CONTINUE
11      CONTINUE
10    CONTINUE
C
C     START MAIN LOOP OVER UPPER AND LOWER WALKS, RESPECTIVELY.
C
      DO 500 IHALF=1,2
        IF(IHALF.EQ.1) THEN
          IVTSTA=1
          IVTEND=1
          LEV1=NLEV
          LEV2=MIDLEV
        ELSE
          IVTSTA=MIDV1
          IVTEND=MIDV2
          LEV1=MIDLEV
          LEV2=0
        END IF
C
C     LOOP OVER VERTICES STARTING AT TOP OF SUBGRAPH
C
        DO 400 IVTOP=IVTSTA,IVTEND
C SET CURRENT LEVEL=TOP LEVEL OF SUBGRAPH:
          LEV=LEV1
          ISCR(IVERT,LEV)=IVTOP
          ISCR(ISYM,LEV)=1
          ISCR(ISTEP,LEV)=-1
100       IF(LEV.GT.LEV1) GOTO 400
C FIND FIRST POSSIBLE UNTRIED ARC DOWN FROM CURRENT VERTEX:
          IVT=ISCR(IVERT,LEV)
          DO 110 ISTP=ISCR(ISTEP,LEV)+1,3
            IVB=IDOWN(IVT,ISTP)
            IF(IVB.EQ.0) GOTO 110
            GOTO 200
110       CONTINUE
C ALT A -- NO SUCH ARC WAS POSSIBLE. GO UP ONE STEP AND TRY AGAIN.
          ISCR(ISTEP,LEV)=-1
          LEV=LEV+1
          GOTO 100
C ALT B -- SUCH AN ARC WAS FOUND. WALK DOWN:
200       ISCR(ISTEP,LEV)=ISTP
          ISML=1
          IF((ISTP.EQ.1).OR.(ISTP.EQ.2)) ISML=ISM(LEV)
          LEV=LEV-1
CFUE      ISCR(ISYM,LEV)=MUL(ISML,ISCR(ISYM,LEV+1))
          ISCR(ISYM,LEV)=1+IEOR(ISML-1,ISCR(ISYM,LEV+1)-1)
          ISCR(IVERT,LEV)=IVB
          ISCR(ISTEP,LEV)=-1
          IF (LEV.GT.LEV2) GOTO 100
C WE HAVE REACHED THE LOWER LEVEL. THE WALK IS COMPLETE.
C MIDVERTEX NUMBER:
          MV=ISCR(IVERT,MIDLEV)+1-MIDV1
C SYMMETRY LABEL OF THIS WALK:
          IWSYM=ISCR(ISYM,LEV2)
C ITS ORDERING NUMBER WITHIN THE SAME BATCH OF (IHALF,IWSYM,MV):
          ILND=1+NOW(IHALF,IWSYM,MV)
          NOW(IHALF,IWSYM,MV)=ILND
C CONSEQUENTLY, THE POSITION IMMEDIATELY BEFORE THIS COMPRESSED WALK:
          IPOS=IOW(IHALF,IWSYM,MV)+(ILND-1)*NIPWLK
C PACK THE STEPS IN GROUPS OF 15 LEVELS PER INTEGER:
          DO 205 LL=LEV2+1,LEV1,15
            IC=0
            DO 204 L=MIN(LL+14,LEV1),LL,-1
              IC=4*IC+ISCR(ISTEP,L)
204         CONTINUE
            IPOS=IPOS+1
            ICASE(IPOS)=IC
205       CONTINUE
C FINISHED WITH THIS WALK. BACK UP ONE LEVEL AND TRY AGAIN:
          LEV=LEV+1
          GOTO 100
400     CONTINUE
500   CONTINUE
C
C
C     EXIT
C
      RETURN
c Avoid unused argument warnings
      IF (.FALSE.) CALL Unused_integer(IPRINT)
      END
