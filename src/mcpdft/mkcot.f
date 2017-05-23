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
      SUBROUTINE MKCOT_m(ISM,IDOWN,NOW,IOW,IOCSF,NOCSF,ISCR,IPRINT)
C     PURPOSE: SET UP COUNTER AND OFFSET TABLES FOR WALKS AND CSFS
C     NOTE:    TO GET GET VARIOUS COUNTER AND OFFSET TABLES
C              THE DOWN-CHAIN TABLE IS SCANNED TO PRODUCE ALL POSSIBLE
C              WALKS. POSSIBLY, THERE ARE MORE EFFICIENT WAYS, BUT
C              SINCE ONLY UPPER AND LOWER WALKS ARE REQUIRED
C              THEIR NUMBER IS VERY LIMITTED, EVEN FOR LARGE CASES.
      IMPLICIT REAL*8 (A-H,O-Z)
C
#include "rasdim.fh"
#include "general.fh"
#include "output_ras.fh"
      Parameter (ROUTINE='MKCOT   ')
#include "gugx.fh"
C
      DIMENSION ISM(NLEV),IDOWN(NVERT,0:3)
      DIMENSION NOW(2,NSYM,NMIDV),IOW(2,NSYM,NMIDV)
      DIMENSION NOCSF(NSYM,NMIDV,NSYM),IOCSF(NSYM,NMIDV,NSYM)
      DIMENSION ISCR(3,0:NLEV)
      PARAMETER(IVERT=1,ISYM=2,ISTEP=3)
C
C     CLEAR ARRAYS IOW AND NOW
C
      DO IHALF=1,2
        DO MV=1,NMIDV
          DO IS=1,NSYM
            NOW(IHALF,IS,MV)=0
            IOW(IHALF,IS,MV)=0
          END DO
        END DO
      END DO
C
C     CLEAR ARRAYS IOCSF AND NOCSF
C
      DO IS=1,NSYM
        DO MV=1,NMIDV
          DO JS=1,NSYM
            IOCSF(JS,MV,IS)=0
            NOCSF(JS,MV,IS)=0
          END DO
        END DO
      END DO
C
C     START MAIN LOOP OVER UPPER AND LOWER WALKS, RESPECTIVELY.
C
      DO IHALF=1,2
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
        DO IVTOP=IVTSTA,IVTEND
C     SET CURRENT LEVEL=TOP LEVEL OF SUBGRAPH
          LEV=LEV1
          ISCR(IVERT,LEV)=IVTOP
          ISCR(ISYM,LEV)=1
          ISCR(ISTEP,LEV)=-1
100       IF(LEV.GT.LEV1) GOTO 400
C     FIND FIRST POSSIBLE UNTRIED ARC DOWN FROM CURRENT VERTEX
          IVT=ISCR(IVERT,LEV)
          DO ISTP=ISCR(ISTEP,LEV)+1,3
            IVB=IDOWN(IVT,ISTP)
            IF(IVB.NE.0) GOTO 200
          END DO
C     NO SUCH ARC WAS POSSIBLE. GO UP ONE STEP AND TRY AGAIN.
          ISCR(ISTEP,LEV)=-1
          LEV=LEV+1
          GOTO 100
C     SUCH AN ARC WAS FOUND. WALK DOWN:
200       ISCR(ISTEP,LEV)=ISTP
          ISML=1
          IF((ISTP.EQ.1).OR.(ISTP.EQ.2)) ISML=ISM(LEV)
          LEV=LEV-1
          ISCR(ISYM,LEV)=MUL(ISML,ISCR(ISYM,LEV+1))
          ISCR(IVERT,LEV)=IVB
          ISCR(ISTEP,LEV)=-1
          IF (LEV.GT.LEV2) GOTO 100
C     WE HAVE REACHED THE BOTTOM LEVEL. THE WALK IS COMPLETE.
C     FIND MIDVERTEX NUMBER ORDERING NUMBER AND SYMMETRY OF THIS WALK
          MV=ISCR(IVERT,MIDLEV)+1-MIDV1
          IWSYM=ISCR(ISYM,LEV2)
          ILND=1+NOW(IHALF,IWSYM,MV)
C     SAVE THE MAX WALK NUMBER FOR GIVEN SYMMETRY AND MIDVERTEX
          NOW(IHALF,IWSYM,MV)=ILND
C     BACK UP ONE LEVEL AND TRY AGAIN:
          LEV=LEV+1
          GOTO 100
400       CONTINUE
        END DO
      END DO
C
C     NOW,CONSTRUCT OFFSET TABLES FOR UPPER AND LOWER WALKS
C     SEPARATED FOR EACH MIDVERTEX AND SYMMETRY
C
      NUW=0
      DO MV=1,NMIDV
        DO IS=1,NSYM
          IOW(1,IS,MV)=NUW*NIPWLK
          NUW=NUW+NOW(1,IS,MV)
        END DO
      END DO
      NWALK=NUW
      DO MV=1,NMIDV
        DO IS=1,NSYM
          IOW(2,IS,MV)=NWALK*NIPWLK
          NWALK=NWALK+NOW(2,IS,MV)
        END DO
      END DO
      NLW=NWALK-NUW
C
C     FINALLY, CONSTRUCT COUNTER AND OFFSET TABLES FOR THE CSFS
C     SEPARATED BY MIDVERTICES AND SYMMETRY.
C     FORM ALSO CONTRACTED SUMS OVER MIDVERTICES.
C
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
      IF (IPRINT.GE.5) THEN
        Write(LF,*)
        Write(LF,*)' TOTAL NR OF WALKS: UPPER ',NUW
        Write(LF,*)'                    LOWER ',NLW
        Write(LF,*)'                     SUM  ',NWALK
        Write(LF,*)
        Write(LF,*)' NR OF CONFIGURATIONS/SYMM:'
        Write(LF,'(8(1X,I8))')(NCSF(IS),IS=1,NSYM)
        Write(LF,*)
        Write(LF,*)
        Write(LF,*)' NR OF WALKS AND CONFIGURATIONS IN NRCOUP'
        Write(LF,*)' BY MIDVERTEX AND SYMMETRY.'
        DO MV=1,NMIDV
          Write(LF,*)
          Write(LF,'(A,I2,A,8I6)') '  MV=',MV,'    UPPER WALKS:',
     &                             (NOW(1,IS,MV),IS=1,NSYM)
          Write(LF,'(A,8I6)') '           LOWER WALKS:',
     &                             (NOW(2,IS,MV),IS=1,NSYM)
          DO IST=1,NSYM
          Write(LF,'(A,I2,A,8I6)') ' IST=',IST,'  CONFIGURATIONS:',
     &                           (NOCSF(IS,MV,IST),IS=1,NSYM)
          END DO
        END DO
      ENDIF

      RETURN
      END
