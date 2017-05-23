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
      SUBROUTINE MKSGNUM_m(IDOWN,IUP,IDAW,IRAW,NOW,IOW,
     *                   IUSGNUM,ILSGNUM,ICASE,IPRINT)
C     PURPOSE: FOR ALL UPPER AND LOWER WALKS
C              COMPUTE THE DIRECT ARC WEIGHT SUM AND THE
C              REVERSE ARC WEIGHT SUM, RESPECTIVELY.
C              STORE THE DATA IN THE TABLES IUSGNUM AND ILSGNUM
C
      IMPLICIT REAL*8 (A-H,O-Z)
#include "rasdim.fh"
#include "general.fh"
#include "gugx.fh"
#include "WrkSpc.fh"
#include "output_ras.fh"
      Parameter(Routine='MKSGNUM')
C
      DIMENSION IDOWN(NVERT,0:3),IUP(NVERT,0:3)
      DIMENSION IDAW(NVERT,0:4),IRAW(NVERT,0:4)
      DIMENSION NOW(2,NSYM,NMIDV),IOW(2,NSYM,NMIDV)
      DIMENSION IUSGNUM(MXUP,NMIDV),ILSGNUM(MXDWN,NMIDV)
      DIMENSION ICASE(NICASE)
      DIMENSION ISTEPVEC(mxact)

      Call qEnter(Routine)
C
C     INITIALIZE NUMBERING TABLES
C
      DO MIDV=1,NMIDV
        DO IUW=1,MXUP
           IUSGNUM(IUW,MIDV)=0
        END DO
        DO ILW=1,MXDWN
           ILSGNUM(ILW,MIDV)=0
        END DO
      END DO
C
C     MAIN LOOP RUNS OVER MIDVERTICES AND SYMMETRIES
C
      ICONF=0
      DO MIDV=1,NMIDV
        DO ISYM=1,NSYM
          IUOFF=1+IOW(1,ISYM,MIDV)
          NUW=NOW(1,ISYM,MIDV)
          JSYM=MUL(ISYM,LSYM)
          ILOFF=1+IOW(2,JSYM,MIDV)
          NLW=NOW(2,JSYM,MIDV)
          IF( NUW.EQ.0 .OR. NLW.EQ.0 ) GOTO 110
C
C         LOOP OVER ALL UPPER WALKS
C
          DO IUW=1,NUW
            IPOS=IUOFF+NIPWLK*(IUW-1)
C     UNPACK THE UPPER WALK STEP VECTOR
            ICODE=ICASE(IPOS)
            JPOS=0
            DO LEV=(MIDLEV+1),NLEV
              JPOS=JPOS+1
              IF( JPOS.EQ.16 ) THEN
                JPOS=1
                IPOS=IPOS+1
                ICODE=ICASE(IPOS)
              ENDIF
              ISTEP=MOD(ICODE,4)
              ISTEPVEC(LEV)=ISTEP
              ICODE=ICODE/4
            END DO
C     GET REVERSE ARC WEIGHT FOR UPPER WALK
            IRAWSUM=1
            LV=1
            DO LEV=NLEV,(MIDLEV+1),-1
              IC=ISTEPVEC(LEV)
              LV=IDOWN(LV,IC)
              IRAWSUM=IRAWSUM+IRAW(LV,IC)
            END DO
            IUSGNUM(IRAWSUM,MIDV)=IUW
          END DO
C
C         LOOP OVER ALL LOWER WALKS
C
          DO ILW=1,NLW
            IPOS=ILOFF+NIPWLK*(ILW-1)
C     UNPACK WALK STEP VECTOR
            ICODE=ICASE(IPOS)
            JPOS=0
            DO LEV=1,MIDLEV
              JPOS=JPOS+1
              IF( JPOS.EQ.16 ) THEN
                JPOS=1
                IPOS=IPOS+1
                ICODE=ICASE(IPOS)
              ENDIF
              ISTEP=MOD(ICODE,4)
              ISTEPVEC(LEV)=ISTEP
              ICODE=ICODE/4
            END DO
C     GET DIRECT ARC WEIGHT FOR THE LOWER WALK
            IDAWSUM=1
            LV=NVERT
            DO LEV=1,MIDLEV
              IC=ISTEPVEC(LEV)
              LV=IUP(LV,IC)
              IDAWSUM=IDAWSUM+IDAW(LV,IC)
            END DO
            ILSGNUM(IDAWSUM,MIDV)=ICONF
            ICONF=ICONF+NUW
          END DO
C
110       CONTINUE
        END DO
      END DO
      IF( IPRINT.GT.5 ) THEN
        Write(LF,*)
        Write(LF,*)' ILSGNUM IN SUBROUTINE MKSGNUM'
        DO MIDV=1,NMIDV
          Write(LF,'(1X,''MIDV='',I3,/,(20I6))')MIDV,
     *         (ILSGNUM(J,MIDV),J=1,MXDWN)
        END DO
        Write(LF,*)
        Write(LF,*)' IUSGNUM IN SUBROUTINE MKSGNUM'
        DO MIDV=1,NMIDV
          Write(LF,'(1X,''MIDV='',I3,/,(20I6))')MIDV,
     *          (IUSGNUM(J,MIDV),J=1,MXUP)
        END DO
        Write(LF,*)
      ENDIF

      Call qExit(Routine)
      RETURN
      END
