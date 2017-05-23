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
      SUBROUTINE MKSGNUM
     &           (LSYM,NSYM,NLEV,NVERT,MIDLEV,NMIDV,
     &            MXUP,MXDWN,NICASE,NIPWLK,
     &            IDOWN,IUP,IDAW,IRAW,NOW,IOW,
     &            IUSGNUM,ILSGNUM,ICASE,IPRINT)
C
C     PURPOSE: FOR ALL UPPER AND LOWER WALKS
C              COMPUTE THE DIRECT ARC WEIGHT SUM AND THE
C              REVERSE ARC WEIGHT SUM, RESPECTIVELY.
C              STORE THE DATA IN THE TABLES IUSGNUM AND ILSGNUM
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      DIMENSION IDOWN(NVERT,0:3),IUP(NVERT,0:3)
      DIMENSION IDAW(NVERT,0:4),IRAW(NVERT,0:4)
      DIMENSION NOW(2,NSYM,NMIDV),IOW(2,NSYM,NMIDV)
      DIMENSION IUSGNUM(MXUP,NMIDV),ILSGNUM(MXDWN,NMIDV)
      DIMENSION ICASE(NICASE)
      DIMENSION ISTEPVEC(50)
C
C
C     INITIALIZE NUMBERING TABLES
C
      DO 10 MIDV=1,NMIDV
        DO 20 IUW=1,MXUP
           IUSGNUM(IUW,MIDV)=0
20      CONTINUE
        DO 30 ILW=1,MXDWN
           ILSGNUM(ILW,MIDV)=0
30      CONTINUE
10    CONTINUE
C
C     MAIN LOOP RUNS OVER MIDVERTICES AND SYMMETRIES
C
      ICONF=0
      DO 100 MIDV=1,NMIDV
        DO 110 ISYM=1,NSYM
          IUOFF=1+IOW(1,ISYM,MIDV)
          NUW=NOW(1,ISYM,MIDV)
          JSYM=1+IEOR(ISYM-1,LSYM-1)
          ILOFF=1+IOW(2,JSYM,MIDV)
          NLW=NOW(2,JSYM,MIDV)
          IF( NUW.EQ.0 .OR. NLW.EQ.0 ) GOTO 110
C
C         LOOP OVER ALL UPPER WALKS
C
          DO 200 IUW=1,NUW
            IPOS=IUOFF+NIPWLK*(IUW-1)
C     UNPACK THE UPPER WALK STEP VECTOR
            ICODE=ICASE(IPOS)
            JPOS=0
            DO 210 LEV=(MIDLEV+1),NLEV
              JPOS=JPOS+1
              IF( JPOS.EQ.16 ) THEN
                JPOS=1
                IPOS=IPOS+1
                ICODE=ICASE(IPOS)
              ENDIF
              ISTEP=MOD(ICODE,4)
              ISTEPVEC(LEV)=ISTEP
              ICODE=ICODE/4
210         CONTINUE
C     GET REVERSE ARC WEIGHT FOR UPPER WALK
            IRAWSUM=1
            LV=1
            DO 220 LEV=NLEV,(MIDLEV+1),-1
              IC=ISTEPVEC(LEV)
              LV=IDOWN(LV,IC)
              IRAWSUM=IRAWSUM+IRAW(LV,IC)
220         CONTINUE
            IUSGNUM(IRAWSUM,MIDV)=IUW
200       CONTINUE
C
C         LOOP OVER ALL LOWER WALKS
C
          DO 300 ILW=1,NLW
            IPOS=ILOFF+NIPWLK*(ILW-1)
C     UNPACK WALK STEP VECTOR
            ICODE=ICASE(IPOS)
            JPOS=0
            DO 310 LEV=1,MIDLEV
              JPOS=JPOS+1
              IF( JPOS.EQ.16 ) THEN
                JPOS=1
                IPOS=IPOS+1
                ICODE=ICASE(IPOS)
              ENDIF
              ISTEP=MOD(ICODE,4)
              ISTEPVEC(LEV)=ISTEP
              ICODE=ICODE/4
310         CONTINUE
C     GET DIRECT ARC WEIGHT FOR THE LOWER WALK
            IDAWSUM=1
            LV=NVERT
            DO 320 LEV=1,MIDLEV
              IC=ISTEPVEC(LEV)
              LV=IUP(LV,IC)
              IDAWSUM=IDAWSUM+IDAW(LV,IC)
320         CONTINUE
            ILSGNUM(IDAWSUM,MIDV)=ICONF
            ICONF=ICONF+NUW
300       CONTINUE
C
110     CONTINUE
100   CONTINUE

      IF( IPRINT.GT.5 ) THEN
        WRITE(6,*)
        WRITE(6,*)' ILSGNUM IN SUBROUTINE MKSGNUM'
        DO 400 MIDV=1,NMIDV
          WRITE(6,'(1X,''MIDV='',I3,/,(20I6))')MIDV,
     *         (ILSGNUM(J,MIDV),J=1,MXDWN)
400     CONTINUE
        WRITE(6,*)
        WRITE(6,*)' IUSGNUM IN SUBROUTINE MKSGNUM'
        DO 410 MIDV=1,NMIDV
          WRITE(6,'(1X,''MIDV='',I3,/,(20I6))')MIDV,
     *          (IUSGNUM(J,MIDV),J=1,MXUP)
410     CONTINUE
        WRITE(6,*)
        WRITE(6,*)
        WRITE(6,*)' NR OF WALKS AND CONFIGURATIONS IN SGNUM'
        WRITE(6,*)' BY MIDVERTEX AND SYMMETRY.'
        DO 810 MV=1,NMIDV
          WRITE(6,*)
          WRITE(6,1234) MV,(NOW(1,IS,MV),IS=1,NSYM)
          WRITE(6,1235)    (NOW(2,IS,MV),IS=1,NSYM)
1234  FORMAT('  MV=',I2,'    UPPER WALKS:',8I6)
1235  FORMAT('           LOWER WALKS:',8I6)
810     CONTINUE
        WRITE(6,*)' OFFSETS IN MKSGNUM'
        WRITE(6,*)' BY MIDVERTEX AND SYMMETRY.'
        DO 820 MV=1,NMIDV
          WRITE(6,*)
          WRITE(6,1234) MV,(IOW(1,IS,MV),IS=1,NSYM)
          WRITE(6,1235)    (IOW(2,IS,MV),IS=1,NSYM)
820     CONTINUE
      ENDIF
C
C
C     EXIT
C
      RETURN
      END
