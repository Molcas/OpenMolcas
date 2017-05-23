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
      SUBROUTINE GETDPREF(DREF,PREF)
      IMPLICIT REAL*8 (A-H,O-Z)
#include "rasdim.fh"
#include "caspt2.fh"
#include "pt2_guga.fh"
#include "output.fh"
#include "WrkSpc.fh"
#include "intgrl.fh"
      REAL*8 DREF(NDREF)
      REAL*8 PREF(NPREF)

      CALL QENTER('GETDPREF')
* Get active 1-density and 2-density matrices GAMMA1 and
* GAMMA2, and construct DREF and PREF which are in a tringular
* storage.

C Remember: NDREF=1 if NASHT=0. Similar NPREF.
      DREF(1)=0.0d0
      PREF(1)=0.0d0
      IF(NASHT.EQ.0) GO TO 99
c Active density 1-matrix:
      CALL GETMEM('LG1','ALLO','REAL',LG1,NG1)
      CALL PT2_GET(NG1,'GAMMA1',WORK(LG1))
      DO I=1,NASHT
       DO J=1,I
        IJ=(I*(I-1))/2+J
        DREF(IJ)=WORK(LG1-1+I+NASHT*(J-1))
       END DO
      END DO
      CALL GETMEM('LG1','FREE','REAL',LG1,NG1)
C CONSTRUCT PREF, 2-ELECTRON DENSITY MATRIX:
      CALL GETMEM('LG2','ALLO','REAL',LG2,NG2)
      CALL PT2_GET(NG2,'GAMMA2',WORK(LG2))
      IJT=0
      IJKLT=0
      N2=NASHT**2
      DO I=1,NASHT
        DO J=1,I
          IJT=IJT+1
          IJ=I+NASHT*(J-1)
          JI=J+NASHT*(I-1)
          KLT=0
          DO K=1,NASHT
            DO L=1,K
              KLT=KLT+1
              IF(KLT.GT.IJT) GOTO 130
              IJKLT=IJKLT+1
              KL=K+NASHT*(L-1)
              LK=L+NASHT*(K-1)

              P1=0.5D0*WORK(LG2-1+IJ+N2*(KL-1))
              P2=0.5D0*WORK(LG2-1+IJ+N2*(LK-1))
              IF(IJ.GE.KL) THEN
                IJKL=(IJ*(IJ-1))/2+KL
              ELSE
                IJKL=(KL*(KL-1))/2+IJ
              END IF
              IF(IJ.GE.LK) THEN
                IJLK=(IJ*(IJ-1))/2+LK
              ELSE
                IJLK=(LK*(LK-1))/2+IJ
              END IF
              JIKL=(JI*(JI-1))/2+KL
              JILK=(JI*(JI-1))/2+LK
              PREF(IJKL)=P1
              PREF(IJLK)=P2
              PREF(JIKL)=P2
              PREF(JILK)=P1
            END DO
          END DO
 130    CONTINUE
        END DO
      END DO

      CALL GETMEM('LG2','FREE','REAL',LG2,NG2)

      IF(IPRGLB.GE.DEBUG) THEN
       WRITE(6,*)' GETDPREF has constructed DREF and PREF.'
       CALL XFLUSH(6)
      END IF

  99  CONTINUE

      CALL QEXIT('GETDPREF')
      RETURN
      END
