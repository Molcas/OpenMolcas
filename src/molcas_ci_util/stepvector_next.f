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
      SUBROUTINE STEPVECTOR_NEXT(MV,IDWN,IUP,
     &                           STEPVECTOR)
      IMPLICIT NONE
      INTEGER :: MV, IDWN, IUP, STEPVECTOR(*)
#include "rasdim.fh"
#include "rasscf.fh"
#include "general.fh"
#include "output_ras.fh"
      PARAMETER(ROUTINE='STEPVECTOR_NEXT')
#include "gugx.fh"
#include "WrkSpc.fh"

C stop when MV is zero
      IF (MV.EQ.0) THEN
        WRITE(6,'(1X,A)') 'stepvector_next has been depleted'
      END IF

      CALL GETSTEPVECTOR(IWORK(LNOW),IWORK(LIOW),
     &                   MV,IDWN,IUP,
     &                   STEPVECTOR)

      END

      SUBROUTINE GETSTEPVECTOR(NOW,IOW,
     &                         MV,IDWN,IUP,
     &                         ICS)
      IMPLICIT REAL*8 (A-H,O-Z)
C
#include "rasdim.fh"
#include "rasscf.fh"
#include "general.fh"
#include "output_ras.fh"
      PARAMETER(ROUTINE='GETSTEPVECTOR')
#include "gugx.fh"
#include "WrkSpc.fh"

      DIMENSION NOW(2,NSYM,NMIDV),IOW(2,NSYM,NMIDV)

      DIMENSION ICS(mxact)
C
C     RECONSTRUCT THE CASE LIST
C
      NICASE=NWALK*NIPWLK
*     NSCR=3*(NLEV+1)
*     CALL GETMEM('SCR1','ALLO','INTEG',LSCR,NSCR)
*     CALL GETMEM('CASE','ALLO','INTEG',LICASE,NICASE)
*     CALL MKCLIST(NSM,IWORK(LDOWN),IWORK(LNOW),IWORK(LIOW),
*    &             IWORK(LICASE),IWORK(LSCR))
*     CALL GETMEM('SCR1','FREE','INTEG',LSCR,NSCR)

C
C     ENTER THE MAIN LOOP IS OVER BLOCKS OF THE ARRAY CI
C     WITH SPECIFIED MIDVERTEX MV, AND UPPERWALK SYMMETRY ISYUP.
C
!     DO MV=1,NMIDV
!       DO ISYUP=1,NSYM
          NUP=NOW(1,1,MV)
          NDWN=NOW(2,1,MV)
          IUW0=LICASE-NIPWLK+IOW(1,1,MV)
          IDW0=LICASE-NIPWLK+IOW(2,1,MV)
*         IDWNSV=0
!         DO IDWN=1,NDWN
!           DO IUP=1,NUP
              ! determine the stepvector
*             IF(IDWNSV.NE.IDWN) THEN
                ICDPOS=IDW0+IDWN*NIPWLK
                ICDWN=IWORK(ICDPOS)
                ! unpack lower walk
                NNN=0
                DO LEV=1,MIDLEV
                  NNN=NNN+1
                  IF(NNN.EQ.16) THEN
                    NNN=1
                    ICDPOS=ICDPOS+1
                    ICDWN=IWORK(ICDPOS)
                  END IF
                  IC1=ICDWN/4
                  ICS(LEV)=ICDWN-4*IC1
                  ICDWN=IC1
                END DO
*               IDWNSV=IDWN
*             END IF
              ICUPOS=IUW0+NIPWLK*IUP
              ICUP=IWORK(ICUPOS)
              ! unpack upper walk
              NNN=0
              DO LEV=MIDLEV+1,NLEV
                NNN=NNN+1
                IF(NNN.EQ.16) THEN
                  NNN=1
                  ICUPOS=ICUPOS+1
                  ICUP=IWORK(ICUPOS)
                END IF
                IC1=ICUP/4
                ICS(LEV)=ICUP-4*IC1
                ICUP=IC1
              END DO
!           END DO
!         END DO
!       END DO
!     END DO

C compute the next set of indices
      IF (IUP.EQ.NUP) THEN
        IF (IDWN.EQ.NDWN) THEN
          IF (MV.EQ.NMIDV) THEN
            MV=0
          ELSE
            MV=MV+1
          END IF
          IDWN=1
        ELSE
          IDWN=IDWN+1
        END IF
        IUP=1
      ELSE
        IUP=IUP+1
      END IF

*     CALL GETMEM('CASE','FREE','INTEG',LICASE,NICASE)
      RETURN
      END SUBROUTINE
