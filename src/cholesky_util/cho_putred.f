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
      SUBROUTINE CHO_PUTRED(IPASS,IRED)
C
C     Purpose: write reduced set indices to disk and set address for
C              next write.
C
      use ChoArr, only: iSP2F
      use ChoSwp, only: nnBstRSh, IndRSh, InfRed, IndRed
#include "implicit.fh"
#include "cholesky.fh"

      CHARACTER*10 SECNAM
      PARAMETER (SECNAM = 'CHO_PUTRED')

      IF (IPASS .GT. MAXRED) THEN
         WRITE(LUPRI,*) SECNAM,': integral pass ',IPASS
         WRITE(LUPRI,*) SECNAM,': max. allowed is ',MAXRED
         WRITE(LUPRI,*) SECNAM,': please increase max. allowed!'
         CALL CHO_QUIT('Too many integral passes in '//SECNAM,104)
      ELSE IF (IPASS .EQ. 1) THEN
         CALL CHO_PUTRED1(INFRED,nnBstRSh,IndRed(1,1),INDRSH,iSP2F,
     &                    MAXRED,NSYM,NNSHL,MMBSTRT,IPASS,1)
         IF (MAXRED .GT. 1) THEN
            INFRED(IPASS+1) = INFRED(IPASS)
     &                      + NSYM*NNSHL + 2*NNBSTRT(1) + NNSHL
         END IF
      ELSE IF (IPASS .EQ. MAXRED) THEN
         CALL CHO_PUTRED1(INFRED,nnBstRSh,IndRed(1,IRED),INDRSH,iSP2F,
     &                    MAXRED,NSYM,NNSHL,MMBSTRT,IPASS,IRED)
      ELSE
         CALL CHO_PUTRED1(INFRED,nnBstRSh,IndRed(1,IRED),INDRSH,iSP2F,
     &                    MAXRED,NSYM,NNSHL,MMBSTRT,IPASS,IRED)
         INFRED(IPASS+1) = INFRED(IPASS)
     &                   + NSYM*NNSHL + NNBSTRT(IRED)
      END IF

      END
