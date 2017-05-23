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
      SUBROUTINE CHO_RDDBUF(DIAG,BUF,IBUF,INDRSH,INDRED,
     &                      LENBUF,LMMBSTRT,NDUMP)
C
C     Purpose: read diagonal from disk and set first reduced set
C              indices.
C
#include "implicit.fh"
      DIMENSION DIAG(*), BUF(LENBUF)
      INTEGER   IBUF(4,LENBUF)
      INTEGER   INDRSH(LMMBSTRT), INDRED(LMMBSTRT,3)
#include "cholesky.fh"
#include "choptr.fh"
#include "WrkSpc.fh"

      CHARACTER*10 SECNAM
      PARAMETER (SECNAM = 'CHO_RDDBUF')

      IIBSTRSH(I,J,K)=IWORK(ip_IIBSTRSH-1+NSYM*NNSHL*(K-1)+NSYM*(J-1)+I)
      ISP2F(I)=IWORK(ip_iSP2F-1+I)

      IF (LENBUF .LT. LBUF) THEN
         WRITE(LUPRI,'(//,1X,A,A)') SECNAM,
     &                              ': LENBUF >= LBUF required!'
         WRITE(LUPRI,'(1X,A,I10)')    'LENBUF = ',LENBUF
         WRITE(LUPRI,'(1X,A,I10,/)')  'LBUF   = ',LBUF
         CALL CHO_QUIT('Buffer error in '//SECNAM,102)
      END IF

      IUNIT = LUSCR
      LUSCR = -1
      REWIND(IUNIT)

      DO IDUMP = 1,NDUMP
         CALL CHO_RDBUF(LENGTH,BUF,IBUF,LBUF,IUNIT)
         IF (IDUMP .EQ. NDUMP) THEN
             CALL CHO_CLOSE(IUNIT,'DELETE')
         END IF
         DO L = 1,LENGTH
            IF (IBUF(2,L) .GT. 0) THEN
               ISHLAB = IBUF(1,L)
               ISYMAB = IBUF(3,L)
               IAB    = IIBSTR(ISYMAB,1) + IIBSTRSH(ISYMAB,ISHLAB,1)
     &                + IBUF(2,L)
               DIAG(IAB) = BUF(L)
               INDRSH(IAB) = ISP2F(ISHLAB)
               INDRED(IAB,1) = IBUF(4,L)
            END IF
         END DO
      END DO

      END
