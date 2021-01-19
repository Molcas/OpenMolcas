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
      SUBROUTINE CHO_GETDIAG1(DIAG,BUF,IBUF,LENBUF,NDUMP)
C
C     Purpose: read diagonal in first reduced set.
C
      use ChoSwp, only: nnBstRSh, iiBstRSh, IndRSh, IndRed
#include "implicit.fh"
      DIMENSION DIAG(*), BUF(LENBUF)
      INTEGER   IBUF(4,LENBUF)
#include "cholesky.fh"
#include "choprint.fh"

      CHARACTER*12 SECNAM
      PARAMETER (SECNAM = 'CHO_GETDIAG1')

      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)

      INTEGER ISYLST(8)

      PARAMETER (INFOD = INF_DIAG)
      PARAMETER (TINY  = 1.0D-14)

C     Read diagonal from file.
C     ------------------------

      IF (RSTDIA) THEN
         IOPT = 2
         CALL CHO_IODIAG(DIAG,IOPT)
      ELSE
         CALL CHO_DZERO(DIAG,NNBSTRT(1))
         CALL CHO_IZERO(INDRSH,NNBSTRT(1))
         CALL CHO_IZERO(INDRED,NNBSTRT(1))
         CALL CHO_RDDBUF(DIAG,BUF,IBUF,INDRSH,INDRED,
     &                   LENBUF,MMBSTRT,NDUMP)
         CALL CHO_GADGOP(DIAG,NNBSTRT(1),'+')
         CALL CHO_GAIGOP(INDRSH,NNBSTRT(1),'+')
         CALL CHO_GAIGOP(INDRED,NNBSTRT(1),'+')
      END IF

C     Copy info to current reduced set (IRED=2).
C     Also set up IRED=3 (although it should be redundant).
C     -----------------------------------------------------

      DO IRS = 2,3
         CALL CHO_RSCOPY(IIBSTRSH,NNBSTRSH,
     &                   INDRED,1,IRS,NSYM,NNSHL,NNBSTRT(1),3)
      END DO

C     Print.
C     ------

      IF (LOCDBG .OR. (IPRINT.GE.INFOD)) THEN
         DO ISYM = 1,NSYM
            ISYLST(ISYM) = ISYM
         END DO
         NSYLST = NSYM
         IRED = 1
         CALL CHO_PRTDIA(DIAG,ISYLST,NSYLST,IRED)
         IF (LOCDBG) THEN
            IRED = 2
            CALL CHO_PRTDIA(DIAG,ISYLST,NSYLST,IRED)
         END IF
      END IF

      END
