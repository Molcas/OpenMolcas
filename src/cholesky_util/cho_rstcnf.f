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
      SUBROUTINE CHO_RSTCNF(NERR)
C
C     Purpose: check restart configuration info.
C
      IMPLICIT NONE
      INTEGER NERR
#include "cholesky.fh"

      INTEGER I, J

      CHARACTER*3 SWITCH(2)
      DATA SWITCH /' ON','OFF'/

      REAL*8 ERRTOL, ERR
      PARAMETER (ERRTOL = 1.0D-14)

      NERR = 0

      IF (CHO_ADRVEC .NE. XCHO_ADRVEC) THEN
         WRITE(LUPRI,'(A,I9,/,A,I9)')
     &   'RESTART: addressing mode for vectors from restart file:',
     &   XCHO_ADRVEC,
     &   '         addressing mode for vectors from input       :',
     &   CHO_ADRVEC
         WRITE(LUPRI,'(A,A)')
     &   '         Restart will fail - please specify correct address ',
     &   'mode.'
         CALL CHO_QUIT('Cholesky restart failure in CHO_RSTCNF',105)
      END IF

      ERR = ABS(THRCOM - XTHRCOM)
      IF (ERR .GT. ERRTOL) THEN
         WRITE(LUPRI,'(A,D16.8,/,A,D16.8)')
     &   'RESTART: decomposition threshold from restart file: ',XTHRCOM,
     &   '         decomposition threshold from input       : ',THRCOM
         NERR = NERR + 1
      END IF

      ERR = ABS(THRDIAG - XTHRDIAG)
      IF (ERR .GT. ERRTOL) THEN
         WRITE(LUPRI,'(A,D16.8,/,A,D16.8)')
     &   'RESTART: init. diag. screening from restart file: ',XTHRDIAG,
     &   '         init. diag. screening from input       : ',THRDIAG
         NERR = NERR + 1
      END IF

      ERR = ABS(DAMP(1) - XDAMP(1))
      IF (ERR .GT. ERRTOL) THEN
         WRITE(LUPRI,'(A,D16.8,/,A,D16.8)')
     &   'RESTART: 1st screening damping from restart file: ',XDAMP(1),
     &   '         1st screening damping from input       : ',DAMP(1)
         NERR = NERR + 1
      END IF

      ERR = ABS(DAMP(2) - XDAMP(2))
      IF (ERR .GT. ERRTOL) THEN
         WRITE(LUPRI,'(A,D16.8,/,A,D16.8)')
     &   'RESTART: 2nd screening damping from restart file: ',XDAMP(2),
     &   '         2nd screening damping from input       : ',DAMP(2)
         NERR = NERR + 1
      END IF

      IF (SCDIAG .NEQV. XSCDIAG) THEN
         IF (XSCDIAG) THEN
            I = 1
         ELSE
            I = 2
         END IF
         IF (SCDIAG) THEN
            J = 1
         ELSE
            J = 2
         END IF
         WRITE(LUPRI,'(A,A,/,A,A)')
     &   'RESTART: diag. screening from restart file: ',SWITCH(I),
     &   '         diag. screening from input       : ',SWITCH(J)
         NERR = NERR + 1
      END IF

      ERR = ABS(THRNEG - XTHRNEG)
      IF (ERR .GT. ERRTOL) THEN
         WRITE(LUPRI,'(A,D16.8,/,A,D16.8)')
     &   'RESTART: neg. diag. threshold from restart file: ',XTHRNEG,
     &   '         neg. diag. threshold from input       : ',THRNEG
         NERR = NERR + 1
      END IF

      ERR = ABS(WARNEG - XWARNEG)
      IF (ERR .GT. ERRTOL) THEN
         WRITE(LUPRI,'(A,D16.8,/,A,D16.8)')
     &   'RESTART: neg. diag. warn thr. from restart file: ',XWARNEG,
     &   '         neg. diag. warn thr. from input       : ',WARNEG
         NERR = NERR + 1
      END IF

      ERR = ABS(TOONEG - XTOONEG)
      IF (ERR .GT. ERRTOL) THEN
         WRITE(LUPRI,'(A,D16.8,/,A,D16.8)')
     &   'RESTART: too neg. diag. thr. from restart file: ',XTOONEG,
     &   '         too neg. diag. thr. from input       : ',TOONEG
         NERR = NERR + 1
      END IF

      ERR = ABS(SPAN - XSPAN)
      IF (ERR .GT. ERRTOL) THEN
         WRITE(LUPRI,'(A,D16.8,/,A,D16.8)')
     &   'RESTART: span factor from restart file: ',XSPAN,
     &   '         span factor from input       : ',SPAN
         NERR = NERR + 1
      END IF

      END
