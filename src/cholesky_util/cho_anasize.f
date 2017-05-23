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
      SUBROUTINE CHO_ANASIZE(VEC,LVEC,BIN,LBIN,LUPRI)
C
C     Purpose: analyse vector (histogram).
C
#include "implicit.fh"
      DIMENSION VEC(LVEC), BIN(LBIN)

      PARAMETER (ZERO = 0.0D0)

      PARAMETER (MBIN = 20)
      INTEGER   ICOUNT(MBIN)
      LOGICAL   FOUND

C     Return if nothing to do.
C     ------------------------

      IF ((LVEC.LT.1) .OR. (LBIN.LT.1)) RETURN

C     Ensure that BIN is in descending order.
C     ---------------------------------------

      IJOB = -1
      CALL CHO_ORDER(BIN,LBIN,IJOB)

C     Test that BIN is positive.
C     --------------------------

      IF (BIN(1) .LE. ZERO) RETURN

C     Analysis.
C     ---------

      NBIN = MIN(LBIN,MBIN)
      CALL CHO_IZERO(ICOUNT,NBIN)
      NLOW = 0
      NZER = 0
      NNEG = 0
      XNEG = ZERO

      DO I = 1,LVEC

         TEST = VEC(I)

         IF (TEST .LT. ZERO) THEN
            NNEG = NNEG + 1
            XNEG = MIN(XNEG,TEST)
         ELSE IF (TEST .EQ. ZERO) THEN
            NZER = NZER + 1
         END IF

         IBIN  = 0
         FOUND = .FALSE.
         DO WHILE ((.NOT.FOUND) .AND. (IBIN.LT.NBIN))
            IBIN = IBIN + 1
            IF (TEST .GE. BIN(IBIN)) THEN
               ICOUNT(IBIN) = ICOUNT(IBIN) + 1
               FOUND = .TRUE.
            END IF
         END DO
         IF (.NOT.FOUND) NLOW = NLOW + 1

      END DO

C     Print.
C     ------

      TOPCT = 1.0D2/DBLE(LVEC)

      JCOUNT = ICOUNT(1)
      WRITE(LUPRI,'(/,1X,A,11X,D11.4,A,I12,1X,F7.2,A,3X,A,F7.2,A)')
     & 'Larger than ',BIN(1),':',ICOUNT(1),
     & DBLE(ICOUNT(1))*TOPCT,'%','Accumulated: ',DBLE(JCOUNT)*TOPCT,'%'
      DO IBIN = 2,NBIN
         JCOUNT = JCOUNT + ICOUNT(IBIN)
         WRITE(LUPRI,'(1X,A,D11.4,A,D11.4,A,I12,1X,F7.2,A,3X,A,F7.2,A)')
     &   'Between ',BIN(IBIN-1),' and ',BIN(IBIN),':',ICOUNT(IBIN),
     &   DBLE(ICOUNT(IBIN))*TOPCT,'%',
     &   'Accumulated: ',DBLE(JCOUNT)*TOPCT,'%'
      END DO
      JCOUNT = JCOUNT + NLOW
      WRITE(LUPRI,'(1X,A,10X,D11.4,A,I12,1X,F7.2,A,3X,A,F7.2,A)')
     & 'Smaller than ',BIN(NBIN),':',NLOW,
     & DBLE(NLOW)*TOPCT,'%','Accumulated: ',DBLE(JCOUNT)*TOPCT,'%'

      WRITE(LUPRI,'(/,1X,A,I12,1X,F7.2,A)')
     & 'Number of elements exactly 0.0D0 :',NZER,
     & DBLE(NZER)*TOPCT,'%'
      WRITE(LUPRI,'(1X,A,I12,1X,F7.2,A)')
     & 'Number of negative elements      :',NNEG,
     & DBLE(NNEG)*TOPCT,'%'
      IF (NNEG .GT. 0) THEN
         WRITE(LUPRI,'(1X,A,D12.4)')
     & ' - numerically largest           :',XNEG
      END IF

      END
