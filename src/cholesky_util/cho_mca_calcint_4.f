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
      SUBROUTINE CHO_MCA_CALCINT_4(XINT,LINT,ISHLCD,ISHLAB)
C
C     Purpose: calculate qualified integral columns from
C              shell quadruple (ISHLC ISHLD|ISHLA ISHLB).
C
C     Version 4: avoid storage of full shell quadruple in interface to
C                seward; get qualified directly as in Version 2 and 3!
C                Changes from Version 3:
C                - only one shell quadruple is computed (not an entire
C                  column).
C
#include "implicit.fh"
      DIMENSION XINT(LINT)
#include "cholesky.fh"
#include "choprint.fh"
#include "choptr.fh"
#include "WrkSpc.fh"

      CHARACTER*17 SECNAM
      PARAMETER (SECNAM = 'CHO_MCA_CALCINT_4')

      INTEGER  CHO_ISUMELM
      EXTERNAL CHO_ISUMELM

      INTEGER NAB(8)

      LOGICAL   LOCDBG
      PARAMETER (LOCDBG = .FALSE.)
      PARAMETER (INFINT = INF_INT, INFIN2 = INF_IN2)

      ISP2F(I)=IWORK(ip_iSP2F-1+I)

C     Set mapping from shell pair AB to qualified columns.
C     ----------------------------------------------------

      IRC  = 0
      ILOC = 2
      CALL CHO_SETSHP2Q_2(IRC,ILOC,ISHLAB,NAB)
      IF (IRC .NE. 0) THEN
         WRITE(LUPRI,*) SECNAM,': CHO_SETSHP2Q_2 returned ',IRC
         CALL CHO_QUIT('Error termination in '//SECNAM,IRC)
      END IF

C     Print.
C     ------

      IF (IPRINT .GE. INF_IN2) THEN
         CALL CHO_INVPCK(ISP2F(ISHLAB),ISHLA,ISHLB,.TRUE.)
         CALL CHO_INVPCK(ISP2F(ISHLCD),ISHLC,ISHLD,.TRUE.)
         NCOLAB = CHO_ISUMELM(NAB,NSYM)
         WRITE(LUPRI,'(/,A,I5,1X,I5,A,I5,1X,I5,A,I9,A)')
     &   'Calculating shell quadruple (',ISHLC,ISHLD,'|',ISHLA,ISHLB,
     &   '):',NCOLAB,' columns have been qualified'
         WRITE(LUPRI,'(89A)') ('=',i=1,89)
      END IF

C     Set mapping from shell pair CD to reduced set.
C     ----------------------------------------------

      IRC  = 0
      ILOC = 2
      CALL CHO_SETSHP2RS_2(IRC,ILOC,ISHLCD,NAB)
      IF (IRC .NE. 0) THEN
         WRITE(LUPRI,*) SECNAM,': CHO_SETSHP2RS_2 returned ',IRC
         CALL CHO_QUIT('Error termination in '//SECNAM,IRC)
      END IF

C     Calculate integrals.
C     --------------------

      CALL CHO_TIMER(C1,W1)
      CALL CHO_MCA_INT_1(ISHLCD,ISHLAB,
     &                   XINT,LINT,
     &                   LOCDBG.OR.(IPRINT.GE.100))
      CALL CHO_TIMER(C2,W2)
      TINTEG(1,1) = TINTEG(1,1) + C2 - C1
      TINTEG(2,1) = TINTEG(2,1) + W2 - W1

      END
