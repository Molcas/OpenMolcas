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
      SUBROUTINE CHO_MCA_DIAGINT(ISHLA,ISHLB,SCR,LSCR)
C
C     Purpose: call Seward to calculate diagonal shell (AB|AB).
C
      IMPLICIT REAL*8 (A-H,O-Z)
      EXTERNAL  Integral_WrOut_Cho_diag
      DIMENSION SCR(LSCR)
#include "itmax.fh"
#include "cholesky.fh"
#if defined (_DEBUGPRINT_)
      CHARACTER*15 SECNAM
      PARAMETER (SECNAM = 'CHO_MCA_DIAGINT')
#endif

      CALL CHO_DZERO(SCR,LSCR)

#if defined (_DEBUGPRINT_)
      CALL CHO_PRESCR(CUTINT1,THRINT1)
#endif

      CALL EVAL_IJKL(ISHLA,ISHLB,ISHLA,ISHLB,SCR,LSCR,
     &               Integral_WrOut_Cho_diag)

#if defined (_DEBUGPRINT_)
      CALL CHO_PRESCR(CUTINT2,THRINT2)
      IF (CUTINT2.NE.CUTINT1 .OR. THRINT2.NE.THRINT1) THEN
         WRITE(LUPRI,*) SECNAM,': CutInt before Eval_Ints_: ',CUTINT1
         WRITE(LUPRI,*) SECNAM,': CutInt after  Eval_Ints_: ',CUTINT2
         WRITE(LUPRI,*) SECNAM,': ThrInt before Eval_Ints_: ',THRINT1
         WRITE(LUPRI,*) SECNAM,': ThrInt after  Eval_Ints_: ',THRINT2
         CALL CHO_QUIT('Integral prescreening error detected in '
     &                 //SECNAM,102)
      END IF
#endif

      END
