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
      Subroutine SetPrLev(LF_IN,IPRGLB_IN,IPRLOC_IN)
      use printlevel, only: USUAL,DEBUG,SILENT
      use output_ras, only: IPRLOC,IPRGLB
      Implicit None
#include "warnings.h"
      Integer LF_IN,IPRGLB_IN,IPRLOC_IN(7)
*
      Logical, External :: REDUCE_PRT
      Intrinsic MAX
      External GETENVF
      Integer I

* The local print levels are the maximum of the requested global and
* local ones, except that if any of IPRGLB or IPRLOC(I) is zero
*  (meaning silence), then IPRLOC(I) is set to zero.
      IPRGLB=IPRGLB_IN
      IF (IPRGLB_IN.EQ.0) THEN
       DO I=1,7
        IPRLOC(I)=0
       END DO
      ELSE
       DO I=1,7
        IPRLOC(I)=0
        IF(IPRLOC_IN(I).GT.0) IPRLOC(I)=MAX(IPRGLB_IN,IPRLOC_IN(I))
       END DO
      END IF
* If inside an optimization loop, set down the print
* level unless we *really* want a lot of output.
      IF (REDUCE_PRT()) THEN
       IPRGLB=MAX(IPRGLB-USUAL,SILENT)
       DO I=1,7
        IPRLOC(I)=MAX(IPRLOC(I)-USUAL,SILENT)
       END DO
      END IF

      IF (IPRLOC(1).GE.DEBUG) THEN
       WRITE(6,*)' SetPrLev: Print levels have been set to'
       write(6,*)'  Global print level IPRGLB=',IPRGLB
       write(6,*)'  Individual sections print levels, IPRLOC:'
       write(6,'(1x,7I5)')(IPRLOC(I),I=1,7)
      END IF

c Avoid unused argument warnings
#ifdef _WARNING_WORKAROUND_
      IF (.FALSE.) CALL Unused_integer(LF_IN)
#endif
      END Subroutine SetPrLev
