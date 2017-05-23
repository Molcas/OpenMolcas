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
      Subroutine SetPrLev_m(LF_IN,IPRGLB_IN,IPRLOC_IN)
      Implicit Real*8 (A-H,O-Z)
#include "warnings.fh"
#include "output_ras.fh"
      Parameter (ROUTINE='SETPRLEV')
      Dimension IPRLOC_IN(7)
*
      Logical REDUCE_PRT
      External REDUCE_PRT
      Intrinsic MAX
      External GETENVF

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

!      IPRLOC(1)=DEBUG !AMS

      IF (IPRLOC(1).GE.DEBUG) THEN
       WRITE(6,*)' SetPrLev: Print levels have been set to'
       write(6,*)'  Global print level IPRGLB=',IPRGLB
       write(6,*)'  Individual sections print levels, IPRLOC:'
       write(6,'(1x,7I5)')(IPRLOC(I),I=1,7)
      END IF

      RETURN
c Avoid unused argument warnings
      IF (.FALSE.) CALL Unused_integer(LF_IN)
      END
