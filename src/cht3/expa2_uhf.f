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
      SUBROUTINE EXPA2_UHF(ARR1,IDM,LI,NSP,ARR2)
      implicit none
      REAL*8 ARR1,ARR2
      integer IDM,LI,NSP, IJ, I,J
C
C THIS SUBROUTINE EXPANDS THE SECOND INDEX OF A MATRIX ARR1
C
      DIMENSION ARR1(IDM,*),ARR2(IDM,LI,*)
      IJ=0
      CALL ZEROMA(ARR2(1,1,1),1,IDM)
      DO  I=2,LI
         DO  J=1,I-1
            IJ=IJ+1
            CALL DCOPY_(IDM,ARR1(1,IJ),1,ARR2(1,I,J),1)
            CALL DCOPY_(IDM,ARR1(1,IJ),1,ARR2(1,J,I),1)
         ENDDO
         CALL ZEROMA(ARR2(1,I,I),1,IDM)
      ENDDO
      IF(NSP.LT.0)THEN
         DO I=1,LI
            CALL VNEG(ARR2(1,1,I),1,ARR2(1,1,I),1,IDM*I)
         enddo
      ENDIF
      RETURN
      END
