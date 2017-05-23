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
      SUBROUTINE EXPA1_UHF(ARR1,IDM,LI,NSP,ARR2)
      implicit none
      REAL*8 ARR1,ARR2
      integer IDM,LI,NSP, IJ, I,J,K
C
C THIS SUBROUTINE EXPANDS THE FIRST INDEX OF A MATRIX ARR1
C
C    LI*(LI+1)/2
      DIMENSION ARR1(*),ARR2(LI,LI,*)
      IF(NSP.GT.0)THEN
         IJ=1
         DO K=1,IDM
            DO I=1,LI
               CALL DCOPY_(I,ARR1(IJ),1,ARR2(I,1,K),LI)
               CALL DCOPY_(I,ARR1(IJ),1,ARR2(1,I,K),1)
               IJ=IJ+I
            enddo
         enddo
      ELSE
         IJ=1
         DO K=1,IDM
            ARR2(1,1,K)=0.D0
            DO I=2,LI
               ARR2(I,I,K)=0.D0
               CALL DCOPY_(I-1,ARR1(IJ),1,ARR2(I,1,K),LI)
               DO J=1,I-1
                  ARR2(J,I,K)=-ARR1(IJ)
                  IJ=IJ+1
               ENDDO
            enddo
         enddo
      ENDIF
      RETURN
      END
