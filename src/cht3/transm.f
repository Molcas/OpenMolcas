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
      SUBROUTINE TRANSM(V1,V2,ID1,ID2)
      implicit none
      integer ID1,ID2, I
      REAL*8 V1,V2,ONE
      PARAMETER (ONE=1.D0)
      DIMENSION V1(ID1,ID2),V2(ID2,ID1)
C usual transposition
      IF(ID1.EQ.0.OR.ID2.EQ.0)RETURN
      DO  I=1,ID1
         CALL DCOPY_(ID2,V1(I,1),ID1,V2(1,I),1)
      ENDDO
      RETURN
C
      ENTRY TRANSM_A(V1,V2,ID1,ID2)
C transposition of the matrix + add to the target matrix
      IF(ID1.EQ.0.OR.ID2.EQ.0)RETURN
      DO  I=1,ID1
         CALL DAXPY_(ID2,ONE,V1(I,1),ID1,V2(1,I),1)
      ENDDO
      RETURN
C
      ENTRY TRANSM_D(V1,V2,ID1,ID2)
C transposition with splitting to "double spacing"
      IF(ID1.EQ.0.OR.ID2.EQ.0)RETURN
      DO  I=1,ID1
         CALL DCOPY_(ID2/2,V1(I,1),ID1,V2(1,I),2)
      ENDDO
      RETURN
      END
