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
      SUBROUTINE ZBASE(NVEC,IVEC,NCLASS)
*
*  Some class division exists with NVEC(ICLASS) members in
*  class ICLASS.
*
*  Construct array IVEC(ICLASS) giving first element of
*  class ICLASS in full adressing
*
      IMPLICIT REAL*8          (A-H,O-Z)
      DIMENSION NVEC(*),IVEC(*)
*
       DO ICLASS = 1,NCLASS
          If (ICLASS.EQ.1) Then
             IVEC(1) = 1
          Else
             IVEC(ICLASS) = IVEC(ICLASS-1)+NVEC(ICLASS-1)
          End If
       End Do
*
      NTEST = 0
      IF( NTEST .NE. 0 ) THEN
        WRITE(6,'(A)') '  ZBASE : NVEC and IVEC '
        WRITE(6,'(A)') '  ===================== '
        CALL IWRTMA(NVEC,1,NCLASS,1,NCLASS)
        CALL IWRTMA(IVEC,1,NCLASS,1,NCLASS)
      END IF
*
      RETURN
      END
