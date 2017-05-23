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
      SUBROUTINE DGETMI (A,ldA,N)
C
C     TRANSPOSE A SQUARE MATRIX (IN-PLACE)
C
      Real*8 A(ldA,*), Temp
C
      If ( N.le.0 ) then
         Write (6,*)
         Write (6,*)
     &   '  *** Error in subroutine DGETMI ***'
         Write (6,*)
     &   '  Invalid dimension of matrix A :'
         Write (6,*)
     &   '  The number of rows/columns, N, must be larger than zero'
         Write (6,*)
      End If
      If ( (ldA.le.0) .or. (ldA.lt.N) ) then
         Write (6,*)
         Write (6,*)
     &   '  *** Error in subroutine DGETMI ***'
         Write (6,*)
     &   '  Invalid leading dimension of matrix A :'
         Write (6,*)
     &   '  ldA must be larger than 0 and larger than N'
         Write (6,*)
      End If
C
      Do i=1,N
         Do j=1,i-1
            Temp=A(j,i)
            A(j,i)=A(i,j)
            A(i,j)=Temp
         End Do
      End Do
C
      RETURN
      END
