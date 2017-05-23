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
      SUBROUTINE DGETMO (A,ldA,M,N,B,ldB)
C
C     TRANSPOSE A REGULAR MATRIX (OUT-OF-PLACE)
C
      REAL*8 A(ldA,*),B(ldB,*)
C
      If ( M.le.0 ) then
         Write (6,*)
         Write (6,*)
     &   '  *** Error in subroutine DGETMO ***'
         Write (6,*)
     &   '  Invalid dimension of matrix A :'
         Write (6,*)
     &   '  The number of columns, M, must be larger than zero'
         Write (6,*)
      End If
      If ( N.le.0 ) then
         Write (6,*)
         Write (6,*)
     &   '  *** Error in subroutine DGETMO ***'
         Write (6,*)
     &   '  Invalid leading dimension of matrix B :'
         Write (6,*)
     &   '  The number of rows, N, must be larger than zero'
         Write (6,*)
      End If
      If ( (ldA.le.0) .or. (ldA.lt.M) ) then
         Write (6,*)
         Write (6,*)
     &   '  *** Error in subroutine DGETMO ***'
         Write (6,*)
     &   '  Invalid leading dimension of matrix A :'
         Write (6,*)
     &   '  ldA must be larger than 0 and larger than M'
         Write (6,*)
      End If
      If ( (ldB.le.0) .or. (ldB.lt.N) ) then
         Write (6,*)
         Write (6,*)
     &   '  *** Error in subroutine DGETMO ***'
         Write (6,*)
     &   '  Invalid leading dimension of matrix B :'
         Write (6,*)
     &   '  ldB must be larger than 0 and larger than N'
         Write (6,*)
      End If
      INC=8
      Do j=1,M,INC
         jj=Min(M-j+1,INC)
         Go To (1,2,3,4,5,6,7,8) jj
         Write (6,*) 'Error in DGETMO!'
 1       Continue
            Do i=1,N
               B(i,j)=A(j,i)
            End Do
         Go To 99
 2       Continue
            Do i=1,N
               B(i,j  )=A(j  ,i)
               B(i,j+1)=A(j+1,i)
            End Do
         Go To 99
 3       Continue
            Do i=1,N
               B(i,j  )=A(j  ,i)
               B(i,j+1)=A(j+1,i)
               B(i,j+2)=A(j+2,i)
            End Do
         Go To 99
 4       Continue
            Do i=1,N
               B(i,j  )=A(j  ,i)
               B(i,j+1)=A(j+1,i)
               B(i,j+2)=A(j+2,i)
               B(i,j+3)=A(j+3,i)
            End Do
         Go To 99
 5       Continue
            Do i=1,N
               B(i,j  )=A(j  ,i)
               B(i,j+1)=A(j+1,i)
               B(i,j+2)=A(j+2,i)
               B(i,j+3)=A(j+3,i)
               B(i,j+4)=A(j+4,i)
            End Do
         Go To 99
 6       Continue
            Do i=1,N
               B(i,j  )=A(j  ,i)
               B(i,j+1)=A(j+1,i)
               B(i,j+2)=A(j+2,i)
               B(i,j+3)=A(j+3,i)
               B(i,j+4)=A(j+4,i)
               B(i,j+5)=A(j+5,i)
            End Do
         Go To 99
 7       Continue
            Do i=1,N
               B(i,j  )=A(j  ,i)
               B(i,j+1)=A(j+1,i)
               B(i,j+2)=A(j+2,i)
               B(i,j+3)=A(j+3,i)
               B(i,j+4)=A(j+4,i)
               B(i,j+5)=A(j+5,i)
               B(i,j+6)=A(j+6,i)
            End Do
         Go To 99
 8       Continue
            Do i=1,N
               B(i,j  )=A(j  ,i)
               B(i,j+1)=A(j+1,i)
               B(i,j+2)=A(j+2,i)
               B(i,j+3)=A(j+3,i)
               B(i,j+4)=A(j+4,i)
               B(i,j+5)=A(j+5,i)
               B(i,j+6)=A(j+6,i)
               B(i,j+7)=A(j+7,i)
            End Do
 99      Continue
      End Do
C
      Return
      End
