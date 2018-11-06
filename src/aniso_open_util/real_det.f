* $ this file belongs to the Molcas repository $
      Real*8 Function FindDetR(matrix, n)
c
c     Function to find the determinant of a square matrix
c
c     Author : Louisda16th a.k.a Ashwith J. Rego
c
c     Description: The Subroutine is based on two key points:
c
c     1]  A determinant is unaltered when row operations are performed:
c         Hence, using this principle, row operations (column operations
c         would work as well) are used to convert the matrix the matrix
c         into upper triangular form
c     2]  The determinant of a triangular matrix is obtained by finding
c         the product of the diagonal elements
c
      Implicit None
      Integer, parameter         :: wp=SELECTED_REAL_KIND(p=15,r=307)
c Calling parameters
      Integer, intent(in)        :: N
      Real(kind=wp)              :: matrix(N,N)
c local variables:
      Real(kind=wp)              :: m, temp
      Integer                    :: i, j, k, l
      Logical                    :: DetExists

      DetExists = .TRUE.
      l = 1
      temp=0
C  Convert to upper triangular form
      Do k = 1, N-1
         If (matrix(k,k) .eq. 0.0_wp) Then
            DetExists = .FALSE.
            Do i = k+1, N
               If (matrix(i,k) .ne. 0.0_wp) Then
                  Do j = 1, N
                           temp  =  matrix(i,j)
                     matrix(i,j) =  matrix(k,j)
                     matrix(k,j) =  temp
                  End Do
                  DetExists = .TRUE.
               End If
            End Do
            If (DetExists .EQV. .FALSE.) Then
               FindDetR = 0
               Return
            End If
         End If
         Do j = k+1, N
            m = matrix(j,k)/matrix(k,k)
            Do i = k+1, N
               matrix(j,i) = matrix(j,i) - m*matrix(k,i)
            End Do
         End Do
      End Do ! k
c  Evaluate determinant by finding product of diagonal elements
      FindDetR = l
      Do i=1, N
         FindDetR = FindDetR * matrix(i,i)
      End Do

      Return
      End
