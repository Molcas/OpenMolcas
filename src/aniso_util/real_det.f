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
      Real*8 Function FindDetR(matrix, n)
c     Function to find the determinant of a square matrix
c
c     Description: The Subroutine is based on two key points:
c
c     1]  A determinant is unaltered when row operations are performed:
c         Hence, using this principle, row operations (column operations
c         would work as well) are used to convert the matrix the matrix
c         into upper triangular form
c     2]  The determinant of a triangular matrix is obtained by finding
c         the product of the diagonal elements
      Implicit None
#include "stdalloc.fh"
      Integer, parameter         :: wp=SELECTED_REAL_KIND(p=15,r=307)
c Calling parameters
      Integer, intent(in)        :: N
      Real(kind=8)              :: matrix(N,N)
c local variables:
      Real(kind=8), allocatable :: w(:), z(:,:)
      Integer                    :: i, info

      info=0
      FindDetR=0.0_wp
      Call mma_allocate(w,n,'eigenvalues')
      Call mma_allocate(z,n,n,'engenvectors')
      Call dcopy_(  n,[0.0_wp],0,w,1)
      Call dcopy_(n*n,[0.0_wp],0,z,1)
      ! diagonalize the matrix:
      Call diag_r2(matrix,n,info,w,z)
      If (info.ne.0) then
         Write(6,*) 'inside FindDetR. diagonalization failed. Info =',
     &               info
         Return
      End If
      ! Evaluate determinant by finding product of diagonal elements
      FindDetR=1.0_wp
      Do i=1, N
         FindDetR = FindDetR * w(i)
      End Do
      Call mma_deallocate(w)
      Call mma_deallocate(z)
      Return
      End
