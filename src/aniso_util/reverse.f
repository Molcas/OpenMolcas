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
      Subroutine REVERSE(A_dir,A_inv,det)
C
C THIS ROUTINE CALCULATES THE INVERSE OF A SQUARE 3x3 MATRIX, AND ITS DETERMINANT.
C

      Implicit None
      Integer, parameter        :: wp=SELECTED_REAL_KIND(p=15,r=307)
      Real(kind=8)  :: A_dir(3,3)
      Real(kind=8)  :: A_inv(3,3)
      Real(kind=8)  :: A(3,3)
      Real(kind=8)  :: B(3,3)
      Real(kind=8)  :: det
      Real(kind=8)  :: FindDetR
      External       :: FindDetR

      det=0.0_wp
      Call dcopy_(3*3,[0.0_wp],0,A,1)
      Call dcopy_(3*3,[0.0_wp],0,B,1)
      Call dcopy_(3*3,[0.0_wp],0,A_inv,1)
      Call dcopy_(3*3,A_dir,1,A,1)

      det=FindDetR(A,3)

      B(1,1)= A(2,2)*A(3,3)-A(3,2)*A(2,3)
      B(1,2)=-A(1,2)*A(3,3)+A(3,2)*A(1,3)
      B(1,3)= A(1,2)*A(2,3)-A(2,2)*A(1,3)
      B(2,1)=-A(2,1)*A(3,3)+A(3,1)*A(2,3)
      B(2,2)= A(1,1)*A(3,3)-A(3,1)*A(1,3)
      B(2,3)=-A(1,1)*A(2,3)+A(2,1)*A(1,3)
      B(3,1)= A(2,1)*A(3,2)-A(3,1)*A(2,2)
      B(3,2)=-A(1,1)*A(3,2)+A(3,1)*A(1,2)
      B(3,3)= A(1,1)*A(2,2)-A(2,1)*A(1,2)

      Call dscal_(3*3, 1.0_wp/det, B, 1 )
      Call dcopy_(3*3, B,1,A_inv,1)

      Return
      End
