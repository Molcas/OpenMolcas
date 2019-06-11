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
      Subroutine dirvect( P1, R1, P2, R2,  vec, dist )
c this Subroutine computes the directional vector of the origins of two points P1 and P2
      Implicit None
      Integer, parameter         :: wp=SELECTED_REAL_KIND(p=15,r=307)
      Real(kind=wp), intent(in)  :: P1(3) ! coords of the first point
      Real(kind=wp), intent(in)  :: P2(3) ! coords of the second point
!     rot. matrix of the  first point to the general coordinate system
      Real(kind=wp), intent(in)  :: R1(3,3)
!     rot. matrix of the second point to the general coordinate system
      Real(kind=wp), intent(in)  :: R2(3,3)
      Real(kind=wp), intent(out) :: vec(3)
      Real(kind=wp), intent(out) :: dist
c local variables
      Integer       :: i, j
      Real(kind=wp) :: C1(3), C2(3)
      Real(kind=wp) :: distance
      External      :: distance

      vec=0.0_wp
      dist=0.0_wp
      C1=0.0_wp
      C2=0.0_wp

      Do i=1,3
         Do j=1,3
            C1(i) = C1(i) + P1(j)*R1(i,j)
            C2(i) = C2(i) + P2(j)*R2(i,j)
         End Do
      End Do
      dist=distance(3, C1, C2 )
      Do i=1,3
        vec(i)=( C1(i)-C2(i) ) / dist
      End Do
      Return
      End

