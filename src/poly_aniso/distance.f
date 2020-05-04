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
      Real*8 function distance(N,C1,C2)
      Implicit None
      Integer, parameter        :: wp=SELECTED_REAL_KIND(p=15,r=307)
      Integer, intent(in)       :: N
      Real(kind=8), intent(in) :: C1(N), C2(N)
      ! local variables
      Integer       :: i
      Real(kind=8) :: X, R
      distance=0.0_wp
      X=0.0_wp
      Do i=1,N
        R=0.0_wp
        R=C1(i)-C2(i)
        X=X+R*R
      End Do
      distance=sqrt(X)
      Return
      End
